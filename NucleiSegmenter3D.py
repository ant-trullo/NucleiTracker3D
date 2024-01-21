"""This function segments nuclei in 3D.

It uses standard image analysis techiques for the segmentation, then thanks
to a classifier previously trained it recognizes nuclei badly segmented and
performes the last modifications.
Input is raw data and training.
"""


import itertools
import numpy as np
from skimage.measure import label, regionprops_table
from skimage.filters import gaussian, threshold_otsu
from skimage.segmentation import watershed, expand_labels
from skimage.feature import peak_local_max
from scipy.ndimage import distance_transform_edt
import joblib
# import pyqtgraph as pg

import ServiceWidgets


def all_possible_subsets(s):
    """This function gives all the possible subset (without repetitions) of the first s numbers."""
    sub_sets  =  []
    for k in range(1, s + 1):
        sub_sets  +=  list(itertools.combinations(np.arange(s), k))
    return sub_sets


class NucleiSegmenter3D:
    """Main class, does all the job."""
    def __init__(self, nucs):

        model     =  joblib.load('finalized_model.sav')                                                                # import the pretrained classifier
        nucs_sgm  =  np.zeros_like(nucs)
        tlen      =  nucs.shape[0]
        pbar      =  ServiceWidgets.ProgressBar(total1=tlen + 1)
        pbar.update_progressbar1(1)
        pbar.show()

        for tt in range(tlen):
            pbar.update_progressbar1(tt + 1)
            aa_g0       =  gaussian(nucs[tt], 2)                                                                            # smoothing with a gaussian filter
            val1        =  threshold_otsu(aa_g0)                                                                            # threshold for detetction
            aze1        =  aa_g0 > val1
            nuc_lbls    =  label(aze1).astype(np.uint16)                                                                    # label connected components
            distance    =  distance_transform_edt(nuc_lbls)                                                                 # distance matrix for watershed
            local_maxi  =  peak_local_max(distance, footprint=np.ones((5, 11, 11)), labels=label(nuc_lbls))                 # find peaks of distance matrix
            markers     =  np.zeros_like(nuc_lbls)
            for k in local_maxi:
                markers[k[0], k[1], k[2]]  =  1                                                                             # peaks in a 3D matrix
            markers     =  label(markers)                                                                                   # label markers
            nuc_fin     =  watershed(-distance, markers, mask=np.sign(nuc_lbls))                                            # watershed

            rgp_tot  =  regionprops_table(nuc_fin, properties=["label", "solidity", "coords"])                              # regionprops of the segmtend nuclei
            uus      =  np.where(np.isinf(rgp_tot["solidity"]))[0]                                                          # identify spots with infinitive solidity (clear errors in segmentation, non-3D objects)
            for uu in uus:
                bff_mtx  =  (expand_labels((nuc_fin == rgp_tot["label"][uu])) - (nuc_fin == rgp_tot["label"][uu]) * 1) * nuc_fin    # expand label of problematic detection
                bff_idx  =  np.unique(bff_mtx[bff_mtx != 0])                                                                        # find the label of touching nuclei
                if bff_idx.size > 0:
                    nuc_fin[rgp_tot["coords"][uu][:, 0], rgp_tot["coords"][uu][:, 1], rgp_tot["coords"][uu][:, 2]]  =  bff_idx[0]  # change the label of the problematic component into the label of touching object (if more than one touching objects is present, just connect with the first -solidity infinitive objects are super small)

            rgp_new              =  regionprops_table(nuc_fin, properties=["label", "coords"])                              # re-run regionprops on the corrected matrix
            msk_brdr             =  np.zeros_like(nuc_fin)                                                                  # mask to identify all the nuclei touching the border
            msk_brdr[:, :2, :]   =  1                                                                                       # put a border on the mask (2 pixels thick)
            msk_brdr[:, -2:, :]  =  1
            msk_brdr[:, :, :2]   =  1
            msk_brdr[:, :, -2:]  =  1
            idx_brd              =  msk_brdr * nuc_fin                                                                      # multiply mask and segmented nuclei retrive the label of touching nuclei
            idx_brd              =  np.unique(idx_brd[idx_brd != 0])                                                        # retrive the labels of touching nuclei
            for ee in idx_brd:
                ee_idx  =  np.where(rgp_new["label"] == ee)[0][0]                                                           # remove nuclei touching the border
                nuc_fin[rgp_new["coords"][ee_idx][:, 0], rgp_new["coords"][ee_idx][:, 1], rgp_new["coords"][ee_idx][:, 2]]  =  0            # use coords since it is faster

            nuc_fin   =  label(nuc_fin)                                                                                      # re-label new segmented nuclei matrix
            # nuc_fin2  =  np.copy(nuc_fin)
            feat_bff  =  np.squeeze(np.asarray(list(regionprops_table(nuc_fin, properties=["label", "area", "area_bbox", "area_convex", "area_filled", "axis_major_length", "axis_minor_length",        # regionprops: takes directly the array of the values
                                                                                           "equivalent_diameter_area", "extent", "feret_diameter_max", "inertia_tensor_eigvals", "solidity", "inertia_tensor"]).values()))).T
            features   =  feat_bff[:, 1:]
            iso_lbls   =  feat_bff[:, 0].astype(np.int64)
            preds      =  model.predict(features)                                                                           # run the classifier (trained before) on the features matrix
            idxs_over  =  np.where(preds == 1)[0]                                                                           # identify the over segmented nuclei (labeled 1)
            over_mtx   =  np.zeros_like(nuc_fin)                                                                            # initialize the matrix of the over segmented nuclei
            for mm in idxs_over:
                over_mtx  +=  iso_lbls[mm] * (nuc_fin == iso_lbls[mm])                                                      # add the nuclei to the over segmented matrix
            # over_mtx2  =  np.copy(over_mtx)
            for idx_over in idxs_over:                                                                                      # for all the over segmentes obejects
                # print(idx_over)
                over_bff  =  (over_mtx == iso_lbls[idx_over])                                                               # take the connected component form the over segmented matrix
                over_bff  =  (expand_labels(over_bff) ^ over_bff) * over_mtx                                                # expand the component, take only the expansion and multiply it with the over_mtx
                bff_idxs  =  np.unique(over_bff[over_bff != 0])                                                             # check all the tags in the expansion (tags of the touching nuclei)
                new_nuc   =  over_mtx == iso_lbls[idx_over]                                                                 # initialize the new nucleus matrix: add first the connected component with idx_over
                if bff_idxs.size == 1:                                                                                      # if there is 1 touching objetc  (if there are 0 objects, it means that the connected component is isolated, nothing to do about it, maybe the classifier made a mistake)
                    new_nuc  +=  (over_mtx == bff_idxs[0])                                                                    # add it to the nucleus
                elif bff_idxs.size > 1:                                                                                     # in case of several objects touching the component
                    # break
                    combs          =  all_possible_subsets(bff_idxs.size)                                                   # all the possible combination of the indexes of touching objects
                    new_nucs_test  =  []                                                                                    # it is possible to combine all the touching objects in several ways, so list of possible new nuclei
                    for uu in combs:                                                                                        # for each of the possible combination
                        new_nuc_bff  =  over_mtx == iso_lbls[idx_over]                                                      # start with the connected component idx_over
                        for sub_uu in uu:                                                                                   # for all the indexes
                            new_nuc_bff  +=  (over_mtx == bff_idxs[sub_uu])                                                 # build one of the possible nuclei
                        new_nucs_test.append(new_nuc_bff)                                                                   # add to the list of possible nuclei
                    features_bff  =  np.zeros((len(new_nucs_test), 22))                                                     # initialize features matrix
                    for cnt, pp in enumerate(new_nucs_test):                                                                # for each of the possible nuclei measure features and add it in the matrix as before
                        features_bff[cnt]  =  np.asarray(list(regionprops_table(pp * 1, properties=["area", "area_bbox", "area_convex", "area_filled", "axis_major_length", "axis_minor_length", "equivalent_diameter_area",
                                                                                                    "extent", "feret_diameter_max", "inertia_tensor_eigvals", "solidity", "inertia_tensor"]).values())).T
                    preds_bff  =  model.predict(features_bff)                                                               # use the trained model to make predictions on the possible nuclei

                    if 0 in preds_bff:                                                                                      # check if there is at least 1 possbile nucleus that the trained model recognizes as a well segmented one.
                        qq         =  np.where(preds_bff == 0)[0][-1]                                                       # if there is, take it (in case they are more than one, take the one with less connected components involved)
                        new_nuc    =  new_nucs_test[qq]                                                                     # new nuc to add is the one selected before

                nuc_fin   *=  (1 - new_nuc)                                                                                 # remove all the involved connected components froml the original matrix
                over_mtx  *=  (1 - new_nuc)                                                                                 # remove the selected objects from the matrix of the over segmented
                nuc_fin   +=  new_nuc * iso_lbls[idx_over]                                                                  # add the reconstructed new nucleus to the final matrix

            nucs_sgm[tt]  =  label(nuc_fin)

        pbar.close()

        self.nucs_sgm  =  nucs_sgm




        # mycmap = np.fromfile("mycmap.bin", "uint16").reshape((10000, 3))  # / 255.0
        # colors4map = []
        # for k in range(mycmap.shape[0]):
        #     colors4map.append(mycmap[k, :])
        # colors4map     =  colors4map
        # colors4map[0] = np.array([0, 0, 0])
        # w = pg.image(nuc_fin, title='Treated')
        # nucs_cmap = pg.ColorMap(np.linspace(0, 1, nuc_fin.max()), color=colors4map)
        # w.setColorMap(nucs_cmap)
        # w2 = pg.image(nuc_fin2, title='Not Treated')
        # nucs_cmap = pg.ColorMap(np.linspace(0, 1, nuc_fin2.max()), color=colors4map)
        # w2.setColorMap(nucs_cmap)
        #
        #
        #
        #
        #
        # under_mtx   =  np.zeros_like(nuc_fin)                                                                           # initialize the matrix of the under segmented matrix
        # idxs_under  =  np.where(preds == 2)[0]                                                                          # identify the under segmented nuclei (labeled 2)
        # for vv in idxs_under:
        #     under_mtx  +=  vv * (nuc_fin == iso_lbls[vv])                                                               # add the nuclei to the under segmented matrix
        #
        # well_mtx   =  np.zeros_like(nuc_fin)                                                                           # initialize the matrix of the under segmented matrix
        # idxs_well  =  np.where(preds == 0)[0]                                                                          # identify the under segmented nuclei (labeled 2)
        # for vv in idxs_well:
        #     well_mtx  +=  vv * (nuc_fin2 == iso_lbls[vv])                                                               # add the nuclei to the under segmented matrix
        #
        # w3 = pg.image(well_mtx, title='Well Done 1st shot')
        # nucs_cmap = pg.ColorMap(np.linspace(0, 1, well_mtx.max()), color=colors4map)
        # w3.setColorMap(nucs_cmap)




        # z_edge, x_edge, y_edge  =  0, 0, 0                                                                              # initialize edges size with 0
        # for k in rgp["coords"]:                                                                                         # for each nucleus, we take its coords and find the edge capable of contain the nucleus along each dimension
        #     z_edge  =  max(z_edge, k[:, 0].max() - k[:, 0].min())                                                       # we can create like this a list of boxes of the same size with the nucleus inside (for classification purpouses)
        #     x_edge  =  max(x_edge, k[:, 1].max() - k[:, 1].min())
        #     y_edge  =  max(y_edge, k[:, 2].max() - k[:, 2].min())
        #
        # z_edge  +=  1                                                                                                   # add 1 to not have the nucs across the border
        # x_edge  +=  1
        # y_edge  +=  1
        #
        # if (z_edge % 2) != 0:                                                                                           # for some reason I need the edge to be even (???????? not a big deal, but...)
        #     z_edge  +=  1
        # if (x_edge % 2) != 0:
        #     x_edge  +=  1
        # if (y_edge % 2) != 0:
        #     y_edge  +=  1




        # mycmap      =  np.fromfile("mycmap.bin", "uint16").reshape((10000, 3))  # / 255.0
        # colors4map  =  []
        # for k in range(mycmap.shape[0]):
        #     colors4map.append(mycmap[k, :])
        # # colors4map     =  colors4map
        # colors4map[0] = np.array([0, 0, 0])
        #
        # w          =  pg.image(nuc_fin)
        # nucs_cmap  =  pg.ColorMap(np.linspace(0, 1, nuc_fin.max()), color=colors4map)
        # w.setColorMap(nucs_cmap)
        #
        # ww          =  pg.image(under_mtx)
        # under_cmap  =  pg.ColorMap(np.linspace(0, 1, under_mtx.max()), color=colors4map)
        # ww.setColorMap(under_cmap)
        #
        # www        =  pg.image(over_mtx)
        # over_cmap  =  pg.ColorMap(np.linspace(0, 1, over_mtx.max()), color=colors4map)
        # www.setColorMap(over_cmap)
