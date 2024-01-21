"""This function tracks 3D nuclei.

Input is the 4D matrix of the segmented nuclei plus a threshold distance.
"""

import numpy as np
from skimage.measure import regionprops_table
from PyQt5 import QtWidgets

import ServiceWidgets


class Nuclei3dTracker:
    """This class performs nuclei tracking in 3D."""
    def __init__(self, nucs_segm, dist_thr, pix_size_XY, pix_size_Z):

        nucs_trck   =  np.zeros_like(nucs_segm)                     # initialize output matrix of tracked nuclei

        tlen       =  nucs_segm.shape[0]                            # number of time frames
        ctrs_mtx   =  np.squeeze(np.asarray(list(regionprops_table(nucs_segm[0], properties=["label", "centroid"]).values())))      # matrix of centroids: label, z, x, y. time is missing
        ctrs_mtx   =  np.append(ctrs_mtx, np.zeros((ctrs_mtx.shape[1]))[None, :], 0)    # 5XN matrix: you have label, z, x, y, time as centroid coordinates for each nucleus. Here a column for time is add. It is all zeros since it is the first time frame

        for tt in range(1, tlen):
            ctrs_bff  =  np.squeeze(np.asarray(list(regionprops_table(nucs_segm[tt], properties=["label", "centroid"]).values())))   # extract label, z, x, y for the other time steps
            ctrs_bff  =  np.append(ctrs_bff, tt * np.ones((ctrs_bff.shape[1]))[None, :], 0)                                          # add time to the matrix (now it is label, z, x, y, time)
            ctrs_mtx  =  np.concatenate((ctrs_mtx, ctrs_bff), axis=1)                                                                # concatenate with hte previous

        pbar      =  ServiceWidgets.ProgressBar(total1=ctrs_mtx.shape[1])
        pbar_ref  =  ctrs_mtx.shape[1]
        pbar.update_progressbar1(0)
        pbar.show()

        new_lbl  =  0                                                               # new label is set to 0 since it is increased by 1
        while ctrs_mtx.shape[1] != 0:                                               # any time some coordinates are taken, they are removed from the matrix. Tracking is done when the ctrs_mtx matrix is empty                       
            pbar.update_progressbar1(pbar_ref - ctrs_mtx.shape[1])
            new_lbl  +=  1                                                          # update label to assign
            ctrs_ref  =  ctrs_mtx[:, 0]                                             # take a new ref ctrs
            ctrs_mtx  =  np.delete(ctrs_mtx, 0, axis=1)                       # delete the ref from the matrix
            nucs_trck[int(ctrs_ref[-1])]  +=  ((nucs_segm[int(ctrs_ref[-1])]  ==  ctrs_ref[0]) * new_lbl).astype(np.uint16)   # add the corresponding nucleus to the new matrix
            for uu in range(int(ctrs_ref[-1]) + 1, tlen):                           # for all the following time frames
                ctrs_ext  =  ctrs_mtx[:, ctrs_mtx[4, :] == uu]                      # isolate the centroids with proper time frame 
                distance  =  ctrs_ext - ctrs_ref[:, None] * np.ones_like(ctrs_ext)  # define a matrix to calculate the distance of the ref from all the centroids in the following frame
                ppp       =  np.sqrt((distance[1, :] * pix_size_Z) ** 2 + (distance[2, :] * pix_size_XY) ** 2 + (distance[3, :] * pix_size_XY) ** 2)    # explicit distance calculation (it is in µm since we multiply for pixel size --xy and z are different)
                # print(ppp.min())
                if ppp.min() < dist_thr:                 # check distance is smaller than distance threshold
                    ll                                =  np.argmin(ppp)                 # coordinate of the min of the distance 
                    nucs_trck[int(ctrs_ext[-1, 0])]  +=  ((nucs_segm[int(ctrs_ext[-1, 0])] == ctrs_ext[0, ll]) * new_lbl).astype(np.uint16)    # add the new nucleus to the final matrix 
                    dd                                =  np.intersect1d(np.where(ctrs_mtx[-1] == ctrs_ext[-1, ll])[0], np.where(ctrs_mtx[0] == ctrs_ext[0, ll])[0])     # identify the ctrs coordinates of the just found nucleus: use the time and the label since for each time frame nuclei labels are unique
                    ctrs_mtx                          =  np.delete(ctrs_mtx, dd, axis=1)            # remove from the total matrix
                else:
                    break                                # if all the distances are higher than the threshold, just go for the next nucleus 

        self.nucs_trck  =  nucs_trck


class Nuclei2dTracker:
    """This class performs nuclei tracking considering only the xy coordinate of the centroids."""
    def __init__(self, nucs_segm, dist_thr, pix_size_XY):

        nucs_trck   =  np.zeros_like(nucs_segm)                     # initialize output matrix of tracked nuclei

        tlen       =  nucs_segm.shape[0]                            # number of time frames
        ctrs_mtx   =  np.squeeze(np.asarray(list(regionprops_table(nucs_segm[0], properties=["label", "centroid"]).values())))      # matrix of centroids: label, z, x, y. time is missing
        ctrs_mtx   =  np.delete(ctrs_mtx, 1, 0)                                                                                     # remove the z coordinate: 3xN matrix label, x, y
        ctrs_mtx   =  np.append(ctrs_mtx, np.zeros((ctrs_mtx.shape[1]))[None, :], 0)    # 4XN matrix: you have label, x, y, time as centroid coordinates for each nucleus. Here a colomn for time is add. It is all zeros since it is the first time frame

        for tt in range(1, tlen):
            ctrs_bff  =  np.squeeze(np.asarray(list(regionprops_table(nucs_segm[tt], properties=["label", "centroid"]).values())))   # extract label, z, x, y for the other time steps
            ctrs_bff  =  np.delete(ctrs_bff, 1, 0)                                                                                   # remove the z coordinate: 3xN matrix label, x, y
            ctrs_bff  =  np.append(ctrs_bff, tt * np.ones((ctrs_bff.shape[1]))[None, :], 0)                                          # add time to the matrix (now it is label, z, x, y, time)
            ctrs_mtx  =  np.concatenate((ctrs_mtx, ctrs_bff), axis=1)                                                                # concatenate with hte previous

        pbar      =  ServiceWidgets.ProgressBar(total1=ctrs_mtx.shape[1])
        pbar_ref  =  ctrs_mtx.shape[1]
        pbar.update_progressbar1(0)
        pbar.show()

        new_lbl  =  0                                                               # new label is set to 0 since it is increased by 1
        while ctrs_mtx.shape[1] != 0:                                               # any time some coordinates are taken, they are removed from the matrix. Tracking is done when the ctrs_mtx matrix is empty
            pbar.update_progressbar1(pbar_ref - ctrs_mtx.shape[1])
            new_lbl  +=  1                                                          # update label to assign
            ctrs_ref  =  ctrs_mtx[:, 0]                                             # take a new ref ctrs
            ctrs_mtx  =  np.delete(ctrs_mtx, 0, axis=1)                             # delete the ref from the matrix
            nucs_trck[int(ctrs_ref[-1])]  +=  ((nucs_segm[int(ctrs_ref[-1])]  ==  ctrs_ref[0]) * new_lbl).astype(np.uint16)   # add the corresponding nucleus to the new matrix
            for uu in range(int(ctrs_ref[-1]) + 1, tlen):                           # for all the following time frames
                ctrs_ext  =  ctrs_mtx[:, ctrs_mtx[3, :] == uu]                      # isolate the centroids with proper time frame
                distance  =  ctrs_ext - ctrs_ref[:, None] * np.ones_like(ctrs_ext)  # define a matrix to calculate the distance of the ref from all the centroids in the following frame
                ppp       =  np.sqrt((distance[1, :] * pix_size_XY) ** 2 + (distance[2, :] * pix_size_XY) ** 2)    # explicit distance calculation (it is in µm since we multiply for pixel size --xy and z are different)
                # print(ppp.min())
                if ppp.min() < dist_thr:                 # check distance is smaller than distance threshold
                    ll                                =  np.argmin(ppp)                 # coordinate of the min of the distance
                    nucs_trck[int(ctrs_ext[-1, 0])]  +=  ((nucs_segm[int(ctrs_ext[-1, 0])] == ctrs_ext[0, ll]) * new_lbl).astype(np.uint16)    # add the new nucleus to the final matrix
                    dd                                =  np.intersect1d(np.where(ctrs_mtx[-1] == ctrs_ext[-1, ll])[0], np.where(ctrs_mtx[0] == ctrs_ext[0, ll])[0])     # identify the ctrs coordinates of the just found nucleus: use the time and the label since for each time frame nuclei labels are unique
                    ctrs_mtx                          =  np.delete(ctrs_mtx, dd, axis=1)            # remove from the total matrix
                else:
                    break                                # if all the distances are higher than the threshold, just go for the next nucleus

        self.nucs_trck  =  nucs_trck


class NucleiOverlTracker:
    """This class performs nuclei tracking projecting the calc of each nucleus in the following frame and check the overlapping one by median."""
    def __init__(self, nucs_segm):

        nucs_trck     =  np.zeros_like(nucs_segm)                                   # initialize output matrix of tracked nuclei
        nucs_trck[0]  =  np.copy(nucs_segm[0])                                      # first frameof tracked can be taken as the first frame of segmented
        tlen          =  nucs_segm.shape[0]                                         # number of time frames

        pbar  =  ServiceWidgets.ProgressBar(total1=tlen)
        pbar.show()
        pbar.update_progressbar1(0)

        for tt in range(1, tlen):                                                       # for all the frames following the first
            pbar.update_progressbar1(tt)
            nucs_bff   =  nucs_trck[tt - 1]                                             # isolate the previous tracked frame matrix
            rgp_bff    =  regionprops_table(nucs_bff, properties=["label", "coords"])   # regionprops for label and coordinates (work with coordinates is faster)
            nucs_tags  =  np.unique(nucs_bff[nucs_bff != 0])                            # array of all the tags in the previous tracked frame
            rgp_next   =  regionprops_table(nucs_segm[tt], properties=["label", "coords"])  # regionprops of the segmented nuclei in the following frames for label and coords
            for nuc_tag in nucs_tags:                                                   # for each nucleus
                jj_bff   =  np.where(rgp_bff["label"] == nuc_tag)[0][0]                 # search the index of the proper label
                nxt_tag  =  nucs_segm[tt, rgp_bff["coords"][jj_bff][:, 0], rgp_bff["coords"][jj_bff][:, 1], rgp_bff["coords"][jj_bff][:, 2]]  # use nucleus coordinates in the following frame to check the tags of the overlapping nuclei
                nxt_tag  =  np.median(nxt_tag[nxt_tag != 0])                            # take the median (the nucleus mostly overlapped)
                if nxt_tag - np.fix(nxt_tag) == 0.0:                                    # eventually, in some special and rare cases, you can have deciaml results; these must be avoided (maybe do something smarter than just remove)
                    jj_next         =  np.where(rgp_next["label"] == nxt_tag)[0][0]     # check the index of the tag in the label dictionary
                    nucs_trck[tt, rgp_next["coords"][jj_next][:, 0], rgp_next["coords"][jj_next][:, 1], rgp_next["coords"][jj_next][:, 2]]  =  nuc_tag   # add the nucleus with the proper tag to the matrix of the tracked nuclei

        pbar.close()
        self.nucs_trck  =  nucs_trck


# class ProgressBar(QtWidgets.QWidget):
#     """Simple progressbar widget."""
#     def __init__(self, parent=None, total1=20):
#         super().__init__(parent)
#         self.name_line1  =  QtWidgets.QLineEdit()
#
#         self.progressbar1  =  QtWidgets.QProgressBar()
#         self.progressbar1.setMinimum(1)
#         self.progressbar1.setMaximum(total1)
#
#         main_layout  =  QtWidgets.QGridLayout()
#         main_layout.addWidget(self.progressbar1, 0, 0)
#
#         self.setLayout(main_layout)
#         self.setWindowTitle("Progress")
#         self.setGeometry(500, 300, 300, 50)
#
#     def update_progressbar1(self, val1):
#         """Update progressbar."""
#         self.progressbar1.setValue(val1)
#         QtWidgets.qApp.processEvents()
