"""This function loads raw data files and gives as output the matrix video of the 4D raw data
and the MIP of them, plus some metadata.

This function assumes that there are not going to be singgle time frame files.
"""


import numpy as np
from aicsimageio import AICSImage
import czifile
from natsort import natsorted

import ServiceWidgets


class LoadRawDataCzi:
    """Only function, does all the job. """
    def __init__(self, fnames, analysis_folder=None):

        fnames           =  natsorted(fnames, key=lambda y: y.lower())                  # natural order for file names
        time_step_value  =  0
        j                =  0
        while time_step_value == 0:                                                     # take into account that the first file can have only 1 time frame
            with czifile.CziFile(str(fnames[j])) as czi:
                for attachment in czi.attachments():
                    if attachment.attachment_entry.name == 'TimeStamps':
                        timestamps = attachment.data()
                        break
                else:
                    raise ValueError('TimeStamps not found')

            try:
                time_step_value  =  np.round(timestamps[1] - timestamps[0], 2)  # time step value
            except IndexError:
                pass
            j  +=  1

        raw  =  AICSImage(fnames[0])                                                            # read and load file
        try:
            pix_size_xy  =  raw.physical_pixel_sizes.X
            pix_size_z   =  raw.physical_pixel_sizes.Z

        except:
            pix_size_xy, pix_size_z  =  ServiceWidgets.InputPixSize.get_vals()

        if analysis_folder is None:
            chs_red_green  =  ServiceWidgets.SetColorChannel.getChannels(raw.channel_names)
        else:
            chs_red_green  =  np.load(analysis_folder + '/chs_red_green.npy')
        red    =  raw.get_image_data("TZXY", C=chs_red_green[0])
        green  =  raw.get_image_data("TZXY", C=chs_red_green[1])

        tlen, zlen, xlen, ylen  =  raw.dims.T, raw.dims.Z, raw.dims.X, raw.dims.Y

        green_mip  =  np.zeros((tlen, xlen, ylen))
        red_mip    =  np.zeros((tlen, xlen, ylen))

        pbar  =  ServiceWidgets.ProgressBar(total1=tlen)
        pbar.show()

        for t in range(tlen):                                                # maximum intensity projection
            pbar.update_progressbar1(t)
            for x in range(xlen):
                red_mip[t, x, :]    =  red[t, :, x, :].max(0)
                green_mip[t, x, :]  =  green[t, :, x, :].max(0)

        pbar.close()

        if len(fnames) > 1:
            for fname in fnames[1:]:
                raw_bff    =  AICSImage(fname)  # read and load file
                red_bff    =  raw_bff.get_image_data("TZXY", C=chs_red_green[0])
                green_bff  =  raw_bff.get_image_data("TZXY", C=chs_red_green[1])

                tlen           =  raw_bff.dims.T
                red            =  np.concatenate((red, red_bff), axis=0)
                green          =  np.concatenate((green, green_bff), axis=0)
                green_bff_mip  =  np.zeros((tlen, xlen, ylen))
                red_bff_mip    =  np.zeros((tlen, xlen, ylen))

                pbar  =  ServiceWidgets.ProgressBar(total1=tlen)
                pbar.show()
                for t in range(tlen):                                                # maximum intensity projection
                    pbar.update_progressbar1(t)
                    for x in range(xlen):
                        red_bff_mip[t, x, :]    =  red_bff[t, :, x, :].max(0)
                        green_bff_mip[t, x, :]  =  green_bff[t, :, x, :].max(0)

                red_mip    =  np.concatenate((red_mip, red_bff_mip), axis=0)
                green_mip  =  np.concatenate((green_mip, green_bff_mip), axis=0)
                pbar.close()

        self.pix_size_xy      =  pix_size_xy
        self.pix_size_z       =  pix_size_z
        self.time_step_value  =  time_step_value
        self.red              =  red
        self.red_mip          =  red_mip
        self.green            =  green
        self.green_mip        =  green_mip
        self.chs_red_green    =  chs_red_green
