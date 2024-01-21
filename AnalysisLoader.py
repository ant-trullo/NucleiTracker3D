"""Ths function loads a previously done analysis.

Input are the analysis folder path and raw data path(s).
"""


import numpy as np
from openpyxl import load_workbook

import LoadRawDataCzi


class RawDataLoader:
    """This class loads raw data."""
    def __init__(self, analysis_folder, fnames):

        raw_data           =  LoadRawDataCzi.LoadRawDataCzi(fnames, analysis_folder)
        first_slctd_frame  =  np.load(analysis_folder + '/first_slctd_frame.npy')
        last_slctd_frame   =  np.load(analysis_folder + '/last_slctd_frame.npy')
        jj_start           =  np.where(np.sum(raw_data.red_mip - first_slctd_frame, axis=(1, 2)) == 0)[0][0]
        jj_end             =  np.where(np.sum(raw_data.red_mip - last_slctd_frame, axis=(1, 2)) == 0)[0][0]

        raw_data.red_mip    =  raw_data.red_mip[jj_start:jj_end + 1]
        raw_data.red        =  raw_data.red[jj_start:jj_end + 1]
        raw_data.green_mip  =  raw_data.green_mip[jj_start:jj_end + 1]
        raw_data.green      =  raw_data.green[jj_start:jj_end + 1]

        self.green            =  raw_data.green
        self.red              =  raw_data.red
        self.green_mip        =  raw_data.green_mip
        self.red_mip          =  raw_data.red_mip
        self.pix_size_xy      =  raw_data.pix_size_xy
        self.pix_size_z       =  raw_data.pix_size_z
        self.time_step_value  =  raw_data.time_step_value
        self.ch_red_green     =  np.load(analysis_folder + '/ch_red_nucs.npy')