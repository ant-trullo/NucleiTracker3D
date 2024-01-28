"""This is the main window of the software that detects and follows 3D nuclei in drosophila embryo.
This is version 1.0, since February 2023.

developer: antonio.trullo@igmm.cnrs.fr

"""

import sys
import os.path
import time
import traceback
from importlib import reload
from natsort import natsorted
import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore

import LoadRawDataCzi
import NucleiSegmenter3D
import Nuclei3dTracker
import AnalysisSaver
import AnalysisLoader
import PopUpTool
import ServiceWidgets


class MainWindow(QtWidgets.QMainWindow):
    """Main windows: coordinates all the actions, algorithms, visualization tools and analysis tools."""
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        mycmap = np.fromfile("mycmap.bin", "uint16").reshape((10000, 3))  # / 255.0
        colors4map = []
        for k in range(mycmap.shape[0]):
            colors4map.append(mycmap[k, :])
        colors4map[0] = np.array([0, 0, 0])

        widget = QtWidgets.QWidget(self)
        self.setCentralWidget(widget)

        load_data_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/load-hi.png'), "&Load Data", self)
        load_data_action.setShortcut("Ctrl+L")
        load_data_action.setStatusTip("Load raw data files with a single channel")
        load_data_action.triggered.connect(self.load_raw_data)

        settings_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/settings.png'), "&Settings", self)
        settings_action.setShortcut("Ctrl+T")
        settings_action.setStatusTip("Changes default settings values")
        settings_action.triggered.connect(self.settings_changes)

        save_analysis_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Save Analysis", self)
        save_analysis_action.setShortcut("Ctrl+S")
        save_analysis_action.setStatusTip("Save analysis")
        save_analysis_action.triggered.connect(self.save_analysis)

        load_analysis_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Load Analysis", self)
        load_analysis_action.setShortcut("Ctrl+A")
        load_analysis_action.setStatusTip("Load analysis")
        load_analysis_action.triggered.connect(self.load_analysis)

        spatial_analysis_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/geo_pin.jpg'), "&Spatial Analysis", self)
        spatial_analysis_action.setShortcut("Ctrl+W")
        spatial_analysis_action.setStatusTip("Spatial Analysis")
        spatial_analysis_action.triggered.connect(self.spatial_analysis)

        exit_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/exit.png'), "&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.close)

        popup_raw_mip_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/popup.jpg'), "&Pop-up Raw Mip Nuclei", self)
        popup_raw_mip_action.setStatusTip("Launch popup tool to visualize mip raw data")
        popup_raw_mip_action.triggered.connect(self.popup_raw_mip)

        menubar   =  self.menuBar()

        file_menu  =  menubar.addMenu("&File")
        file_menu.addAction(load_data_action)
        file_menu.addAction(save_analysis_action)
        file_menu.addAction(load_analysis_action)
        file_menu.addAction(settings_action)
        file_menu.addAction(exit_action)

        file_menu  =  menubar.addMenu("&View")
        file_menu.addAction(popup_raw_mip_action)

        file_menu  =  menubar.addMenu("&Post-Processing")
        file_menu.addAction(spatial_analysis_action)

        fname_raw_edt  =  QtWidgets.QLineEdit(" ", self)
        fname_raw_edt.setToolTip("Names of the files you are working on")

        frame_raw_red  =  pg.ImageView(self, name="RawRed")
        frame_raw_red.ui.roiBtn.hide()
        frame_raw_red.ui.menuBtn.hide()
        frame_raw_red.view.setXLink("RawGreen")
        frame_raw_red.view.setYLink("RawGreen")

        frame_raw_green  =  pg.ImageView(self, name="RawGreen")
        frame_raw_green.ui.roiBtn.hide()
        frame_raw_green.ui.menuBtn.hide()
        frame_raw_green.view.setXLink("Segm")
        frame_raw_green.view.setYLink("Segm")

        frame_sgm  =  pg.ImageView(self, name="Segm")
        frame_sgm.ui.roiBtn.hide()
        frame_sgm.ui.menuBtn.hide()
        frame_sgm.view.setXLink("RawRed")
        frame_sgm.view.setYLink("RawRed")
        frame_sgm.getImageItem().mouseClickEvent  =  self.click

        framepp1  =  pg.PlotWidget(self)

        frame_sgm_box  =  QtWidgets.QHBoxLayout()
        frame_sgm_box.addWidget(frame_sgm)

        tab_sgm   =  QtWidgets.QTabWidget()
        tab1_sgm  =  QtWidgets.QWidget()
        tab_sgm.addTab(tab1_sgm, "Segmented")

        tab1_sgm.setLayout(frame_sgm_box)

        tabs_raw  =  QtWidgets.QTabWidget()
        tab1_raw  =  QtWidgets.QWidget()
        tab2_raw  =  QtWidgets.QWidget()

        frame_raw_red_box  =  QtWidgets.QHBoxLayout()
        frame_raw_red_box.addWidget(frame_raw_red)

        frame_raw_green_box  =  QtWidgets.QHBoxLayout()
        frame_raw_green_box.addWidget(frame_raw_green)

        tab1_raw.setLayout(frame_raw_red_box)
        tab2_raw.setLayout(frame_raw_green_box)

        tabs_raw.addTab(tab1_raw, "Raw Red")
        tabs_raw.addTab(tab2_raw, "Raw Green")

        frames_box  =  QtWidgets.QHBoxLayout()
        frames_box.addWidget(tab_sgm)
        frames_box.addWidget(tabs_raw)

        sld_time  =  QtWidgets.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_time.valueChanged.connect(self.sld_time_change)

        sld_zed  =  QtWidgets.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_zed.valueChanged.connect(self.sld_zed_change)

        busy_lbl  =  QtWidgets.QLabel("Ready")
        busy_lbl.setStyleSheet("color: green")

        pixsize_x_lbl  =  QtWidgets.QLabel("pix size XY =;")
        pixsize_z_lbl  =  QtWidgets.QLabel("Z step =;")
        time_step_lbl  =  QtWidgets.QLabel("Time step =")

        bottom_labels_box  =  QtWidgets.QHBoxLayout()
        bottom_labels_box.addWidget(busy_lbl)
        bottom_labels_box.addStretch()
        bottom_labels_box.addWidget(pixsize_x_lbl)
        bottom_labels_box.addWidget(pixsize_z_lbl)
        bottom_labels_box.addWidget(time_step_lbl)

        first_frame_btn  =  QtWidgets.QPushButton("Start", self)
        first_frame_btn.setToolTip("First frame to analyse")
        first_frame_btn.setFixedSize(int(ksf_w * 75), int(ksf_h * 25))
        first_frame_btn.clicked.connect(self.first_frame)

        last_frame_btn  =  QtWidgets.QPushButton("End", self)
        last_frame_btn.setToolTip("Last frame to analyse")
        last_frame_btn.setFixedSize(int(ksf_w * 75), int(ksf_h * 25))
        last_frame_btn.clicked.connect(self.last_frame)

        first_last_box  =  QtWidgets.QHBoxLayout()
        first_last_box.addWidget(first_frame_btn)
        first_last_box.addWidget(last_frame_btn)

        segment_btn  =  QtWidgets.QPushButton("Segment", self)
        segment_btn.setToolTip("Segment nuclei in 3D")
        segment_btn.setFixedSize(int(ksf_w * 150), int(ksf_h * 25))
        segment_btn.clicked.connect(self.segment)

        track_btn  =  QtWidgets.QPushButton("Track", self)
        track_btn.setToolTip("Track the segmented nuclei")
        track_btn.setFixedSize(int(ksf_w * 150), int(ksf_h * 25))
        track_btn.clicked.connect(self.track)

        time_lbl  =  QtWidgets.QLabel("time     " + '0', self)
        time_lbl.setFixedSize(int(ksf_w * 150), int(ksf_h * 25))

        zed_lbl  =  QtWidgets.QLabel("Z           " + '0', self)
        zed_lbl.setFixedSize(int(ksf_w * 150), int(ksf_h * 25))

        frame_numb_lbl = QtWidgets.QLabel("frame  " + '0', self)
        frame_numb_lbl.setFixedSize(int(ksf_w * 150), int(ksf_h * 25))

        hor_line  =  QtWidgets.QFrame()
        hor_line.setFrameStyle(QtWidgets.QFrame.HLine)
        hor_line.setFixedWidth(int(ksf_w * 150))

        selectpop_flag_checkbox  =  QtWidgets.QCheckBox(self)
        selectpop_flag_checkbox.setFixedSize(int(ksf_h * 150), int(ksf_h * 25))
        selectpop_flag_checkbox.setText("Traces")
        selectpop_flag_checkbox.setToolTip("Activate the option to plot nucleus trace in both channels")

        commands  =  QtWidgets.QVBoxLayout()
        commands.addLayout(first_last_box)
        commands.addWidget(segment_btn)
        commands.addWidget(track_btn)
        commands.addStretch()
        commands.addWidget(selectpop_flag_checkbox)
        commands.addStretch()
        commands.addWidget(zed_lbl)
        commands.addWidget(time_lbl)
        commands.addWidget(frame_numb_lbl)
        commands.addWidget(hor_line)

        vert_box  =  QtWidgets.QVBoxLayout()
        vert_box.addWidget(fname_raw_edt)
        vert_box.addLayout(frames_box)
        vert_box.addWidget(sld_time)
        vert_box.addWidget(sld_zed)
        vert_box.addWidget(framepp1)

        hor_box  =  QtWidgets.QHBoxLayout()
        hor_box.addLayout(vert_box)
        hor_box.addLayout(commands)

        layout  =  QtWidgets.QVBoxLayout(widget)
        layout.addLayout(hor_box)
        layout.addLayout(bottom_labels_box)

        self.soft_version             =  "NucleiTracker3D_v1.0"
        self.busy_lbl                 =  busy_lbl
        self.sld_time                 =  sld_time
        self.sld_zed                  =  sld_zed
        self.fname_raw_edt            =  fname_raw_edt
        self.pixsize_x_lbl            =  pixsize_x_lbl
        self.pixsize_z_lbl            =  pixsize_z_lbl
        self.time_step_lbl            =  time_step_lbl
        self.time_lbl                 =  time_lbl
        self.frame_numb_lbl           =  frame_numb_lbl
        self.framepp1                 =  framepp1
        self.zed_lbl                  =  zed_lbl
        self.frame_raw_green          =  frame_raw_green
        self.frame_raw_red            =  frame_raw_red
        self.frame_sgm                =  frame_sgm
        self.colors4map               =  colors4map
        self.flag_trck                =  False
        self.selectpop_flag_checkbox  =  selectpop_flag_checkbox

        self.setGeometry(800, 100, 700, 500)
        self.setWindowTitle(self.soft_version)
        self.setWindowIcon(QtGui.QIcon('Icons/DrosophilaIcon.png'))
        self.show()

    def closeEvent(self, event):
        """Close the GUI, asking confirmation."""
        quit_msg  =  "Are you sure you want to exit the program?"
        reply     =  QtWidgets.QMessageBox.question(self, 'Message', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def busy_indicator(self):
        """Write a red text (BUSY) as a label on the GUI (bottom left)."""
        self.busy_lbl.setText("Busy")
        self.busy_lbl.setStyleSheet('color: red')

    def ready_indicator(self):
        """Write a green text (READY) as a label on the GUI (bottom left)."""
        self.busy_lbl.setText("Ready")
        self.busy_lbl.setStyleSheet('color: green')

    def settings_changes(self):
        """Change settings."""
        self.mpp2 = SettingsChanges()
        self.mpp2.show()
        self.mpp2.procStart.connect(self.settings_update)

    def settings_update(self):
        """Restart the GUI to make changes in button size effective."""
        self.mpp2.close()
        os.execl(sys.executable, sys.executable, *sys.argv)

    def sld_time_change(self):
        """Update frames when time cursor changes position."""
        self.frame_raw_red.setImage(self.raw_data.red[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        self.frame_raw_green.setImage(self.raw_data.green[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        self.frame_numb_lbl.setText("Frame " + str(self.sld_time.value()))
        self.time_lbl.setText("time  " + time.strftime("%M:%S", time.gmtime(self.sld_time.value() * self.raw_data.time_step_value)))
        try:
            if not self.flag_trck:
                self.frame_sgm.setImage(self.nucs_segm[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
            elif self.flag_trck:
                self.frame_sgm.setImage(self.nucs_trck[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        except AttributeError:
            pass

    def sld_zed_change(self):
        """Update frames when zed cursor changes position."""
        self.frame_raw_red.setImage(self.raw_data.red[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        self.frame_raw_green.setImage(self.raw_data.green[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        self.zed_lbl.setText("Z         " + str(self.sld_zed.value()))
        try:
            if not self.flag_trck:
                self.frame_sgm.setImage(self.nucs_segm[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
            elif self.flag_trck:
                self.frame_sgm.setImage(self.nucs_trck[self.sld_time.value(), self.sld_zed.value()], autoRange=False, autoLevels=False)
        except AttributeError:
            pass

    def first_frame(self):
        """Define the first frame to analyze."""
        self.raw_data.green      =  self.raw_data.green[self.sld_time.value():]
        self.raw_data.red        =  self.raw_data.red[self.sld_time.value():]
        self.raw_data.green_mip  =  self.raw_data.green_mip[self.sld_time.value():]
        self.raw_data.red_mip    =  self.raw_data.red_mip[self.sld_time.value():]
        self.sld_time.setMaximum(self.raw_data.red.shape[0] - 1)
        self.sld_time.setValue(0)

    def last_frame(self):
        """Define the last frame to analyze."""
        self.raw_data.green      =  self.raw_data.green[:self.sld_time.value()]
        self.raw_data.red        =  self.raw_data.red[:self.sld_time.value()]
        self.raw_data.green_mip  =  self.raw_data.green_mip[:self.sld_time.value()]
        self.raw_data.red_mip    =  self.raw_data.red_mip[:self.sld_time.value()]
        self.sld_time.setMaximum(self.raw_data.red.shape[0] - 1)
        self.sld_time.setValue(self.raw_data.red.shape[0] - 1)

    def popup_raw_mip(self):
        """Popup mip raw data."""
        PopUpTool.PopUpTool(self.raw_data.red_mip, "Red Raw Mip")
        PopUpTool.PopUpTool(self.raw_data.green_mip, "Green Raw Mip")

    def load_raw_data(self):
        """Load raw data."""
        reload(LoadRawDataCzi)
        self.fnames         =  natsorted(QtWidgets.QFileDialog.getOpenFileNames(None, "Select czi (or lsm) data files to concatenate...", filter="*.lsm *.czi *.tif *.lif")[0])
        # self.fnames         =  ['/home/atrullo/Dropbox/Virginia_Anto/snail-llama/09012023_E1/09012023_snaLlama_HisRFP_E1f.czi']
        self.fname4journal  =  []
        joined_fnames       =  ' '
        for fname in self.fnames:
            joined_fnames  +=  str(fname[fname.rfind('/') + 1:]) +  ' ~~~ '
            self.fname4journal.append(str(fname[fname.rfind('/') + 1:]))
        self.fname_raw_edt.setText(joined_fnames)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:
            self.frame_sgm.clear()
            self.frame_raw_red.clear()
            self.frame_raw_green.clear()
            self.raw_data       =  LoadRawDataCzi.LoadRawDataCzi(self.fnames)
            self.frame_raw_red.setImage(self.raw_data.red[0, 0], autoRange=False, autoLevels=True)
            self.frame_raw_green.setImage(self.raw_data.green[0, 0], autoRange=False, autoLevels=True)
            self.pixsize_x_lbl.setText("pix size XY = " + str(np.round(self.raw_data.pix_size_xy, decimals=4)) + "µm;")
            self.pixsize_z_lbl.setText("Z step = " + str(np.round(self.raw_data.pix_size_z, decimals=4)) + "µm;")
            self.time_step_lbl.setText("Time Step = " + str(np.round(self.raw_data.time_step_value, decimals=4)) + "s")
            self.sld_time.setMaximum(self.raw_data.red.shape[0] - 1)
            self.sld_time.setValue(0)
            self.sld_zed.setMaximum(self.raw_data.red.shape[1] - 1)
            self.sld_zed.setValue(0)

        except Exception:
            traceback.print_exc()
        self.ready_indicator()

    def segment(self):
        """Launch 3D nuclei segmentation."""
        reload(NucleiSegmenter3D)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:
            self.frame_sgm.clear()
            self.nucs_segm  =  NucleiSegmenter3D.NucleiSegmenter3D(self.raw_data.red).nucs_sgm
            self.frame_sgm.setImage(self.nucs_segm[self.sld_time.value(), self.sld_zed.value()], autoRange=False)
            nucs_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_segm.max()), color=self.colors4map)
            self.frame_sgm.setColorMap(nucs_cmap)
            self.flag_trck  =  False
        except Exception:
            traceback.print_exc()
        self.ready_indicator()

    def track(self):
        """Launch nuclei tracking."""
        reload(Nuclei3dTracker)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.frame_sgm.clear()
            # self.nucs_trck  =  Nuclei3dTracker.Nuclei3dTracker(self.nucs_segm, self.dist_thr_value, self.raw_data.pix_size_xy, self.raw_data.pix_size_z).nucs_trck
            self.nucs_trck  =  Nuclei3dTracker.NucleiOverlTracker(self.nucs_segm).nucs_trck
            self.frame_sgm.setImage(self.nucs_trck[self.sld_time.value(), self.sld_zed.value()], autoRange=False)
            nucs_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_trck.max()), color=self.colors4map)
            self.frame_sgm.setColorMap(nucs_cmap)
            self.flag_trck  =  True
        except Exception:
            traceback.print_exc()
        self.ready_indicator()

    def save_analysis(self):
        """Save analysis."""
        reload(AnalysisSaver)
        folder2write  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select or Define a Folder to Store the analysis results"))
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            # AnalysisSaver.AnalysisSaver(self.nucs_segm, self.nucs_trck, self.raw_data, folder2write, self.dist_thr_value, self.chs_red_green, self.soft_version, self.fname4journal)
            AnalysisSaver.AnalysisSaver(self.nucs_segm, self.nucs_trck, self.raw_data, folder2write, self.raw_data.chs_red_green, self.soft_version, self.fname4journal)
        except Exception:
            traceback.print_exc()
        self.ready_indicator()

    def load_analysis(self):
        """Load analysis."""
        reload(AnalysisLoader)
        analysis_folder     =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the analysis Folder"))
        self.fnames         =  natsorted(QtWidgets.QFileDialog.getOpenFileNames(None, "Select czi (or lsm) data files to concatenate...", filter="*.lsm *.czi *.tif *.lif")[0])
        self.fname4journal  =  []
        joined_fnames       =  ' '
        for fname in self.fnames:
            joined_fnames  +=  str(fname[fname.rfind('/') + 1:]) +  ' ~~~ '
            self.fname4journal.append(str(fname[fname.rfind('/') + 1:]))

        self.fname_raw_edt.setText(joined_fnames)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.nucs_trck      =  np.load(analysis_folder + '/nucs_trck.npy')
            self.nucs_segm      =  np.load(analysis_folder + '/nucs_segm.npy')
            # self.chs_red_green  =  np.load(analysis_folder + '/chs_red_green.npy')
            self.frame_sgm.clear()
            self.frame_sgm.setImage(self.nucs_trck[self.sld_time.value(), self.sld_zed.value()], autoRange=True, autoLevels=True)
            nucs_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_trck.max()), color=self.colors4map)
            self.frame_sgm.setColorMap(nucs_cmap)
            self.frame_sgm.setLevels(min=0, max=self.nucs_trck.max())
            self.flag_trck  =  True
            self.raw_data   =  AnalysisLoader.RawDataLoader(analysis_folder, self.fnames)
            self.frame_raw_red.setImage(self.raw_data.red[0, 0], autoRange=False, autoLevels=True)
            self.frame_raw_green.setImage(self.raw_data.green[0, 0], autoRange=False, autoLevels=True)
            self.pixsize_x_lbl.setText("pix size XY = " + str(np.round(self.raw_data.pix_size_xy, decimals=4)) + "µm;")
            self.pixsize_z_lbl.setText("Z step = " + str(np.round(self.raw_data.pix_size_z, decimals=4)) + "µm;")
            self.time_step_lbl.setText("Time Step = " + str(np.round(self.raw_data.time_step_value, decimals=4)) + "s")
            self.crop_roi  =  [0, 0, self.raw_data.red_mip.shape[1], self.raw_data.red_mip.shape[2]]
            self.sld_time.setMaximum(self.raw_data.red.shape[0] - 1)
            self.sld_time.setValue(0)
            self.sld_zed.setMaximum(self.raw_data.red.shape[1] - 1)
            self.sld_zed.setValue(0)
        except Exception:
            traceback.print_exc()
        self.ready_indicator()

    def click(self, event):
        """Plot time traces of a selected nucleus in both channels."""
        reload(AnalysisSaver)
        event.accept()
        pos        =  event.pos()
        modifiers  =  QtWidgets.QApplication.keyboardModifiers()

        if self.selectpop_flag_checkbox.isChecked():
            if modifiers  ==  QtCore.Qt.ShiftModifier:
                nucs_tag        =  self.nucs_trck[self.sld_time.value(), self.sld_zed.value(), int(pos.x()), int(pos.y())]
                red_green2plot  =  AnalysisSaver.ShowTraces(self.nucs_trck, nucs_tag, self.raw_data)
                self.framepp1.clear()
                self.framepp1.plot(red_green2plot.red_avints, pen='r', symbol='o')
                self.framepp1.plot(red_green2plot.green_avints, pen='g', symbol='x')

    def spatial_analysis(self):
        """Run the spatial analysis tool."""
        analysis_folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the analysis Folder"))
        fnames           =  natsorted(QtWidgets.QFileDialog.getOpenFileNames(None, "Select czi (or lsm) data files to concatenate...", filter="*.lsm *.czi *.tif *.lif")[0])
        fname4journal  =  []
        joined_fnames       =  ' '
        for fname in fnames:
            joined_fnames  +=  str(fname[fname.rfind('/') + 1:]) +  ' ~~~ '
            fname4journal.append(str(fname[fname.rfind('/') + 1:]))
        raw_data         =  AnalysisLoader.RawDataLoader(analysis_folder, fnames)
        self.mpp3        =  SpatialAnalysisTool(raw_data, analysis_folder, self.soft_version, fname4journal)
        self.mpp3.show()


class SettingsChanges(QtWidgets.QWidget):
    """Tool to change visualization and analysis parameters."""
    procStart  =  QtCore.pyqtSignal()

    def __init__(self):
        QtWidgets.QWidget.__init__(self)

        self.ksf_h         =  np.load('keys_size_factor.npy')[0]
        self.ksf_w         =  np.load('keys_size_factor.npy')[1]
        self.red_green_ch  =  np.load('red_green_ch.npy')

        ksf_h_lbl  =  QtWidgets.QLabel("Keys Scale Factor W")

        ksf_h_edt  =  QtWidgets.QLineEdit(self)
        ksf_h_edt.textChanged[str].connect(self.ksf_h_var)
        ksf_h_edt.setToolTip("Sets keys scale size (width)")
        ksf_h_edt.setFixedSize(int(self.ksf_h * 50), int(self.ksf_w * 25))
        ksf_h_edt.setText(str(self.ksf_h))

        ksf_w_lbl  =  QtWidgets.QLabel("Keys Scale Factor H")

        ksf_w_edt  =  QtWidgets.QLineEdit(self)
        ksf_w_edt.textChanged[str].connect(self.ksf_w_var)
        ksf_w_edt.setToolTip("Sets keys scale size (heigth)")
        ksf_w_edt.setFixedSize(int(self.ksf_h * 50), int(self.ksf_w * 25))
        ksf_w_edt.setText(str(self.ksf_w))

        red_ch_numb_lbl  =  QtWidgets.QLabel("Red channel")

        red_ch_numb_combo  =  QtWidgets.QComboBox()
        red_ch_numb_combo.addItem("1")
        red_ch_numb_combo.addItem("2")
        red_ch_numb_combo.setToolTip('Set the default channel number of the red')
        red_ch_numb_combo.setFixedSize(int(self.ksf_h * 65), int(self.ksf_w * 25))
        red_ch_numb_combo.activated[str].connect(self.red_ch_numb_var)
        red_ch_numb_combo.setCurrentIndex(self.red_green_ch[0])
        self.red_ch_numb_var(red_ch_numb_combo.currentText())

        green_ch_numb_lbl  =  QtWidgets.QLabel("Green channel")

        green_ch_numb_combo  =  QtWidgets.QComboBox()
        green_ch_numb_combo.addItem("1")
        green_ch_numb_combo.addItem("2")
        green_ch_numb_combo.setToolTip('Set the default channel number of the green')
        green_ch_numb_combo.setFixedSize(int(self.ksf_h * 65), int(self.ksf_w * 25))
        green_ch_numb_combo.activated[str].connect(self.green_ch_numb_var)
        green_ch_numb_combo.setCurrentIndex(self.red_green_ch[1])
        self.green_ch_numb_var(green_ch_numb_combo.currentText())

        save_btn  =  QtWidgets.QPushButton("Save", self)
        save_btn.clicked.connect(self.save_vars)
        save_btn.setToolTip('Make default the choseen parameters')
        save_btn.setFixedSize(int(self.ksf_h * 55), int(self.ksf_w * 25))

        close_btn  =  QtWidgets.QPushButton("Close", self)
        close_btn.clicked.connect(self.close_)
        close_btn.setToolTip('Close Widget')
        close_btn.setFixedSize(int(self.ksf_h * 55), int(self.ksf_w * 25))

        restart_btn  =  QtWidgets.QPushButton("Refresh", self)
        restart_btn.clicked.connect(self.restart)
        restart_btn.setToolTip('Refresh GUI')
        restart_btn.setFixedSize(int(self.ksf_h * 75), int(self.ksf_w * 25))

        btns_box  =  QtWidgets.QHBoxLayout()
        btns_box.addWidget(save_btn)
        btns_box.addWidget(restart_btn)
        btns_box.addWidget(close_btn)

        layout_grid  =  QtWidgets.QGridLayout()
        layout_grid.addWidget(ksf_h_lbl, 0, 0)
        layout_grid.addWidget(ksf_h_edt, 0, 1)
        layout_grid.addWidget(ksf_w_lbl, 1, 0)
        layout_grid.addWidget(ksf_w_edt, 1, 1)
        layout_grid.addWidget(red_ch_numb_lbl, 2, 0)
        layout_grid.addWidget(red_ch_numb_combo, 2, 1)
        layout_grid.addWidget(green_ch_numb_lbl, 3, 0)
        layout_grid.addWidget(green_ch_numb_combo, 3, 1)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(layout_grid)
        layout.addStretch()
        layout.addLayout(btns_box)

        self.setLayout(layout)
        self.setGeometry(300, 300, 60, 60)
        self.setWindowTitle("Settings Tool")

    def ksf_h_var(self, text):
        """Set keys size factor value (hight)."""
        self.ksf_h  =  np.float64(text)

    def ksf_w_var(self, text):
        """Set keys size factor value (width)."""
        self.ksf_w  =  np.float64(text)

    def red_ch_numb_var(self, text):
        """Set the red channel number."""
        self.red_ch_numb_value  =  int(text) - 1

    def green_ch_numb_var(self, text):
        """Set the green channel number."""
        self.green_ch_numb_value  =  int(text) - 1

    def save_vars(self):
        """Save new settings."""
        np.save('keys_size_factor.npy', [self.ksf_h, self.ksf_w])
        np.save('red_green_ch.npy', [self.red_ch_numb_value, self.green_ch_numb_value])

    def close_(self):
        """Close the widget."""
        self.close()

    @QtCore.pyqtSlot()
    def restart(self):
        """Send message to main GUI."""
        self.procStart.emit()


class SpatialAnalysisTool(QtWidgets.QWidget):
    """Tool to change visualization and analysis parameters."""
    def __init__(self, raw_data, analysis_folder, soft_version, fname4journal):
        QtWidgets.QWidget.__init__(self)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        raw_data_3c          =  np.zeros(raw_data.red_mip.shape + (3,), dtype=raw_data.red_mip.dtype)
        raw_data_3c[..., 0]  =  raw_data.red_mip
        raw_data_3c[..., 1]  =  raw_data.green_mip

        analysis_folder_lbl  =  QtWidgets.QLabel("Folder: " +  analysis_folder[analysis_folder.rfind("/") + 1:])
        analysis_folder_lbl.setToolTip("Name of the analysis folder you are working on")

        frame_raw_3ch  =  pg.ImageView(self)
        frame_raw_3ch.ui.roiBtn.hide()
        frame_raw_3ch.ui.menuBtn.hide()
        frame_raw_3ch.timeLine.sigPositionChanged.connect(self.update_frame_counter)
        frame_raw_3ch.setImage(raw_data_3c)

        roi_line  =  pg.LinearRegionItem(orientation=True, brush=[0, 0, 0, 0])
        frame_raw_3ch.addItem(roi_line)

        frame_lbl_box  =  QtWidgets.QVBoxLayout()
        frame_lbl_box.addWidget(analysis_folder_lbl)
        frame_lbl_box.addWidget(frame_raw_3ch)

        save_btn  =  QtWidgets.QPushButton("Save", self)
        save_btn.clicked.connect(self.save_spatial)
        save_btn.setToolTip("Save data spatially organized")
        save_btn.setFixedSize(int(ksf_h * 140), int(ksf_w * 25))

        frame_numb_nucs_lbl  =  QtWidgets.QLabel("Frame ", self)
        frame_numb_nucs_lbl.setFixedSize(int(ksf_h * 140), int(ksf_h * 25))

        commands  =  QtWidgets.QVBoxLayout()
        commands.addWidget(save_btn)
        commands.addStretch()
        commands.addWidget(frame_numb_nucs_lbl)

        layout  =  QtWidgets.QHBoxLayout()
        layout.addLayout(frame_lbl_box)
        layout.addLayout(commands)

        self.frame_raw_3ch        =  frame_raw_3ch
        self.frame_numb_nucs_lbl  =  frame_numb_nucs_lbl
        self.roi_line             =  roi_line
        self.raw_data             =  raw_data
        self.analysis_folder      =  analysis_folder
        self.soft_version         =  soft_version
        self.fname4journal        =  fname4journal

        self.setLayout(layout)
        self.setGeometry(300, 300, 500, 500)
        self.setWindowTitle("Spatial Analysis Tool")

    def update_frame_counter(self):
        """Update frame number counter value."""
        self.frame_numb_nucs_lbl.setText("Frame " + str(self.frame_raw_3ch.currentIndex))

    def save_spatial(self):
        """Save results spatially organized."""
        reload(AnalysisSaver)
        AnalysisSaver.SaveTracesSpatial(self.analysis_folder, self.raw_data, np.load(self.analysis_folder + '/nucs_segm.npy'), self.roi_line, self.soft_version, self.fname4journal)

def main():
    app         =  QtWidgets.QApplication(sys.argv)
    splash_pix  =  QtGui.QPixmap('Icons/DrosophilaIcon.png')
    splash      =  QtWidgets.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()
    ex  =  MainWindow()
    splash.finish(ex)
    sys.exit(app.exec_())


def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)


if __name__ == '__main__':

    main()
    sys.excepthook = except_hook
