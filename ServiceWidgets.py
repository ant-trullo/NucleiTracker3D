"""Fun ction to use useful service widgets, like progressbars etc...

"""

import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt


class ProgressBar(QtWidgets.QWidget):
    """Simple progressbar widget."""
    def __init__(self, parent=None, total1=20):
        super().__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar1(self, val1):
        """Update progressbar."""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()


class ProgressBarDouble(QtWidgets.QWidget):
    """Double Progressbar widget."""
    def __init__(self, parent=None, total1=20, total2=20):
        super().__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        self.progressbar2  =  QtWidgets.QProgressBar()
        self.progressbar2.setMinimum(1)
        self.progressbar2.setMaximum(total2)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)
        main_layout.addWidget(self.progressbar2, 1, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar1(self, val1):
        """Update progressbar 1."""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()

    def update_progressbar2(self, val2):
        """Update progressbar 2."""
        self.progressbar2.setValue(val2)
        QtWidgets.qApp.processEvents()


class ChannelNumber(QtWidgets.QDialog):
    """Popup tool to input the channel number to work on"""
    def __init__(self, parent=None):
        super().__init__(parent)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        ch_numb_lbl  =  QtWidgets.QLabel("Spots channel", self)
        ch_numb_lbl.setFixedSize(int(ksf_h * 140), int(ksf_w * 25))

        choose_channel_combo  =  QtWidgets.QComboBox(self)
        choose_channel_combo.addItem("1")
        choose_channel_combo.addItem("2")
        choose_channel_combo.addItem("3")
        # choose_channel_combo.activated[str].connect(self.choose_channel_combo_var)
        choose_channel_combo.setCurrentIndex(0)
        choose_channel_combo.setFixedSize(int(ksf_h * 75), int(ksf_w * 25))

        input_close_btn  =  QtWidgets.QPushButton("Ok", self)
        input_close_btn.clicked.connect(self.input_close)
        input_close_btn.setToolTip('Input values')
        input_close_btn.setFixedSize(int(ksf_h * 50), int(ksf_w * 25))

        ch_numb_lbl_edit_box  =  QtWidgets.QHBoxLayout()
        ch_numb_lbl_edit_box.addWidget(ch_numb_lbl)
        ch_numb_lbl_edit_box.addWidget(choose_channel_combo)

        input_close_box  =  QtWidgets.QHBoxLayout()
        input_close_box.addStretch()
        input_close_box.addWidget(input_close_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(ch_numb_lbl_edit_box)
        layout.addLayout(input_close_box)

        self.choose_channel_combo  =  choose_channel_combo

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, 70, 50)
        self.setWindowTitle("Spot Nuc Max Distance")

    def input_close(self):
        """Close"""
        self.close()

    def ch_numb(self):
        """Return the channel number."""
        return int(self.choose_channel_combo.currentText()) - 1

    @staticmethod
    def getNumb(parent=None):
        """For signal sending."""
        dialog   =  ChannelNumber(parent)
        result   =  dialog.exec_()
        ch_numb  =  dialog.ch_numb()
        return ch_numb


class InputPixSize(QtWidgets.QDialog):
    """Popup tool to input pixel sizes in xy and z."""
    def __init__(self, parent=None):
        super().__init__(parent)

        ksf_h  =  np.load('keys_size_factor.npy')[0]
        ksf_w  =  np.load('keys_size_factor.npy')[1]

        pix_sizexy_lbl  =  QtWidgets.QLabel("Pix size x-y", self)
        pix_sizexy_lbl.setFixedSize(int(ksf_h * 140), int(ksf_w * 25))

        pix_sizez_lbl  =  QtWidgets.QLabel("Pix size z", self)
        pix_sizez_lbl.setFixedSize(int(ksf_h * 140), int(ksf_w * 25))

        pix_sizexy_edt  =  QtWidgets.QLineEdit(self)
        pix_sizexy_edt.setToolTip("Set pixels size in x-y")
        pix_sizexy_edt.setFixedSize(int(ksf_h * 30), int(ksf_w * 22))
        pix_sizexy_edt.textChanged[str].connect(self.pix_sizexy_var)

        pix_sizez_edt  =  QtWidgets.QLineEdit(self)
        pix_sizez_edt.setToolTip("Set pixels size in z")
        pix_sizez_edt.setFixedSize(int(ksf_h * 30), int(ksf_w * 22))
        pix_sizez_edt.textChanged[str].connect(self.pix_sizez_var)

        input_close_btn  =  QtWidgets.QPushButton("Ok", self)
        input_close_btn.clicked.connect(self.input_close)
        input_close_btn.setToolTip('Input values')
        input_close_btn.setFixedSize(int(ksf_h * 50), int(ksf_w * 25))

        pix_sizexy_lbl_edit_box  =  QtWidgets.QHBoxLayout()
        pix_sizexy_lbl_edit_box.addWidget(pix_sizexy_lbl)
        pix_sizexy_lbl_edit_box.addWidget(pix_sizexy_edt)

        pix_sizez_lbl_edit_box  =  QtWidgets.QHBoxLayout()
        pix_sizez_lbl_edit_box.addWidget(pix_sizez_lbl)
        pix_sizez_lbl_edit_box.addWidget(pix_sizez_edt)

        input_close_box  =  QtWidgets.QHBoxLayout()
        input_close_box.addStretch()
        input_close_box.addWidget(input_close_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(pix_sizexy_lbl_edit_box)
        layout.addLayout(pix_sizez_lbl_edit_box)
        layout.addLayout(input_close_box)

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, 70, 50)
        self.setWindowTitle("Spot Nuc Max Distance")

    def pix_sizexy_var(self, text):
        """Input the pixel size in xy."""
        self.pix_sizexy_value  =  float(text)

    def pix_sizez_var(self, text):
        """Input the pixel size in z."""
        self.pix_sizez_value  =  float(text)

    def input_close(self):
        """Input both pixel size values."""
        self.close()

    def pix_sizes(self):
        """Output the pixel size values."""
        return [self.pix_sizexy_value, self.pix_sizez_value]

    @staticmethod
    def get_vals(parent=None):
        """Send the output."""
        dialog     =  InputPixSize(parent)
        result     =  dialog.exec_()
        pix_sizes  =  dialog.pix_sizes()
        return pix_sizes


class SetColorChannel(QtWidgets.QDialog):
    """Set the color channels of the raw data to put in the gui."""
    def __init__(self, chnames_list, parent=None):
        super().__init__(parent)

        ksf_h         =  np.load('keys_size_factor.npy')[0]
        ksf_w         =  np.load('keys_size_factor.npy')[1]
        red_green_ch  =  np.load('red_green_ch.npy')

        red_channel_lbl  =  QtWidgets.QLabel("Red Channel", self)
        red_channel_lbl.setFixedSize(int(ksf_h * 120), int(ksf_w * 22))

        green_channel_lbl  =  QtWidgets.QLabel("Green Channel", self)
        green_channel_lbl.setFixedSize(int(ksf_h * 120), int(ksf_w * 22))

        red_channel_combo  =  QtWidgets.QComboBox(self)
        for chname in chnames_list:
            red_channel_combo.addItem(chname)
        red_channel_combo.activated.connect(self.red_channel_var)
        red_channel_combo.setCurrentIndex(red_green_ch[0])
        red_channel_combo.setFixedSize(int(ksf_h * 120), int(ksf_w * 25))

        green_channel_combo  =  QtWidgets.QComboBox(self)
        for chname2 in chnames_list:
            green_channel_combo.addItem(chname2)
        green_channel_combo.activated.connect(self.green_channel_var)
        green_channel_combo.setCurrentIndex(red_green_ch[1])
        green_channel_combo.setFixedSize(int(ksf_h * 120), int(ksf_w * 25))

        enter_values_btn  =  QtWidgets.QPushButton("OK", self)
        enter_values_btn.setToolTip("Set Channels Number")
        enter_values_btn.setFixedSize(int(ksf_h * 60), int(ksf_w * 25))
        enter_values_btn.clicked.connect(self.enter_values)

        red_box  =  QtWidgets.QHBoxLayout()
        red_box.addWidget(red_channel_lbl)
        red_box.addWidget(red_channel_combo)

        green_box  =  QtWidgets.QHBoxLayout()
        green_box.addWidget(green_channel_lbl)
        green_box.addWidget(green_channel_combo)

        enter_box  =  QtWidgets.QHBoxLayout()
        enter_box.addStretch()
        enter_box.addWidget(enter_values_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(red_box)
        layout.addLayout(green_box)
        layout.addLayout(enter_box)

        self.red_green_ch         =  np.copy(red_green_ch)
        self.red_channel_combo    =  red_channel_combo
        self.green_channel_combo  =  green_channel_combo

        self.red_channel_var()
        self.green_channel_var()

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, 250, 150)
        self.setWindowTitle("Set Channels")

    def red_channel_var(self):
        """Set red channel."""
        self.red_channel_value  =  self.red_channel_combo.currentIndex()

    def green_channel_var(self):
        """Set green channel."""
        self.green_channel_value  =  self.green_channel_combo.currentIndex()

    def enter_values(self):
        """Organizing channels info for the output."""
        self.close()

    def params(self):
        """Function to send results."""
        return [self.red_channel_value, self.green_channel_value]

    @staticmethod
    def getChannels(parent=None):
        """Send results."""
        dialog  =  SetColorChannel(parent)
        result  =  dialog.exec_()
        flag    =  dialog.params()
        return flag
