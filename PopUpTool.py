"""This function provides only a pg.image frame with frame number.


"""


import pyqtgraph as pg


class PopUpTool:
    """Update frame number text on a image without colormap"""
    def __init__(self, img, img_title):

        w  =  pg.image(img, title=img_title)

        lx  =  []
        for t in range(img.shape[0]):
            lx.append(pg.TextItem(str(t), color='r', anchor=(0, 1)))
        w.addItem(lx[0])
        self.current  =  0

        def upadateframe():
            w.removeItem(lx[self.current])
            w.addItem(lx[w.currentIndex])
            self.current  =  w.currentIndex

        w.timeLine.sigPositionChanged.connect(upadateframe)


class PopUpToolWithMap:
    """Update frame number text on a image with colormap"""
    def __init__(self, img, img_title, cmap):

        w  =  pg.image(img, title=img_title)
        w.setColorMap(cmap)

        lx = []
        for t in range(img.shape[0]):
            lx.append(pg.TextItem(str(t), color='r', anchor=(0, 1)))
        w.addItem(lx[0])
        self.current  =  0

        def upadateframe():
            w.removeItem(lx[self.current])
            w.addItem(lx[w.currentIndex])
            self.current  =  w.currentIndex

        w.timeLine.sigPositionChanged.connect(upadateframe)
