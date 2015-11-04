import re
from PyQt4.QtGui import (QGraphicsRectItem, QGraphicsLineItem,
                         QGraphicsEllipseItem, QPen, QColor, QBrush,
                         QFont, QPixmap, QFontMetrics, QPainter,
                         QGraphicsSimpleTextItem, QGraphicsItem)

from PyQt4.QtCore import Qt,  QPointF, QRect, QRectF
from ete3 import faces
from ete3.treeview.main import *
import math


class PPieChartFace(faces.StaticItemFace):
    """
    .. versionadded:: 2.2
    :param percents: a list of values summing up 100.
    :param width: width of the piechart
    :param height: height of the piechart
    :param colors: a list of colors (same length as percents)
    :param line_color: color used to render the border of the piechart (None=transparent)
    """
    def __init__(self, percents, width, height, colors=None, line_color=None, label_size=2, labels=None, is_percent=False):
        faces.Face.__init__(self)


        if(labels==[]):
            labels = [str(x) for x in percents]

        if not is_percent:
            s = sum(percents)
            percents = map(lambda x: x*100./s, percents)
       
        if round(sum(percents)) > 100:
            raise ValueError("PPieChartItem: percentage values > 100")

        if(len(percents)!= len(labels)):
            raise ValueError("PPieChartItem : label and data have differents size")

        self.type = "item"
        self.item = None
        self.percents = percents
        self.labels = labels
        if not colors:
            colors = COLOR_SCHEMES["paired"]
        self.colors =  colors
        self.width = width
        self.height = height
        self.label_size = label_size
        self.line_color = line_color

    def update_items(self):
        self.item = _PieChartItem(self.percents, self.width,
                    self.height, self.colors, self.line_color, self.labels, self.label_size)

    def _width(self):
        return self.item.rect().width()

    def _height(self):
        return self.item.rect().height()


class _PieChartItem(QGraphicsRectItem):
    def __init__(self, percents, width, height, colors, line_color=None, labels=None, label_size=2):
        QGraphicsRectItem.__init__(self, 0, 0, width, height)
        self.percents = percents
        self.labels = labels
        self.colors = colors
        self.line_color = line_color
        self.label_size = label_size

    def paint(self, painter, option, widget):
        a = 5760 # 16 * 360, this is a full circle
        angle_start = 0
        radius = max(self.rect().width(), self.rect().height())/2.

        if not self.line_color:
            painter.setPen(Qt.NoPen)
        else:
            painter.setPen(QColor(self.line_color))
            painter.setPen(QColor('#000000'))


        for i, p in enumerate(self.percents):
            col = self.colors[i]
            painter.setBrush(QBrush(QColor(col)))
            angle_span = (p/100.) * a
            painter.drawPie(self.rect(), angle_start, angle_span )
            current_angle = ((angle_start + angle_span/2.) /16.)*(math.pi/180.)
            angle_start += angle_span

            if (self.labels):
                center = self.boundingRect().center()
                text = QGraphicsSimpleTextItem(self.labels[i])
                text.setParentItem(self)
                text.setPen(QColor('#000000'))
                #text.setBrush(Qt.NoBrush)

                text.setFont(QFont("Helvetica", self.label_size))
                tw = text.boundingRect().width()/2.
                th = text.boundingRect().height()/2.

                x = math.cos(current_angle) * (radius - math.sqrt(th**2 + tw**2)) + center.x() 
                y = - math.sin(current_angle) * (radius - math.sqrt(th**2 + tw**2)) + center.y() 
                # Center text according to masterItem size
                text.setPos(x, y)
                #painter.drawRect(self.boundingRect())


class LineFace(faces.Face):
    """
    Creates a Line face.
    """
    def __init__(self, width, height, fgcolor):
        faces.Face.__init__(self)
        self.width = width
        self.height = height
        self.fgcolor = fgcolor
        self.type = "item"
        self.rotable = True

    def update_items(self):
        self.item = _LineItem(self.width, self.height, self.fgcolor)

    def _width(self):
        return self.width

    def _height(self):
        return self.height


class _LineItem(QGraphicsLineItem):
    def __init__(self, w, h, fgcolor):
        QGraphicsLineItem.__init__(self)
        self.setLine(w/2., 0, w/2., h)
        if fgcolor:
            self.setPen(QPen(QColor(fgcolor)))
        else:
            self.setPen(QPen(QColor('#000000')))