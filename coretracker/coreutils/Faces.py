import math

from PyQt4.QtCore import Qt, QPointF
from PyQt4.QtGui import (QGraphicsRectItem, QGraphicsLineItem, QPen,
                         QColor, QBrush, QFont, QGraphicsSimpleTextItem, QLinearGradient)
from ete3 import faces
from ete3.treeview.main import COLOR_SCHEMES


_aafgcolors = {
    'A': "#000000",
    'R': "#000000",
    'N': "#000000",
    'D': "#000000",
    'C': "#000000",
    'Q': "#000000",
    'E': "#000000",
    'G': "#000000",
    'H': "#000000",
    'I': "#000000",
    'L': "#000000",
    'K': "#000000",
    'M': "#000000",
    'F': "#000000",
    'P': "#000000",
    'S': "#000000",
    'T': "#000000",
    'W': "#000000",
    'Y': "#000000",
    'V': "#000000",
    'B': "#000000",
    'Z': "#000000",
    'X': "#000000",
    '.': "#000000",
    '-': "#000000",
    '*': "#000000",
}

_aabgcolors = {
    'A': "#C8C8C8",
    'R': "#145AFF",
    'N': "#00DCDC",
    'D': "#E60A0A",
    'C': "#E6E600",
    'Q': "#00DCDC",
    'E': "#E60A0A",
    'G': "#DBDBDB",
    'H': "#8282D2",
    'I': "#0F820F",
    'L': "#0F820F",
    'K': "#145AFF",
    'M': "#E6E600",
    'F': "#3232AA",
    'P': "#DC9682",
    'S': "#FA9600",
    'T': "#FA9600",
    'W': "#B45AB4",
    'Y': "#3232AA",
    'V': "#0F820F",
    'B': "#FF69B4",
    'Z': "#FF69B4",
    'X': "#BEA06E",
    '.': "#FFFFFF",
    '-': "#FFFFFF",
    '*': "#FFFFFF",
}

_ntfgcolors = {
    'A': '#000000',
    'G': '#000000',
    'I': '#000000',
    'C': '#000000',
    'T': '#000000',
    'U': '#000000',
    '.': "#000000",
    '-': "#000000",
    ' ': "#000000"
}

_ntbgcolors = {
    'A': '#A0A0FF',
    'G': '#FF7070',
    'I': '#80FFFF',
    'C': '#FF8C4B',
    'T': '#A0FFA0',
    'U': '#FF8080',
    '.': "#FFFFFF",
    '-': "#FFFFFF",
    ' ': "#FFFFFF"
}


def _get_codon_fgcolors(codontable, cible_aa):
    """Get colon foreground color"""
    return dict((k, (_aafgcolors[v] if v != cible_aa else '#FFFFFF')) for (k, v) in codontable.items())


def _get_codon_bgcolors(codontable, cible_aa, spec_codon_col):
    """Get colon background color"""
    return dict((k, spec_codon_col.get(k, ("#FFFFFF" if v != cible_aa else '#000000'))) for (k, v) in codontable.items())


class PPieChartFace(faces.StaticItemFace):
    """
    .. versionadded:: 2.2
    :param percents: a list of values summing up 100.
    :param width: width of the piechart
    :param height: height of the piechart
    :param colors: a list of colors (same length as percents)
    :param line_color: color used to render the border of the piechart (None=transparent)
    """

    def __init__(self, percents, width, height, colors=None, line_color=None, label_size=6, show_label=False, is_percent=False):
        faces.Face.__init__(self)
        self.labels = None
        if show_label:
            self.labels = sum(percents)

        if not is_percent:
            s = sum(percents)
            percents = map(lambda x: x * 100. / s, percents)

        if round(sum(percents)) > 100:
            raise ValueError("PPieChartItem: percentage values > 100")

        self.type = "item"
        self.item = None
        self.percents = percents
        if not colors:
            colors = COLOR_SCHEMES["paired"]
        self.colors = colors
        self.width = width
        self.height = height
        self.label_size = label_size
        self.line_color = line_color

    def update_items(self):
        self.item = _PieChartItem(self.percents, self.width,
                                  self.height, self.colors, self.line_color)
        self.add_text()

    def _width(self):
        return self.item.rect().width()

    def _height(self):
        return self.item.rect().height()

    def add_text(self):
        if (self.labels):
            center = self.item.boundingRect().center()
            text = QGraphicsSimpleTextItem(str(self.labels))
            text.setFont(QFont("Arial", self.label_size))
            text.setParentItem(self.item)
            text.setBrush(QBrush(QColor('#ddd')))
            tw = text.boundingRect().width() / 2.
            th = text.boundingRect().height() / 2.
            x = -tw + center.x()
            y = -th + center.y()
            # Center text according to masterItem size
            text.setPos(x, y)


class _PieChartItem(QGraphicsRectItem):

    def __init__(self, percents, width, height, colors, line_color=None):
        QGraphicsRectItem.__init__(self, 0, 0, width, height)
        self.percents = percents
        self.colors = colors
        self.line_color = line_color

    def paint(self, painter, option, widget):
        a = 5760  # 16 * 360, this is a full circle
        angle_start = 0
        # radius = max(self.rect().width(), self.rect().height()) / 2.

        if not self.line_color:
            painter.setPen(Qt.NoPen)
        else:
            painter.setPen(QColor(self.line_color))
            painter.setPen(QColor('#000000'))

        for i, p in enumerate(self.percents):
            col = self.colors[i]
            painter.setBrush(QBrush(QColor(col)))
            angle_span = (p / 100.) * a
            painter.drawPie(self.rect(), angle_start, angle_span)
            current_angle = ((angle_start + angle_span / 2.) /
                             16.) * (math.pi / 180.)
            angle_start += angle_span
            # find middle radius of the arc span
            # add a new line from the intersection point with the arc
            # add a second line from the end of that line in a
            # check if vertical then add a new vertical path
            # check angle orientation before selecting the orientation of the second line
            # probably do this with a path object and multiple point
            # then add a text object


class _RingChartItem(_PieChartItem):

    def __init__(self, percents, width, height, colors, line_color=None):
        QGraphicsRectItem.__init__(self, 0, 0, width, height)
        self.percents = percents
        self.colors = colors
        self.line_color = line_color

    def paint(self, painter, option, widget):
        a = 5760  # 16 * 360, this is a full circle
        angle_start = 0
        # radius = max(self.rect().width(), self.rect().height()) / 2.

        if not self.line_color:
            painter.setPen(Qt.NoPen)
        else:
            painter.setPen(QColor(self.line_color))
            painter.setPen(QColor('#000000'))

        for i, p in enumerate(self.percents):
            col = self.colors[i]
            painter.setBrush(QBrush(QColor(col)))
            angle_span = (p / 100.) * a
            painter.drawPie(self.rect(), angle_start, angle_span)
            current_angle = ((angle_start + angle_span / 2.) /
                             16.) * (math.pi / 180.)
            angle_start += angle_span


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
        self.setLine(w / 2., 0, w / 2., h)
        if fgcolor:
            self.setPen(QPen(QColor(fgcolor)))
        else:
            self.setPen(QPen(QColor('#000000')))


class SequenceFace(faces.StaticItemFace):
    """
    Creates a new molecular sequence face object.
    :param seq: Sequence string to be drawn
    :param seqtype: Type of sequence: "nt" or "aa"
    :param fsize: Font size, (default=10)
    You can set custom colors for amino-acids or nucleotides:
    :param None codon: a string that corresponds to the reverse
      translation of the amino-acid sequence
    :param None col_w: width of the column (if col_w is lower than
      font size, letter wont be displayed)
    :param None fg_colors: dictionary of colors for foreground, with
      as keys each possible character in sequences, and as value the
      colors
    :param None bg_colors: dictionary of colors for background, with
      as keys each possible character in sequences, and as value the
      colors
    :param 3 alt_col_w: works together with special_col option,
      defines the width of given columns
    :param None special_col: list of lists containing the bounds
      of columns to be displayed with alt_col_w as width
      """

    def __init__(self, seq, cible_aa, seqtype="aa", fsize=10,
                 fg_colors=None, bg_colors=None, codon=None,
                 col_w=None, alt_col_w=3, special_col=None,
                 spec_codon_col=None, codontable={}):
        self.seq = seq
        self.codon = codon
        self.fsize = fsize
        self.style = seqtype
        self.col_w = float(self.fsize + 1) if col_w is None else float(col_w)
        self.alt_col_w = float(alt_col_w)
        self.special_col = special_col if special_col else []
        self.width = 0  # will store the width of the whole sequence

        if self.style == "aa":
            if not fg_colors:
                fg_colors = _aafgcolors
            if not bg_colors:
                bg_colors = _aabgcolors

        elif self.style == 'codon':
            self.col_w *= 3
            if not isinstance(self.seq, list):
                # only consider the position where 3 nuc can be obtained
                self.seq = [self.seq[i:i + 3] for i in xrange(0,
                                                              len(self.seq) - len(self.seq) % 3, 3)]
            if not fg_colors:
                fg_colors = _get_codon_fgcolors(codontable, cible_aa)
            if not bg_colors:
                bg_colors = _get_codon_bgcolors(
                    codontable, cible_aa, spec_codon_col)

        else:
            if not fg_colors:
                fg_colors = _ntfgcolors
            if not bg_colors:
                bg_colors = _ntbgcolors

        def __init_col(color_dic, keys, col='#000000'):
            """to speed up the drawing of colored rectangles and characters"""
            new_color_dic = {}
            for car in keys:
                # use default black color if color not found
                new_color_dic[car] = QBrush(QColor(color_dic.get(car, col)))
            return new_color_dic

        self.fg_col = __init_col(fg_colors, set(self.seq))
        self.bg_col = __init_col(bg_colors, set(self.seq), '#FFFFFF')

        # for future?
        self.row_h = 13.0
        super(SequenceFace, self).__init__(None)

    def update_items(self):
        rect_cls = QGraphicsRectItem
        self.item = rect_cls(0, 0, self.width, self.row_h)
        seq_width = 0
        nopen = QPen(Qt.NoPen)
        font = QFont("Courier", self.fsize)
        for i, letter in enumerate(self.seq):
            width = self.col_w
            for reg in self.special_col:
                if reg[0] < i <= reg[1]:
                    width = self.alt_col_w
                    break
            rectitem = rect_cls(0, 0, width, self.row_h, parent=self.item)
            rectitem.setX(seq_width)  # to give correct X to children item
            rectitem.setBrush(self.bg_col[letter])
            rectitem.setPen(nopen)

            # write letter if enough space
            if width >= self.fsize:
                text = QGraphicsSimpleTextItem(letter, parent=rectitem)
                text.setFont(font)
                text.setBrush(self.fg_col[letter])
                # Center text according to rectitem size
                txtw = text.boundingRect().width()
                txth = text.boundingRect().height()
                text.setPos((width - txtw) / 2, (self.row_h - txth) / 2)
            seq_width += width
        self.width = seq_width


class List90Face(faces.StaticItemFace):
    """Static text Face object
    :param l:        List of element to be drawn
    :param fsize:    Font size, e.g. 10,12,6, (default=10)
    :param fgcolor:  Foreground font color. RGB code or color name in :data:`SVG_COLORS`
    :param bgcolor:  Background font color. RGB code or color name in :data:`SVG_COLORS`
    """

    def __init__(self, l, ftype="Courier", fstyle="normal", fsize=10,
                 fgcolor="black", bgcolor="white", col_w=14.0, rotation=90):
        self.liste = l
        self.ftype = ftype
        self.fgcolor = fgcolor
        self.bgcolor = bgcolor
        self.fsize = fsize
        self.row_h = float(self.fsize + 1)
        self.col_w = col_w
        self.width = 0
        self.rot = rotation
        self.fstyle = fstyle
        self.coeff_h = max([len(str(x)) for x in self.liste])

        super(List90Face, self).__init__(None)

    def __repr__(self):
        return "Text Face [%s] (%s)" % (self._text, hex(self.__hash__()))

    def get_text(self):
        return self._text

    def update_items(self):
        self.item = QGraphicsRectItem(
            0, 0, self.width, self.row_h * self.coeff_h)
        seq_width = 0
        nopen = QPen(Qt.NoPen)
        self.item.setPen(nopen)
        font = QFont(self.ftype, self.fsize)
        if self.fstyle == "italic":
            font.setStyle(QFont.StyleItalic)
        elif self.fstyle == "oblique":
            font.setStyle(QFont.StyleOblique)
        rect_cls = QGraphicsRectItem
        for i, val in enumerate(self.liste):
            width = self.col_w
            height = self.row_h * len(str(val)) + 1
            rectitem = rect_cls(0, 0, width, height, parent=self.item)
            rectitem.setX(seq_width)  # to give correct X to children item
            rectitem.setBrush(QBrush(QColor(self.bgcolor)))
            rectitem.setPen(nopen)

            # write letter if enough space in height
            if height >= self.fsize:
                text = QGraphicsSimpleTextItem(str(val), parent=rectitem)
                text.setFont(font)
                text.setBrush(QBrush(QColor(self.fgcolor)))
                # Center text according to rectitem size
                # txtw = text.boundingRect().width()
                txth = text.boundingRect().height()
                text.setRotation(self.rot)
                text.setX(txth)
            seq_width += width
        self.width = seq_width


class ReaRectFace(faces.StaticItemFace):
    """Create a list of rect face for each node
    in order to plot reassignment in both extant and
    ancestral sequence
    """

    def __init__(self, aalist, readict, is_leaf=True, spacer=1, height=12,
                 ffamily='Courier', ncodons=None, fsize=12, col_w=None, margin_left=0):
        self.readict = readict
        if self.readict is None:
            raise ValueError("Readict is required to not be None")

        faces.StaticItemFace.__init__(self, None)
        self.spacer = spacer
        self.h = height
        self.maxcodon = ncodons
        self.ffamily = ffamily
        self.fsize = fsize
        self.w = float(self.fsize + 1) if col_w is None else float(col_w)
        self.w *= 3
        self.width = 0
        self.aa_list = sorted(aalist)
        self.mgl = margin_left
        self.is_leaf = is_leaf

        def _init_colors(bgtype=True, fgcolor='#000000'):
            color = {}
            for aa in readict.keys():
                c = _aabgcolors[aa.upper()] if bgtype else fgcolor
                color[aa.upper()] = QBrush(QColor(c))
            return color
        self.fgcolor = _init_colors(False)
        self.bgcolor = _init_colors()

    def update_items(self):
        try:
            max_codons = math.ceil(
                max([len(x) for x in self.readict.values()]) / 2.0) * 2
        except:
            max_codons = 1
        if self.maxcodon:
            max_codons = max(max_codons, self.maxcodon)

        max_h = max_codons * (self.h + self.spacer)
        rect_cls = QGraphicsRectItem
        self.item = rect_cls()
        nopen = QPen(QColor('#EEEEEE'))
        nobrush = QBrush(Qt.NoBrush)
        width = self.mgl
        font = QFont(self.ffamily, self.fsize)
        for aa in self.aa_list:
            codons = self.readict.get(aa, [])
            tot_codons = len(codons)
            hpos = (self.h + self.spacer) * (max_codons - tot_codons) / 2.0

            for cod in codons:
                rectitem = rect_cls(0, 0, self.w, self.h, parent=self.item)
                rectitem.setX(width)
                rectitem.setY(hpos)
                rectitem.setBrush(self.bgcolor[aa])
                rectitem.setPen(nopen)
                hpos += (self.h + self.spacer)
                # write letter if enough space
                if self.w >= self.fsize:
                    text = QGraphicsSimpleTextItem(cod, parent=rectitem)
                    text.setFont(font)
                    text.setBrush(self.fgcolor[aa])
                    # Center text according to rectitem size
                    txtw = text.boundingRect().width()
                    txth = text.boundingRect().height()
                    text.setPos((self.w - txtw) / 2, (self.h - txth) / 2)

            # this only happen if codon reassignment not found for aa
            # we do not need a spacer if it's an internal node (I hope)
            if hpos == 0 and self.is_leaf:
                rectitem = rect_cls(0, 0, self.w, self.h, parent=self.item)
                rectitem.setX(width)
                rectitem.setY(hpos)
                rectitem.setBrush(nobrush)
                rectitem.setPen(nopen)

            width += self.w + self.spacer

        self.width = width
        self.item.setPen(nopen)
        self.item.setRect(0, 0, self.width, max_h)


class SummaryRectFace(faces.StaticItemFace):
    """Rectangular faces that contains genome count for each  a list of rect face for each node
    in order to plot reassignment in both extant and
    ancestral sequence
    """

    def __init__(self, count, aalist, height=12, width=100, margin=1, colormap={}, ffamily='Courier', fsize=12, fgcolor='#000000'):

        faces.StaticItemFace.__init__(self, None)
        self.h = height
        self.w = width
        self.margin = margin

        self.font = QFont(ffamily, fsize)
        self.aa_list = aalist
        self.codons = count
        self.width = margin

        self.fgcolor = QColor(fgcolor)
        colormap = colormap if colormap else _aabgcolors
        self.bgcolor = {}
        for cod, aas in self.aa_list.items():
            for aa in aas:
                aa = aa.upper()
                if aa not in self.bgcolor:
                    self.bgcolor[aa] = QColor(
                        colormap.get(aa, _aabgcolors[aa]))

    def update_items(self):
        rect_cls = QGraphicsRectItem
        nobrush = QBrush(Qt.NoBrush)
        nopen = QPen(QColor('#FFFFFF'))
        grad = Qt.NoBrush
        self.item = rect_cls()

        for codon in self.codons.keys():
            rectitem = rect_cls(self.width, self.margin,
                                self.w, self.h, parent=self.item)
            total_rea = len(self.aa_list[codon])
            if total_rea > 0:

                grad = QLinearGradient(QPointF(0, 0.5), QPointF(1, 0.5))
                grad.setCoordinateMode(QLinearGradient.ObjectBoundingMode)
                pointpos = 1.0 / total_rea
                starting = -0.001
                for reqcol, aa in enumerate(self.aa_list[codon]):
                    curcol = self.bgcolor[aa]
                    # use the same color twice to mark start and end
                    grad.setColorAt(starting + 0.001, curcol)
                    starting += pointpos
                    grad.setColorAt(starting, curcol)

                    # grad.setColorAt(starting, QColor(curcol))
            # put small rec in big rec
            # Comment award of the year !
            brush = QBrush(QColor('#CCCCCC'))
            pen = QPen(QColor('#BBBBBB'))
            if self.codons[codon]:
                brush = QBrush(grad)
                pen = QPen(QColor('#000000'))

            rectitem.setBrush(brush)
            rectitem.setPen(pen)
            # Center text according to rectitem size
            text = QGraphicsSimpleTextItem(self.codons[codon], parent=rectitem)
            text.setFont(self.font)
            text.setBrush(self.fgcolor)
            center = rectitem.boundingRect().center()
            txtw = text.boundingRect().width()
            txth = text.boundingRect().height()
            text.setPos(center.x() - txtw / 2, center.y() - txth / 2)

            self.width += self.w + self.margin

        # now plot the big rect (with margin)
        self.item.setPen(nopen)
        self.item.setBrush(nobrush)
        self.item.setRect(0, 0, self.width, self.h + (2 * self.margin))
