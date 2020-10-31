import pya
from math import pi
from pya import DPoint, DSimplePolygon, SimplePolygon
from pya import Trans, DTrans, DVector, DPath

from classLib.baseClasses import ElementBase, ComplexBase
from classLib.shapes import Circle, Kolbaska
from classLib.coplanars import CPWParameters, CPW, CPW_RL_Path

from collections import namedtuple

from typing import Union

# nice solution for parameters, man. Had really appreciated this one and put to use already (SH).
SquidParams = namedtuple("SquidParams",
                         "pad_side pad_r pads_distance p_ext_width p_ext_r sq_len sq_area j_width intermediate_width b_ext j_length n bridge")
AsymSquidParams = namedtuple("AsymSquidParams",
                             "pad_r pads_distance p_ext_width p_ext_r sq_len sq_area"
                             " j_width_1 j_width_2 intermediate_width b_ext j_length n bridge"
                             " j_length_2",
                             defaults=[None])  # btw, you're free to rename `intermediate_width`


class AsymSquid(ComplexBase):
    def __init__(self, origin, params, side=0, trans_in=None):
        '''
        Class to draw a symmetrical squid with outer positioning of the junctions.

        The notation 'length' is the dimension along the line which connects the contact pads,
        'width' is for the perpendicular direction.

        Parameters
        ----------
        origin : DPoint
            the geometrical center between centers of contact pads
        params : Union[tuple, AsymSquidParams]
        side : int
            only creates single JJ.
            `side == -1` - only left junction created
            `side == 1` - only right junction created
            `side == 0` - both junctions created (default)
        trans_in : Union[DCplxTrans, ICplxTrans]
            initial transformation

        Notes
        ------
        pad_r: float
            A radious of the circle pad.
        pads_distance:
            The distance between centers of contact pads.
        p_ext_width: float
            The width of DPath leads which connect
            contact pads and junctions.
        p_ext_r: float
            Curved extension of DPath leads along the leads direction
        sq_len: float
            The length of the squid, along leads.
        sq_area: float
            The total area of the squid.
            (does not count the reduction of area due to
             shadow angle evaporation).
        j_width_1: float
            The width of the upper small lead on left
            side (straight) and also a width of
            the junction
        j_width_2: float
            The width of the upper small lead on right
            side (straight) and also a width of
            the junction
        intermediate_width: float
            The width of the lower small bended leads
            before bending.
        b_ext: float
            The extension of bended leads after bending
        j_length: float
            The length of the jj and the width of bended parts of the lower leads.
        n: int
            The number of angle in regular polygon which serves as a large contact pad
        bridge: float
            The value of the gap between two parts of junction in the design
        j_length_2 : float
            optional
            if present, `j_length` is interpreted as
            y-dimensions of left small bottom horizontal lead,
            and `j_length_2` is interpreted as
            y-dimensions of the right small bottom bended lead
        '''
        # To draw only a half of a squid use 'side'
        # side = -1 is left, 1 is right, 0 is both (default)
        self.params: AsymSquidParams = params  # See description of AsymSquidParams tuple and the comment above
        if (self.params.intermediate_width < self.params.j_width_1) or \
                (self.params.intermediate_width < self.params.j_width_2):
            raise ValueError("AsymSquid constructor:\n"
                             "intermediate lead width assumed be bigger than width of each jj")
        if self.params.j_length_2 is None:
            # workaround because namedtuple is immutable
            # there is `recordtype` that is mutable, but
            # in is not included in default KLayout build.
            self.params = AsymSquidParams(*(self.params[:-1] + (self.params.j_length,)))
        self.side = side
        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        pars = self.params
        up_pad_center = origin + DVector(0, pars.pads_distance / 2)
        down_pad_center = origin + DVector(0, -pars.pads_distance / 2)
        self.pad_down = Circle(
            down_pad_center, pars.pad_r,
            n_pts=pars.n, offset_angle=pi / 2
        )
        self.primitives["pad_down"] = self.pad_down

        self.p_ext_down = Kolbaska(
            down_pad_center, origin + DPoint(0, -pars.sq_len / 2),
            pars.p_ext_width, pars.p_ext_r
        )
        self.primitives["p_ext_down"] = self.p_ext_down

        self.pad_up = Circle(
            up_pad_center, pars.pad_r,
            n_pts=pars.n, offset_angle=-pi / 2
        )
        self.primitives["pad_up"] = self.pad_up

        self.p_ext_up = Kolbaska(
            up_pad_center, origin + DVector(0, pars.sq_len / 2),
            pars.p_ext_width, pars.p_ext_r
        )
        self.primitives["p_ext_up"] = self.p_ext_up

        origin = DPoint(0, 0)
        if self.side == 0:
            self.init_half(origin, side=1)  # right
            self.init_half(origin, side=-1)  # left
        else:
            self.init_half(origin, side=self.side)

    def init_half(self, origin, side=-1):
        # side = -1 is left, 1 is right
        pars = self.params
        j_width = pars.j_width_1 if side < 0 else pars.j_width_2
        j_length = pars.j_length if side < 0 else pars.j_length_2
        suff = "_left" if side < 0 else "_right"
        up_st_gap = pars.sq_area / (2 * pars.sq_len)

        # exact correction in first row
        # additional extension to isolate jj's from intermediate bottom polygons
        # without correction horizontal faces of jj's will be adjacent
        # to thick intemediate polygons
        low_st_gap = up_st_gap + ((pars.j_length + pars.j_length_2) / 2 + pars.intermediate_width) + \
                     2 * pars.intermediate_width

        ### upper and lower vertical intermediate and jj leads ###
        ## top leads ##
        upper_leads_extension = j_width / 4
        # upper intemediate lead
        up_st_start = self.primitives["p_ext_up"].connections[1] + \
                      DVector(side * up_st_gap / 2, 0)
        up_st_stop = origin + \
                     DVector(side * up_st_gap / 2, pars.bridge / 2 + upper_leads_extension)
        # top jj lead
        self.primitives["upp_st" + suff] = Kolbaska(
            up_st_start, up_st_stop, j_width, upper_leads_extension
        )
        # top intermediate lead
        upper_thin_part_len = 4 * pars.bridge + pars.intermediate_width / 2
        self.primitives["upp_st_thick" + suff] = Kolbaska(
            up_st_start, up_st_stop + DPoint(0, upper_thin_part_len),
            pars.intermediate_width, pars.intermediate_width / 2
        )

        ## bottom leads ##
        low_st_start = self.primitives["p_ext_down"].connections[1] + \
                       DVector(side * low_st_gap / 2, 0)
        low_st_stop = origin + \
                      DVector(side * (low_st_gap / 2 + 2 * pars.intermediate_width),
                              -pars.bridge / 2 - pars.intermediate_width / 2)
        len_ly = (low_st_stop - low_st_start).y
        # bottom intermediate lead
        self.primitives["low_st" + suff] = CPW_RL_Path(
            low_st_start, 'LR', CPWParameters(pars.intermediate_width, 0),
            pars.intermediate_width / 2, [len_ly], [side * pi / 2], trans_in=DTrans.R90)
        # bottom jj lead (horizontal)
        low_st_end = self.primitives["low_st" + suff].connections[1]
        low_st_jj_start = low_st_end + DPoint(0, -j_length / 2 + pars.intermediate_width / 2)
        low_st_jj_stop = low_st_jj_start + DPoint(-side * pars.b_ext, 0)
        self.primitives["low_st_jj" + suff] = Kolbaska(low_st_jj_start, low_st_jj_stop, j_length,
                                                       j_length / 2)


class Squid(AsymSquid):
    '''
    Class to draw a symmetrical squid with outer positioning of the junctions.

    The notation 'length' is the dimension along the line which connects the contact pads,
    'width' is for the perpendicular direction.

    Parameters:
        pad_side: float
            A length of the side of triangle pad.
        (??) pad_r: float
            The outer_r of round angle of the contact pad.
        pads_distance:
            The distance between triangle contact pads.
        p_ext_width: float
            The width of curved rectangle leads which connect triangle contact pads and junctions.
        p_ext_r: float
            The angle outer_r of the pad extension
        sq_len: float
            The length of the squid, along leads.
        sq_width: float
            The total area of the squid.
            (does not count the reduction of area due to shadow angle evaporation).
        j_width: float
            The width of the upper small leads (straight) and also a width of the junction
        intermediate_width: float
            The width of the lower small bended leads before bending
        b_ext: float
            The extension of bended leads after bending
        j_length: float
            The length of the jj and the width of bended parts of the lower leads.
        n: int
            The number of angle in regular polygon which serves as a large contact pad
        bridge: float
            The value of the gap between two parts of junction in the design
        trans_in: Trans
            Initial transformation
    '''

    def __init__(self, origin, params, side=0, trans_in=None):
        # To draw only a half of a squid use 'side'
        # side = -1 is left, 1 is right, 0 is both (default)
        asymparams = AsymSquidParams(*params[:-3], *params[-2:])
        super().__init__(self, origin, asymparams, side, trans_in)


class Line_N_JJCross(ElementBase):
    def __init__(self, origin, params, trans_in=None):
        self.params = params
        self.a = params[0]
        self.b = params[1]
        self.jos1_b = params[2]
        self.jos1_a = params[3]
        self.f1 = params[4]
        self.d1 = params[5]
        self.jos2_b = params[6]
        self.jos2_a = params[7]
        self.f2 = params[8]
        self.d2 = params[9]
        self.w = params[10]

        self.poly1 = self._make_polygon(self.b, self.w, self.d1, self.f1, self.d2)

        super().__init__(origin, trans_in)

    def _make_polygon(self, length, w, d, f, overlapping):
        polygon = DSimplePolygon
        p1 = DPoint(0, 0)
        p2 = p1 + DPoint(length, 0)
        p3 = p2 + DPoint(0, w)
        p4 = p3 - DPoint(overlapping, 0)
        p5 = p4 - DPoint(0, d)
        p6 = p5 - DPoint(f, 0)
        p7 = p6 + DPoint(0, d)
        p8 = p1 + DPoint(0, w)

        polygon = DSimplePolygon([p1, p2, p3, p4, p5, p6, p7, p8])
        return polygon

    def init_regions(self):
        self.metal_region.insert(SimplePolygon().from_dpoly(self.poly1))
