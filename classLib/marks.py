import pya
from math import sqrt, cos, sin, atan2, pi, copysign
from pya import Point, DPoint, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Vector, DVector
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from classLib.baseClasses import ComplexBase
from classLib.shapes import *
from classLib.coplanars import CPWParameters, CPW_arc


class Mark1(ComplexBase):
    def __init__(self, origin, trans_in=None):
        self.cross_in_a = 20e3
        self.cross_out_a = 40e3
        self.leaf_inner = 50e3
        self.leaf_outer = 150e3
        self.leaf_angle = 70
        self.leafs_N = 4
        self.empty_rings_N = 2
        self.empty_rings_width = 10e3

        self.empty_leaf_angle = (360 - self.leaf_angle * self.leafs_N) / self.leafs_N
        self.avg_r = (self.leaf_inner + self.leaf_outer) / 2

        if (self.cross_out_a < self.cross_in_a):
            print("cross inner square must be larger than outer")
        super(self).__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)

        # square box that clears out the area of the mark
        self.primitives["empty_box"] = Rectangle(DPoint(-self.leaf_outer, -self.leaf_outer), 2 * self.leaf_outer,
                                                 2 * self.leaf_outer, inverse=True)

        self.primitives["cross"] = Cross(DPoint(-self.cross_out_a / 2, -self.cross_out_a / 2), self.cross_in_a,
                                         self.cross_out_a)

        Z = CPWParameters(self.leaf_outer - self.leaf_inner, 0)

        for i in range(self.leafs_N):
            start_angle = i * (self.leaf_angle + self.empty_leaf_angle) - self.leaf_angle / 2
            trans = DCplxTrans(1, start_angle + 90, False, 0, -self.avg_r)
            self.primitives["leaf_" + str(i)] = CPW_arc(Z, origin, self.avg_r, self.leaf_angle * pi / 180,
                                                        trans_in=trans)

        Z_empty = CPWParameters(0, self.empty_rings_width / 2)
        for i in range(1, self.empty_rings_N + 1):
            r = self.leaf_inner + (self.leaf_outer - self.leaf_inner) / (self.empty_rings_N + 1) * i
            self.primitives["empty_ring_" + str(i)] = CPW_arc(Z_empty, DVector(0, -r), r, 2 * pi)


class Mark2(ComplexBase):
    def __init__(self, origin, cross_thickness=1e3, cross_size=3e3, trans_in=None):
        self.ring1_outer_r = 300e3
        self.ring1_thickness = 30e3
        self.ring2_outer_r = 250e3
        self.ring2_thickness = 30e3
        self.inner_circle_radius = 200e3
        self.trap_h = 150e3
        self.trap_t = 100e3
        self.cross_thickness = cross_thickness
        self.cross_size = cross_size
        self.trap_b = self.cross_thickness
        self.trap_dist = 2*self.cross_thickness  # distance from the center
        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        self.primitives["empty_ring1"] = Ring(origin, self.ring1_outer_r, self.ring1_thickness, inverse=True)
        self.primitives["empty_ring2"] = Ring(origin, self.ring2_outer_r, self.ring2_thickness, inverse=True)
        self.primitives["inner_circle"] = Circle(origin, self.inner_circle_radius, inverse=True)
        self.primitives["trap_top"] = IsoTrapezoid(
            origin + DPoint(-self.trap_b / 2, self.trap_dist),
            self.trap_h, self.trap_b, self.trap_t
        )
        self.primitives["trap_left"] = IsoTrapezoid(
            origin + DPoint(-self.trap_dist, -self.trap_b / 2),
            self.trap_h, self.trap_b, self.trap_t, trans_in=Trans.R90
        )
        self.primitives["trap_bottom"] = IsoTrapezoid(
            origin + DPoint(self.trap_b / 2, -self.trap_dist),
            self.trap_h, self.trap_b, self.trap_t, trans_in=Trans.R180
        )
        self.primitives["trap_right"] = IsoTrapezoid(
            origin + DPoint(self.trap_dist, self.trap_b / 2),
            self.trap_h, self.trap_b, self.trap_t, trans_in=Trans.R270
        )
        self.primitives["cross"] = Cross2(origin, self.cross_thickness, self.cross_size)
