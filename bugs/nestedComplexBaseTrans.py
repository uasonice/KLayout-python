import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib

reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Circle
from classLib.baseClasses import ComplexBase, ElementBase


class myBoxElement(ElementBase):
    def __init__(self, origin, a, trans_in=None):
        self.a = a
        super().__init__(origin, trans_in)

    def init_regions(self):
        origin = DPoint(0, 0)
        p1 = origin + DPoint(-self.a / 2, -self.a / 2)
        p2 = origin + DPoint(self.a / 2, self.a / 2)
        self.metal_region.insert(pya.Box().from_dbox(pya.DBox(p1, p2)))


class myBoxComplex(ComplexBase):
    def __init__(self, origin, b, trans_in=None):
        self.b = b
        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        self.box = myBoxElement(origin, self.b)
        self.primitives["box"] = self.box


class A(ComplexBase):
    def __init__(self, origin, a, b, trans_in=None):
        self.a = a
        self.b = b
        super().__init__(origin, trans_in)

    def init_primitives(self):
        origin = DPoint(0, 0)
        self.boxElement = myBoxElement(origin, self.a)
        self.primitives["boxElement"] = self.boxElement

        self.boxComplex = myBoxComplex(origin + DPoint(self.a + self.b, self.a + self.b), self.b)
        self.primitives["boxComplex"] = self.boxComplex


class MyDesign(ChipDesign):
    def __init__(self, name="testScript"):
        super().__init__(name)
        self.a_class: A = None

    def draw(self):
        origin = DPoint(100e3, 100e3)
        self.a_class = A(origin, 10e3, 2e3)
        self.a_class.place(self.region_ph)


import numpy as np

### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()
    a_class = my_design.a_class
    my_design.cell.shapes(my_design.layer_el).insert(a_class.metal_region)
