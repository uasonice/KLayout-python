import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib
reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Circle


class MyDesign(ChipDesign):
    def draw(self):
        origin = DPoint(0, 0)
        circ = Circle(origin, 200e3, n_pts=2000, inverse=False)
        circ.place(self.region_ph)


import numpy as np
### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()



