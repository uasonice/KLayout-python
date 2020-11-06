import pya
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region, Box, DBox
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

import classLib
reload(classLib)
from classLib.chipDesign import ChipDesign
from classLib.shapes import Circle


def split_polygons(obj, max_pts=200, print_tree=False):
    global step_in_str
    if isinstance(obj, pya.Region):
        result_reg = Region()
        reg_to_split = obj
        for poly in reg_to_split:
            result_reg.insert(split_polygons([poly], max_pts))
        return result_reg
    elif isinstance(obj, list):
        if isinstance(obj[0], pya.Polygon):
            # if this is list of polygons
            polygons_list = obj
            resulting_polygons_list = []
            for i, poly in enumerate(polygons_list):
                if poly.num_points() < max_pts:
                    # recursion base (if all polygons satisfy this condition)
                    resulting_polygons_list.append(poly)
                else:
                    resulting_polygons_list.extend(split_polygons(poly.split()))
            return resulting_polygons_list
        else:
            raise ValueError("`split_polygons` function: List is supplied as argument, but this is not "
                             "list of `pya.Polygon` objects")
    else:
        raise ValueError("`split_polygons` function: Unknown argument received as `obj`"
                         "only `pya.Region` or `list[pya.Polygons]` are supported")


class MyDesign(ChipDesign):
    def draw(self):
        origin = DPoint(0, 0)
        circ = Circle(origin, 200e3, n_pts=400, inverse=False)
        circ.place(self.region_ph)

        self.region_ph = split_polygons(self.region_ph, max_pts=100)


import numpy as np
### MAIN FUNCTION ###
if __name__ == "__main__":
    my_design = MyDesign("testScript")
    my_design.draw()
    my_design.show()



