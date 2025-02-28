import pya
from math import sqrt, cos, sin, tan, atan2, pi, copysign
from pya import Point, DPoint, Vector, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

import itertools

from classLib.baseClasses import ComplexBase
from classLib.shapes import Rectangle
from classLib.coplanars import CPWParameters, CPW, CPW_arc
from classLib.contactPads import ContactPad

from typing import Union, List

import copy


class FABRICATION:
    """
        Metal polygons edges are overetched by this value expressen in nm.
        Overetching results in broadening of gaps between polygons.
        In other words, every polygon edge is shifted along direction
    perpendicular to the edge itself from empty space
    to polygon's body by FABRICATION.OVERETCHING distance in nm.
        To account for overetching polygons has to be constructed
    in a way that results in software design with polygons "widened" by
    FABRICATIO.OVERETCHING value. For e.g. witdth of the coplanar
    waveguide central conductor has to be "widened" by 2*FABRICATION.OVERETCHING
    while preseving symmetry along center of the wavegiude.

        Correponding adjustments have to be made to every design element
    that undergoues overetching during fabrication.
    In addition, different areas of the sample can undergo different
    overetching, depending on the design and fabrication process.
    """
    OVERETCHING = 0.0e3


class Chip5x10_with_contactPads(ComplexBase):
    '''
    This object is implementing chip surface with 
    some default contact pads that are already present
    '''

    def __init__(self, origin, Z_params, trans_in=None):
        '''
        @params:
            origin: DPoint
                position of the left-buttom corner of the chip
            Z_params: Coplanars.CPWParameters class instance
                parameters (one value or an array of 8 values) of the coplanar waveguide used as inner end of the
                contact pads.
        '''
        self.chip_x = 10e6
        self.chip_y = 5e6
        self.center = DPoint(self.chip_x / 2, self.chip_y / 2)
        self.Z_params = [Z_params] * 8 if not isinstance(Z_params, list) else Z_params
        super().__init__(origin, trans_in)
        self.center = self.connections[-1]

    def init_primitives(self):
        origin = DPoint(0, 0)

        # drawing chip        
        self.chip = Rectangle(origin, self.chip_x, self.chip_y)
        self.primitives["chip"] = self.chip

        # contact pads
        self.contact_pad_left = ContactPad(origin + DPoint(0, self.chip_y / 2),
                                           chip_cpw_params=self.Z_params[0])
        self.primitives["cp_left"] = self.contact_pad_left

        self.contact_pad_right = ContactPad(origin + DPoint(self.chip_x, self.chip_y / 2),
                                            chip_cpw_params=self.Z_params[1],
                                            trans_in=Trans.R180)
        self.primitives["cp_right"] = self.contact_pad_right

        # top and bottom pads
        N_pads = 3
        self.connections = [self.contact_pad_left.end, self.contact_pad_right.end]
        self.angle_connections = [self.contact_pad_left.angle_connections[1],
                                  self.contact_pad_right.angle_connections[1]]

        self.contact_pads_top = [ContactPad(origin + DPoint(self.chip_x / (N_pads + 1) * (i + 1), self.chip_y),
                                            chip_cpw_params=self.Z_params[i+2],
                                            trans_in=Trans.R270) for i in range(0, N_pads)]
        for i in range(0, N_pads):
            self.primitives["cp_top_" + str(i)] = self.contact_pads_top[i]
            self.connections.append(self.contact_pads_top[i].end)
            self.angle_connections.append(self.contact_pads_top[i].angle_connections[1])

        self.contact_pads_bottom = [ContactPad(origin + DPoint(self.chip_x / (N_pads + 1) * (i + 1), 0),
                                               chip_cpw_params=self.Z_params[5+i],
                                               trans_in=Trans.R90) for i in range(0, N_pads)]
        for i in range(0, N_pads):
            self.primitives["cp_bot_" + str(i)] = self.contact_pads_bottom[i]
            self.connections.append(self.contact_pads_bottom[i].end)
            self.angle_connections.append(self.contact_pads_bottom[i].angle_connections[1])

        self.connections.append(DPoint(self.chip_x / 2, self.chip_y / 2))
        self.angle_connections.append(0)


class CHIP_10x10_12pads:
    """
    10x10 mm chip
    PCB design located here:
    https://drive.google.com/drive/folders/1TGjD5wwC28ZiLln_W8M6gFJpl6MoqZWF?usp=sharing
    """
    dx = 10e6
    dy = 10e6

    pcb_width = 260e3  # 0.26 mm
    pcb_gap = 190e3  # (0.64 - 0.26) / 2 = 0.19 mm
    pcb_feedline_d = 2500e3  # 2.5 mm
    pcb_Z = CPWParameters(pcb_width, pcb_gap)

    chip_cpw_width = 24e3
    chip_cpw_gap = 13e3
    chip_Z = CPWParameters(chip_cpw_width, chip_cpw_gap)

    @staticmethod
    def get_contact_pads(chip_Z_list: List[Union[CPWParameters, CPW, CPW_arc]]=None,
                         overetching: float =0.0e3):
        """
        Constructs objects that represent contact pads. Each pad
        consists of cpw that matches PCB cpw dimension, then trapeziod
        transition region that ends with dimensions corresponding to
        on-chip cpw.

        Parameters
        ----------
        chip_Z_list : List[Union[CPWParams, CPW, CPW_arc]]
            list of 12 structures containing dimensions of the coplanar
            waveguides on chip-side of every contact pad.
            Order starts from top-left (index 0) in counter_clockwise direction:
            top contact pad at the left side of the chip has index 0.
        overetching : float
            parameter that is used to correct contact pad's dimension
            according to fabrication process

        Returns
        -------
        list[ContactPad]
            List of contact pad objects indexed starting from top of the left corner
            in counter-clockwise direction.
        """
        if chip_Z_list is None:
            chip_Z_list = [
                CPWParameters(
                    CHIP_10x10_12pads.chip_cpw_width,
                    CHIP_10x10_12pads.chip_cpw_gap
                )
                for i in range(12)
            ]
        elif len(chip_Z_list) != 12:
            raise ValueError("`cpw_params_list` length is not equal to number of pads (12).")
        else:
            chip_Z_list = [
                CPWParameters(
                    chip_z.width + 2*overetching,
                    chip_z.gap - 2*overetching
                )
                for chip_z in chip_Z_list
            ]

        dx = CHIP_10x10_12pads.dx
        dy = CHIP_10x10_12pads.dy
        pcb_feedline_d = CHIP_10x10_12pads.pcb_feedline_d
        pcb_Z = CHIP_10x10_12pads.pcb_Z
        back_metal_gap = 100e3

        k = 0
        contact_pads_left = [
            ContactPad(
                DPoint(0, dy - pcb_feedline_d * (i + 1)),
                pcb_cpw_params=pcb_Z,
                chip_cpw_params=chip_Z_list[k + i],
                back_metal_width=50e3,
                back_metal_gap=back_metal_gap
            ) for i in range(3)
        ]
        k += 3

        contact_pads_bottom = [
            ContactPad(
                DPoint(pcb_feedline_d * (i + 1), 0), pcb_Z, chip_Z_list[k + i], back_metal_width=50e3,
                back_metal_gap=back_metal_gap,
                trans_in=Trans.R90
            ) for i in range(3)
        ]
        k += 3

        contact_pads_right = [
            ContactPad(
                DPoint(dx, pcb_feedline_d*(i+1)), pcb_Z, chip_Z_list[k + i], back_metal_width=50e3,
                back_metal_gap=back_metal_gap,
                trans_in=Trans.R180
            ) for i in range(3)
        ]
        k += 3

        contact_pads_top = [
            ContactPad(
                DPoint(dx - pcb_feedline_d * (i + 1), dy), pcb_Z, chip_Z_list[k + i], back_metal_width=50e3,
                back_metal_gap=back_metal_gap,
                trans_in=Trans.R270
            ) for i in range(3)
        ]

        # contact pads are ordered starting with top-left corner in counter-clockwise direction
        contact_pads = itertools.chain(
            contact_pads_left, contact_pads_bottom,
            contact_pads_right, contact_pads_top
        )

        return list(contact_pads)

    origin = DPoint(0, 0)
    box = pya.DBox(origin, origin + DPoint(dx, dy))

    @staticmethod
    def get_geometry_params_dict(prefix="", postfix=""):
        from collections import OrderedDict
        geometry_params = OrderedDict(
            [
                ("dx, um", CHIP_10x10_12pads.dx / 1e3),
                ("dy, um", CHIP_10x10_12pads.dy / 1e3),
                ("nX", CHIP_10x10_12pads.nX),
                ("nY", CHIP_10x10_12pads.nY)
            ]
        )
        modified_dict = OrderedDict()
        for key, val in geometry_params.items():
            modified_dict[prefix + key + postfix] = val
        return modified_dict
