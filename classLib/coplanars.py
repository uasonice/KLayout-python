import pya
from math import sqrt, cos, sin, atan2, pi, copysign, tan
from numpy import sign
from pya import Point, DPoint, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from typing import Union, List
from collections import OrderedDict
import itertools
import copy

from classLib.baseClasses import ElementBase, ComplexBase
from classLib.bridgedCoplanars import BridgedCPW, BridgedCPWArc


class CPWParameters(ElementBase):
    def __init__(self, width, gap):
        self.width = width
        self.gap = gap
        self.b = 2 * gap + width


class CPW(ElementBase):
    """@brief: class represents single coplanar waveguide
        @params:  float width - represents width of the central conductor
                        float gap - spacing between central conductor and ground planes
                        float gndWidth - width of ground plane to be drawed
                        DPoint start - center aligned point, determines the start point of the coplanar segment
                        DPoint end - center aligned point, determines the end point of the coplanar segment
    """

    def __init__(self, width=None, gap=None, start=DPoint(0, 0), end=DPoint(0, 0), gndWidth=-1, trans_in=None,
                 cpw_params=None):
        if (cpw_params is None):
            self.width = width
            self.gap = gap
            self.b = 2 * gap + width
        else:
            self.width = cpw_params.width
            self.gap = cpw_params.gap
            self.b = 2 * self.gap + self.width
        self.gndWidth = gndWidth
        self.end = end
        self.start = start
        self.dr = end - start
        super().__init__(start, trans_in)

        self._geometry_parameters = OrderedDict(
            [
                ("width, um", self.width / 1e3),
                ("gap, um", self.gap / 1e3),
                ("start.x, um", self.start.x / 1e3),
                ("start.y, um", self.start.y / 1e3),
                ("end.x", self.end.x / 1e3),
                ("end.y", self.end.y / 1e3)
            ]
        )

    def init_regions(self):
        self.connections = [DPoint(0, 0), self.dr]
        self.start = DPoint(0, 0)
        self.end = self.start + self.dr
        alpha = atan2(self.dr.y, self.dr.x)
        self.angle_connections = [alpha, alpha]
        alpha_trans = ICplxTrans().from_dtrans(DCplxTrans(1, alpha * 180 / pi, False, self.start))
        metal_poly = DSimplePolygon([DPoint(0, -self.width / 2),
                                     DPoint(self.dr.abs(), -self.width / 2),
                                     DPoint(self.dr.abs(), self.width / 2),
                                     DPoint(0, self.width / 2)])
        self.connection_edges = [3, 1]
        self.metal_region.insert(pya.SimplePolygon().from_dpoly(metal_poly))
        if (self.gap != 0): self.empty_region.insert(pya.Box(Point().from_dpoint(DPoint(0, self.width / 2)),
                                                             Point().from_dpoint(
                                                                 DPoint(self.dr.abs(), self.width / 2 + self.gap))))
        self.empty_region.insert(pya.Box(Point().from_dpoint(DPoint(0, -self.width / 2 - self.gap)),
                                         Point().from_dpoint(DPoint(self.dr.abs(), -self.width / 2))))
        self.metal_region.transform(alpha_trans)
        self.empty_region.transform(alpha_trans)

    def _refresh_named_connections(self):
        self.end = self.connections[1]
        self.start = self.connections[0]
        self.dr = self.end - self.start

    def _refresh_named_angles(self):
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]


class CPW_arc(ElementBase):
    def __init__(self, Z0, start, R, delta_alpha, gndWidth=-1, trans_in=None):
        self.R = R
        self.start = start
        self.center = start + DPoint(0, self.R)
        self.end = self.center + DPoint(sin(delta_alpha), -cos(delta_alpha)) * self.R
        self.dr = self.end - self.start

        self.width = Z0.width
        self.gap = Z0.gap
        self.b = self.width + 2*self.gap

        self.delta_alpha = delta_alpha
        self.alpha_start = 0
        self.alpha_end = self.delta_alpha

        super().__init__(start, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.center = self.connections[2]

        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def _get_solid_arc(self, center, R, width, alpha_start, alpha_end, n_inner, n_outer):
        pts = []
        #        print(alpha_start/pi, alpha_end/pi, cos( alpha_start ), cos( alpha_end ),
        #                         sin(alpha_start), sin(alpha_end))

        if alpha_end > alpha_start:
            alpha_start = alpha_start - 1e-3
            alpha_end = alpha_end + 1e-3
        else:
            alpha_start = alpha_start + 1e-3
            alpha_end = alpha_end - 1e-3

        d_alpha_inner = (alpha_end - alpha_start) / (n_inner - 1)
        d_alpha_outer = -d_alpha_inner

        #        print("Center:", center)
        for i in range(0, n_inner):
            alpha = alpha_start + d_alpha_inner * i
            pts.append(center + DPoint(cos(alpha), sin(alpha)) * (R - width / 2))
        for i in range(0, n_outer):
            alpha = alpha_end + d_alpha_outer * i
            pts.append(center + DPoint(cos(alpha), sin(alpha)) * (R + width / 2))
        #        print("Points:", pts[:n_inner],"\n       ", pts[n_inner:], "\n")
        return DSimplePolygon(pts)

    def init_regions(self):
        self.connections = [DPoint(0, 0), self.dr, DPoint(0, self.R)]
        self.angle_connections = [self.alpha_start, self.alpha_end]
        self.start = DPoint(0, 0)
        self.end = self.dr
        self.center = DPoint(0, self.R)

        from ._PROG_SETTINGS import PROGRAM
        n_inner = PROGRAM.ARC_PTS_N
        n_outer = PROGRAM.ARC_PTS_N

        metal_arc = self._get_solid_arc(self.center, self.R, self.width,
                                        self.alpha_start - pi / 2, self.alpha_end - pi / 2, n_inner, n_outer)
        self.connection_edges = [n_inner + n_outer, n_inner]
        empty_arc1 = self._get_solid_arc(self.center, self.R - (self.width + self.gap) / 2,
                                         self.gap, self.alpha_start - pi / 2, self.alpha_end - pi / 2, n_inner, n_outer)

        empty_arc2 = self._get_solid_arc(self.center, self.R + (self.width + self.gap) / 2,
                                         self.gap, self.alpha_start - pi / 2, self.alpha_end - pi / 2, n_inner, n_outer)
        self.metal_region.insert(SimplePolygon().from_dpoly(metal_arc))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc1))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc2))

    def _refresh_named_connections(self):
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.center = self.connections[2]

    def _refresh_named_angles(self):
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]


class CPW2CPW(ElementBase):
    def __init__(self, Z0, Z1, start, end, trans_in=None):
        self.Z0 = Z0
        self.Z1 = Z1
        self.start = start
        self.end = end
        self.dr = self.end - self.start
        super().__init__(start, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def init_regions(self):
        self.connections = [DPoint(0, 0), DPoint(self.dr.abs(), 0)]
        self.angle_connections = [0, 0]
        alpha = atan2(self.dr.y, self.dr.x)
        self.angle_connections = [alpha, alpha]
        alpha_trans = DCplxTrans(1, alpha * 180 / pi, False, 0, 0)

        m_poly = DSimplePolygon([DPoint(0, -self.Z0.width / 2), DPoint(self.dr.abs(), -self.Z1.width / 2),
                                 DPoint(self.dr.abs(), self.Z1.width / 2), DPoint(0, self.Z0.width / 2)])
        e_poly1 = DSimplePolygon([DPoint(0, -self.Z0.b / 2), DPoint(self.dr.abs(), -self.Z1.b / 2),
                                  DPoint(self.dr.abs(), -self.Z1.width / 2), DPoint(0, -self.Z0.width / 2)])
        e_poly2 = DSimplePolygon([DPoint(0, self.Z0.b / 2), DPoint(self.dr.abs(), self.Z1.b / 2),
                                  DPoint(self.dr.abs(), self.Z1.width / 2), DPoint(0, self.Z0.width / 2)])

        m_poly.transform(alpha_trans)
        e_poly1.transform(alpha_trans)
        e_poly2.transform(alpha_trans)

        self.metal_region.insert(SimplePolygon.from_dpoly(m_poly))
        self.empty_region.insert(SimplePolygon.from_dpoly(e_poly1))
        self.empty_region.insert(SimplePolygon.from_dpoly(e_poly2))

    def _refresh_named_connections(self):
        self.start = self.connections[0]
        self.end = self.connections[1]

    def _refresh_named_angles(self):
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]


class Coil_type_1(ComplexBase):
    def __init__(self, Z0, start, L1, r, L2, trans_in=None):
        self.Z0 = Z0
        self.L1 = L1
        self.r = r
        self.L2 = L2
        super().__init__(start, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def init_primitives(self):
        self.cop1 = CPW(self.Z0.width, self.Z0.gap, DPoint(0, 0), DPoint(self.L1, 0))
        self.arc1 = CPW_arc(self.Z0, self.cop1.end, -self.r, -pi)
        self.cop2 = CPW(self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end - DPoint(self.L2, 0))
        self.arc2 = CPW_arc(self.Z0, self.cop2.end, -self.r, pi)

        self.connections = [self.cop1.start, self.arc2.end]
        self.angle_connections = [self.cop1.alpha_start, self.arc2.alpha_end]
        self.primitives = {"cop1": self.cop1, "arc1": self.arc1, "cop2": self.cop2, "arc2": self.arc2}


from collections import Counter


class CPW_RL_Path(ComplexBase):

    def __init__(self, origin, shape, cpw_parameters, turn_radiuses,
                 segment_lengths, turn_angles, trans_in=None, bridged=False):
        """
        A piecewise-linear coplanar waveguide with rounded turns.

        Segment lengths are treated as the lengths of the segments of
        a line with turn_raduises = 0. Changing turning raduises
        will not alter the position of the end of the line.

        TODO: 180 deg turns

        Parameters
        ----------
        origin : DPoint
            The point where the line should start
        shape : str
            String in format "RLLRL" where an R means a turn
            and an L means a straight part
        cpw_parameters : Union[CPWParameters, List[CPWParameters]]
            Parameters of the CPW or an array-like with parameters
            for each peace (R or L)
        turn_radiuses : Union[float, List[float]]
            Radius of the turns or an array-like with radiuses for
            each turn
        segment_lengths: list[float]
            Lengths of the straight parts of the equivalent
            piecewise-linear line with no corner rounding
        turn_angles: list[float]
            Angles for each turn of the line in radians
            !!! 180 turns are not yet supported !!!
        trans_in: DTrans
            Transformation of the line as a whole

        Returns
        -------

        """
        self._shape_string = shape
        self._N_elements = len(shape)
        self._shape_string_counter = Counter(shape)
        self._bridged = bridged

        self._N_turns = self._shape_string_counter['R']
        self._N_straights = self._shape_string_counter['L']
        if hasattr(cpw_parameters, "__len__"):
            if len(cpw_parameters) != self._N_elements:
                raise ValueError("CPW parameters dimension mismatch")
            else:
                self._cpw_parameters = copy.deepcopy(cpw_parameters)
        else:
            self._cpw_parameters = [cpw_parameters] * self._N_elements

        if hasattr(turn_radiuses, "__len__"):
            if len(turn_radiuses) != self._N_turns:
                raise ValueError("Turn raduises dimension mismatch")
            else:
                self._turn_radiuses = copy.deepcopy(turn_radiuses)
        else:
            self._turn_radiuses = [turn_radiuses] * self._N_turns
        if hasattr(segment_lengths, "__len__"):
            if len(segment_lengths) != self._N_straights:
                raise ValueError("Straight segments dimension mismatch")
            else:
                self._segment_lengths = copy.deepcopy(segment_lengths)
        else:
            self._segment_lengths = [segment_lengths] * self._N_straights

        if hasattr(turn_angles, "__len__"):
            if len(turn_angles) != self._N_turns:
                raise ValueError("Turn angles dimension mismatch")
            self._turn_angles = turn_angles
        else:
            self._turn_angles = [turn_angles] * self._N_turns
        super().__init__(origin, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def init_primitives(self):
        R_index = 0
        L_index = 0
        origin = DPoint(0, 0)

        prev_primitive_end = origin
        prev_primitive_end_angle = 0

        for i, symbol in enumerate(self._shape_string):
            if symbol == 'R':
                if (self._turn_angles[R_index] > 0):
                    turn_radius = self._turn_radiuses[R_index]
                else:
                    turn_radius = -self._turn_radiuses[R_index]

                if self._bridged:
                    cpw_arc = BridgedCPWArc(self._cpw_parameters[i], prev_primitive_end,
                                            turn_radius, self._turn_angles[R_index],
                                            200e3,
                                            trans_in=DCplxTrans(1, prev_primitive_end_angle * 180 / pi, False, 0, 0))
                else:
                    cpw_arc = CPW_arc(self._cpw_parameters[i], prev_primitive_end,
                                      turn_radius, self._turn_angles[R_index],
                                      trans_in=DCplxTrans(1, prev_primitive_end_angle * 180 / pi, False, 0, 0))

                self.primitives["arc_" + str(R_index)] = cpw_arc
                R_index += 1
            elif symbol == 'L':
                # Turns are reducing segments' lengths so as if there were no roundings at all
                # next 'R' segment if exists
                if (i + 1 < self._N_elements
                        and self._shape_string[i + 1] == 'R'
                        and abs(self._turn_angles[R_index]) < pi):
                    coeff = abs(tan(self._turn_angles[R_index] / 2))
                    self._segment_lengths[L_index] -= self._turn_radiuses[R_index] * coeff
                # previous 'R' segment if exists
                if (i - 1 > 0
                        and self._shape_string[i - 1] == 'R'
                        and abs(self._turn_angles[R_index - 1]) < pi):
                    coeff = abs(tan(self._turn_angles[R_index - 1] / 2))
                    self._segment_lengths[L_index] -= self._turn_radiuses[R_index - 1] * coeff

                if self._bridged:
                    cpw = BridgedCPW(self._cpw_parameters[i].width, self._cpw_parameters[i].gap, 200e3,
                                     prev_primitive_end, prev_primitive_end + DPoint(self._segment_lengths[L_index], 0),
                                     trans_in=DCplxTrans(1, prev_primitive_end_angle * 180 / pi, False, 0, 0))
                else:
                    # if( self._segment_lengths[L_index] < 0 ):
                    #    print(self._segment_lengths[L_index])
                    #    print("CPW_RL_Path warning: segment length is less than zero")
                    #    print("L_index = {}".format(L_index))
                    cpw = CPW(self._cpw_parameters[i].width, self._cpw_parameters[i].gap,
                              prev_primitive_end, prev_primitive_end + DPoint(self._segment_lengths[L_index], 0),
                              trans_in=DCplxTrans(1, prev_primitive_end_angle * 180 / pi, False, 0, 0))

                self.primitives["cpw_" + str(L_index)] = cpw
                L_index += 1

            prev_primitive = list(self.primitives.values())[-1]
            prev_primitive_end = prev_primitive.end
            prev_primitive_end_angle = prev_primitive.alpha_end

        self.connections = [list(self.primitives.values())[0].start,
                            list(self.primitives.values())[-1].end]
        self.angle_connections = [list(self.primitives.values())[0].alpha_start,
                                  list(self.primitives.values())[-1].alpha_end]

    def _refresh_named_connections(self):
        self.start = self.connections[0]
        self.end = self.connections[1]

    def _refresh_named_angles(self):
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def get_total_length(self):
        return sum(self._segment_lengths) + \
               sum([abs(R * alpha) for R, alpha in zip(self._turn_radiuses, self._turn_angles)])


class Bridge1(ElementBase):
    """
        Class implements bridges that are used to suppress
        non-TEM modes in coplanar or other types of waveguides.
        based on this design:
        https://drive.google.com/file/d/1nHM9lJNT9sBIWH9isRc_zKL6hUPwhhnP/view?usp=sharing
    """
    bridge_width = 20e3
    surround_gap = 8e3
    gnd_touch_dx = 20e3
    gnd_touch_dy = 10e3
    transition_len = 12e3
    gnd2gnd_dy = 70e3

    def __init__(self, center, gnd_touch_dx=20e3, trans_in=None):
        self.center = center
        self.gnd_touch_dx = gnd_touch_dx
        self.angle = 0
        super().__init__(center, trans_in)

        self._geometry_parameters = OrderedDict(
            [
                # TODO: add other members
                ("gnd_touch_dx, um", self.gnd_touch_dx / 1e3)
            ]
        )

    def init_regions(self):
        self.metal_regions["bridges_1"] = Region()  # region with ground contacts
        self.empty_regions["bridges_1"] = Region()  # remains empty

        self.metal_regions["bridges_2"] = Region()  # remains empty
        self.empty_regions["bridges_2"] = Region()  # region with erased bridge area

        center = DPoint(0, 0)
        self.connections = [center]
        self.angle_connections = [0]

        # init metal region of ground touching layer
        top_gnd_center = center + DPoint(0, self.gnd2gnd_dy / 2 + self.gnd_touch_dy / 2)
        p1 = top_gnd_center + DPoint(-self.gnd_touch_dx / 2, -self.gnd_touch_dy / 2)
        p2 = p1 + DVector(self.gnd_touch_dx, self.gnd_touch_dy)
        top_gnd_touch_box = pya.DBox(p1, p2)
        self.metal_regions["bridges_1"].insert(pya.Box().from_dbox(top_gnd_touch_box))

        bot_gnd_center = center + DPoint(0, -(self.gnd2gnd_dy / 2 + self.gnd_touch_dy / 2))
        p1 = bot_gnd_center + DPoint(-self.gnd_touch_dx / 2, -self.gnd_touch_dy / 2)
        p2 = p1 + DVector(self.gnd_touch_dx, self.gnd_touch_dy)
        bot_gnd_touch_box = pya.DBox(p1, p2)
        self.metal_regions["bridges_1"].insert(pya.Box().from_dbox(bot_gnd_touch_box))

        # init empty region for second layout layer
        # points start from left-bottom corner and goes in clockwise direction
        p1 = bot_gnd_touch_box.p1 + DPoint(-self.surround_gap, -self.surround_gap)
        p2 = p1 + DPoint(0, self.surround_gap + self.gnd_touch_dy +
                         self.transition_len - self.surround_gap)
        # top left corner + `surrounding_gap` + `transition_length`
        p3 = bot_gnd_touch_box.p1 + DPoint(0, bot_gnd_touch_box.height()) + \
             DPoint(-(20e3-self.gnd_touch_dx)/2, self.transition_len)
        bl_pts_list = [p1, p2, p3]  # bl stands for bottom-left
        ''' exploiting symmetry of reflection at x and y axes. '''
        # reflecting at x-axis
        tl_pts_list = list(map(lambda x: DTrans.M0 * x, bl_pts_list))  # tl stands for top-left
        # preserving order
        tl_pts_list = reversed(list(tl_pts_list))  # preserving clockwise points order
        # converting iterator to list
        l_pts_list = list(itertools.chain(bl_pts_list, tl_pts_list))  # l stands for left

        # reflecting all points at y-axis
        r_pts_list = list(map(lambda x: DTrans.M90 * x, l_pts_list))
        r_pts_list = list(reversed(r_pts_list))  # preserving clockwise points order

        # gathering points
        pts_list = l_pts_list + r_pts_list  # concatenating proper ordered lists

        empty_polygon = DSimplePolygon(pts_list)
        self.empty_regions["bridges_2"].insert(SimplePolygon.from_dpoly(empty_polygon))

    def _refresh_named_connections(self):
        self.center = self.connections[0]

    def _refresh_named_angles(self):
        self.angle = self.angle_connections[0]

    @staticmethod
    def bridgify_CPW(cpw, bridges_step, dest=None, bridge_layer1=-1,
                     bridge_layer2=-1, dest2=None,
                     avoid_points=[], avoid_distance=0):
        """
            Function puts bridge patterns to fabricate bridges on coplanar waveguide
        `cpw` with bridges having period of `bridges_step` along coplanar's wave
        propagation direction.
            Bridges are distributed over coplanar starting with its center.

        Parameters
        ----------
        cpw : Union[CPW, CPW_arc, CPW_RL_Path]
            instance of coplanar class to be bridged during fabrication
        bridges_step : float
            distance between centers of bridges in nm
        dest : pya.Cell
            cell to place bridge polygons at
        bridge_layer1 : int
            index of the layer in the `cell` with ground touching polygons
        bridge_layer2 : int
            index of the layer in the `cell` with empty polygons
        avoid_points : list[Union[DPoint,Point,Vector, DVector]]
            list points that you wish to keep bridges away
        avoid_distance : float
            distance in nm where there will be no bridges
            near the `avoid_points`
        Returns
        -------
        None
        """
        bridge_tmp = Bridge1(DPoint(0, 0))
        bridge_tmp.__bridgify_CPW(
            cpw, bridges_step,
            dest=dest, bridge_layer1=bridge_layer1,
            bridge_layer2=bridge_layer2, dest2=dest2,
            avoid_points=avoid_points, avoid_distance=avoid_distance
        )

    def __bridgify_CPW(self, cpw, bridges_step, dest=None,
                       bridge_layer1=-1, bridge_layer2=-1, dest2=None,
                       avoid_points=[], avoid_distance=0):
        """
            Function puts bridge patterns to fabricate bridges on coplanar waveguide
        `cpw` with bridges having period of `bridges_step` along coplanar's wave
        propagation direction.
            Bridges are distributed over coplanar starting with its center.

        Parameters
        ----------
        cpw : Union[CPW, CPW_arc, CPW_RL_Path]
            instance of coplanar class to be bridged during fabrication
        bridges_step : float
            distance between centers of bridges in nm
        dest : pya.Cell
            cell to place bridge polygons at
        bridge_layer1 : int
            index of the layer in the `cell` with ground touching polygons
        bridge_layer2 : int
            index of the layer in the `cell` with empty polygons
        avoid_points : list[Union[DPoint,Point,Vector, DVector]]
            list points that you wish to keep bridges away
        avoid_distance : float
            distance in nm where there will be no bridges
            near the `avoid_points`
        Returns
        -------
        None
        """
        if isinstance(cpw, CPW):
            # recursion base
            alpha = atan2(cpw.dr.y, cpw.dr.x)
            cpw_len = cpw.dr.abs()
            if cpw_len < (self.bridge_width + self.surround_gap):
                return

            cpw_dir_unit_vector = cpw.dr / cpw.dr.abs()

            # bridge with some initial dimensions
            tmp_bridge = Bridge1(DPoint(0, 0))
            bridge_width = tmp_bridge.gnd_touch_dx + 2 * tmp_bridge.surround_gap

            # number of additional bridges on either side of center
            additional_bridges_n = int((cpw_len / 2 - bridge_width / 2) // bridges_step)
            bridge_centers = []
            for i in range(-additional_bridges_n, additional_bridges_n + 1):
                bridge_center = cpw.start + (cpw_len / 2 + i * bridges_step) * cpw_dir_unit_vector

                avoid = False
                for avoid_point in avoid_points:
                    if (avoid_point - bridge_center).abs() < avoid_distance:
                        avoid = True
                        break

                if not avoid:
                    bridge_centers.append(bridge_center)

            bridges = []
            for center in bridge_centers:
                bridges.append(
                    Bridge1(
                        center,
                        trans_in=DCplxTrans(1, alpha / pi * 180, False, 0, 0)
                    )
                )
            for bridge in bridges:
                bridge.place(dest=dest, layer_i=bridge_layer1, region_name="bridges_1")
                if dest2 is not None:
                    bridge.place(dest=dest2, layer_i=bridge_layer2, region_name="bridges_2")
                else:
                    bridge.place(dest=dest, layer_i=bridge_layer2, region_name="bridges_2")
        elif isinstance(cpw, CPW_arc):
            # recursion base
            # to be implemented
            pass
        elif isinstance(cpw, CPW_RL_Path) or isinstance(cpw, Coil_type_1):
            for name, primitive in cpw.primitives.items():
                if isinstance(primitive, CPW):
                    Bridge1.bridgify_CPW(
                        primitive, bridges_step,
                        dest, bridge_layer1, bridge_layer2, dest2=dest2,
                        avoid_points=avoid_points, avoid_distance=avoid_distance
                    )
        else:
            # do nothing for other shapes
            return
