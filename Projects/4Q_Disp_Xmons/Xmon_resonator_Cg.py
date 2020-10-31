# Enter your Python code here
import pya
from math import cos, sin, atan2, pi
from pya import Point, DPoint, DVector, DSimplePolygon, SimplePolygon, DPolygon, Polygon, Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from classLib.baseClasses import ComplexBase
from classLib.coplanars import CPW, CPW_arc, CPWParameters
from classLib.resonators import Coil_type_1
from classLib.shapes import XmonCross
from classLib.resonators import EMResonatorTL3QbitWormRLTailXmonFork

from sonnetSim.sonnetLab import SonnetLab, SonnetPort, SimulationBox


class CHIP:
    dx = 0.6e6
    dy = 1.8e6
    L1 = 2.5e6
    gap = 150.e3
    width = 260.e3
    b = 2 * gap + width
    origin = DPoint(0, 0)
    center = DPoint(dx/2, dy/2)
    box = pya.DBox(origin, origin + DPoint(dx, dy))
    # only 4 connections programmed by now
    connections = [box.p1 + DPoint(L1 + b / 2, 0), box.p1 + DPoint(dx - (L1 + b / 2), 0),
                   box.p2 - DPoint(L1 + b / 2, 0), box.p1 + DPoint(L1 + b / 2, dy)]


if __name__ == "__main__":
    # getting main references of the application
    app = pya.Application.instance()
    mw = app.main_window()
    lv = mw.current_view()
    cv = None

    # this insures that lv and cv are valid objects
    if lv is None:
        cv = mw.create_layout(1)
        lv = mw.current_view()
    else:
        cv = lv.active_cellview()

    # find or create the desired by programmer cell and layer
    layout = cv.layout()
    layout.dbu = 0.001
    if (layout.has_cell("testScript")):
        pass
    else:
        cell = layout.create_cell("testScript")

    layer_info_photo = pya.LayerInfo(10, 0)
    layer_info_el = pya.LayerInfo(1, 0)
    layer_photo = layout.layer(layer_info_photo)
    layer_el = layout.layer(layer_info_el)

    # setting layout view  
    lv.select_cell(cell.cell_index(), 0)
    lv.add_missing_layers()

    ## DRAWING SECTION START ##
    chip_box = pya.Box(Point(0, 0), Point(CHIP.dx, CHIP.dy))
    origin = DPoint(0, 0)

    contact_L = 1e6

    # main drive line coplanar
    width = 24.1e3
    gap = 12.95e3
    x = 0.25 * CHIP.dx
    y = 0.9 * CHIP.dy
    p1 = DPoint(0, y)
    p2 = DPoint(CHIP.dx, y)
    Z0 = CPW(width, gap, p1, p2)

    L_coupling_list = [1e3 * x for x in [255, 250, 250, 240, 230]]
    # corresponding to resonanse freq is linspaced in interval [6,9) GHz
    L0 = 1600e3
    L1_list = [1e3 * x for x in [37.7039, 67.6553, 90.925, 81.5881, 39.9021]]
    r = 60e3
    N = 5
    L2_list = [r] * len(L1_list)
    L4_list = [r] * len(L1_list)
    width_res = 20e3
    gap_res = 10e3
    Z_res = CPWParameters(width_res, gap_res)
    to_line_list = [53e3] * len(L1_list)
    fork_metal_width = 25e3
    fork_gnd_gap = 20e3
    xmon_dys_Cg_coupling = [1e3*x for x in [10, 12, 14, 16, 18]]
    # resonator-fork parameters
    # -20e3 for Xmons in upper sweet-spot
    # -10e3 for Xmons in lower sweet-spot
    xmon_fork_penetration_list = [-25e3]

    # xmon cross parameters
    cross_len_x = 180e3
    cross_width_x = 60e3
    cross_gnd_gap_x = 20e3
    cross_len_y = 60e3
    cross_width_y = 60e3
    cross_gnd_gap_y = 20e3
    xmon_x_distance = 485e3

    ### RESONATORS TAILS CALCULATIONS SECTION START ###
    # key to the calculations can be found in hand-written format here:
    # https://drive.google.com/file/d/1wFmv5YmHAMTqYyeGfiqz79a9kL1MtZHu/view?usp=sharing

    # distance between nearest resonators central conductors centers
    resonators_d = 430e3
    # x span between left long vertical line and
    # right-most center of central conductors
    resonators_widths = [2 * r + L_coupling for L_coupling in L_coupling_list]
    x1 = sum(resonators_widths[:2]) + 2 * resonators_d + resonators_widths[3] / 2 - 2 * xmon_x_distance
    x2 = x1 + xmon_x_distance - (resonators_widths[0] + resonators_d)
    x3 = sum(resonators_widths[:3]) + 3 * resonators_d - (x1 + 3 * xmon_x_distance)
    x4 = sum(resonators_widths[:4]) + 4 * resonators_d - (x1 + 4 * xmon_x_distance)

    res_tail_shape = "LRLRL"
    tail_turn_radiuses = r

    L2_list[0] += 6 * Z_res.b
    L2_list[1] += 0
    L2_list[3] += 3 * Z_res.b
    L2_list[4] += 6 * Z_res.b

    L4_list[1] += 6 * Z_res.b
    L4_list[2] += 6 * Z_res.b
    L4_list[3] += 3 * Z_res.b
    tail_segment_lengths_list = [
        [L2_list[0], x1, L4_list[0]],
        [L2_list[1], x2, L4_list[1]],
        [L2_list[2], (L_coupling_list[2] + 2 * r) / 2, L4_list[2]],
        [L2_list[3], x3, L4_list[3]],
        [L2_list[4], x4, L4_list[4]]
    ]
    tail_turn_angles_list = [
        [pi / 2, -pi / 2],
        [pi / 2, -pi / 2],
        [pi / 2, -pi / 2],
        [-pi / 2, pi / 2],
        [-pi / 2, pi / 2],
    ]
    tail_trans_in_list = [
        Trans.R270,
        Trans.R270,
        Trans.R270,
        Trans.R270,
        Trans.R270
    ]
    ### RESONATORS TAILS CALCULATIONS SECTION END ###

    pars = list(
        zip(
            L1_list, to_line_list, L_coupling_list,
            xmon_fork_penetration_list,
            tail_segment_lengths_list, tail_turn_angles_list, tail_trans_in_list,
            xmon_fork_gnd_gaps
        )
    )
    resonator_index = 0
    for res_idx, params in [list(enumerate(pars))[resonator_index]]:
        # parameters exctraction
        L1 = params[0]
        to_line = params[1]
        L_coupling = params[2]
        xmon_fork_penetration = params[3]
        tail_segment_lengths = params[4]
        tail_turn_angles = params[5]
        tail_trans_in = params[6]
        xmon_fork_gnd_gap = params[7]
        # clear this cell and layer
        cell.clear()

        fork_x_span = cross_width_x + 2 * (xmon_fork_gnd_gap + fork_metal_width)
        # `fork_y_span` based on coupling modulated with
        # xmon_fork_penetration from `xmon_fork_penetration`
        fork_y_span = xmon_fork_penetration + xmon_fork_gnd_gap

        worm = EMResonatorTL3QbitWormRLTailXmonFork(
                Z_res, origin, L_coupling, L0, L1, r, N,
                tail_shape=res_tail_shape, tail_turn_radiuses=tail_turn_radiuses,
                tail_segment_lengths=tail_segment_lengths,
                tail_turn_angles=tail_turn_angles, tail_trans_in=tail_trans_in,
                fork_x_span=fork_x_span, fork_y_span=fork_y_span,
                fork_metal_width=fork_metal_width, fork_gnd_gap=fork_gnd_gap
        )
        xmon_center = (worm.fork_y_cpw1.end + worm.fork_y_cpw2.end)/2
        xmon_center += DPoint(0, -(cross_len_y + cross_width_y / 2) + xmon_fork_penetration)
        xmonCross = XmonCross(
            xmon_center,
            sideX_length=cross_len_x, sideX_width=cross_width_x, sideX_gnd_gap=cross_gnd_gap_x,
            sideY_length=cross_len_y, sideY_width=cross_width_y, sideY_gnd_gap=cross_gnd_gap_y
        )
        # cutting everything out except for the tail that will be capcitavely coupled
        # to Xmon
        for key in list(worm.primitives.keys()):
            if key.startswith("fork") or key.startswith("erased"):
                pass
            elif (key == "cpw_end_open_RLPath"):
                for key_lvl2 in list(worm.primitives[key].primitives.keys())[:-2]:
                    del worm.primitives[key].primitives[key_lvl2]
            else:
                del worm.primitives[key]

        # calculate bounding box dimensions
        tmp_reg = Region()
        worm.place(tmp_reg)
        xmonCross.place(tmp_reg)
        bbox = tmp_reg.bbox()

        # calculate simulation box dimensions
        CHIP.dx = 3*bbox.width()
        CHIP.dy = 2*bbox.height()
        chip_p1 = DPoint(bbox.center()) + DPoint(-0.5*CHIP.dx, -0.5*CHIP.dy)
        chip_p2 = DPoint(bbox.center()) + DPoint(0.5*CHIP.dx, 0.5*CHIP.dy)
        chip_box = pya.Box(DPoint(0, 0), DPoint(CHIP.dx, CHIP.dy))

        # placing all objects in proper order and translating them to origin
        cell.shapes(layer_photo).insert(chip_box)

        # translating all objects so that chip.p1 at coordinate origin
        translation_trans = DCplxTrans(-DVector(chip_p1))

        xmonCross.make_trans(translation_trans)
        xmonCross.place(cell, layer_photo)

        worm.make_trans(translation_trans)
        worm.place(cell, layer_photo)

        # forming resonator's tail upper side open-circuited for simulation
        p1 = list(worm.cpw_end_open_RLPath.primitives.values())[0].start
        print(list(worm.cpw_end_open_RLPath.primitives.values()))
        p2 = p1 + DPoint(-20e3, 0)
        empty_cpw_tail_start = CPW(0, Z_res.b / 2, p1, p2)
        # empty_cpw_tail_start.make_trans(translation_trans)
        print(empty_cpw_tail_start.start)
        empty_cpw_tail_start.place(cell, layer_photo)

        xmonCross_corrected = XmonCross(
            xmon_center,
            sideX_length=cross_len_x, sideX_width=cross_width_x, sideX_gnd_gap=cross_gnd_gap_x,
            sideY_length=cross_len_y, sideY_width=cross_width_y, sideY_gnd_gap=min(xmon_fork_gnd_gap, cross_gnd_gap_y)
        )
        xmonCross_corrected.make_trans(translation_trans)
        xmonCross_corrected.place(cell, layer_photo)

        ## DRAWING SECTION END ##
        lv.zoom_fit()

        ### MATLAB COMMANDER SECTION START ###
        ml_terminal = SonnetLab()
        print("starting connection...")
        from sonnetSim.cMD import CMD

        ml_terminal._send(CMD.SAY_HELLO)
        ml_terminal.clear()
        simBox = SimulationBox(CHIP.dx, CHIP.dy, 900, 1200)
        ml_terminal.set_boxProps(simBox)
        print("sending cell and layer")
        from sonnetSim.pORT_TYPES import PORT_TYPES

        ports = [
            SonnetPort(p1, PORT_TYPES.AUTOGROUNDED),
            SonnetPort(xmonCross.cpw_b.end, PORT_TYPES.AUTOGROUNDED)
        ]
        ml_terminal.set_ports(ports)

        ml_terminal.send_polygons(cell, layer_photo)
        ml_terminal.set_linspace_sweep(0.01, 0.01, 1)
        print("simulating...")
        result_path = ml_terminal.start_simulation(wait=True)
        ml_terminal.release()

        # get the .csv result file and exctract capcity of island in fF
        import os
        import csv
        project_dir = os.path.dirname(__file__)

        C_qr = None
        C1 = None

        with open(result_path.decode("ascii"), "r") as csv_file:
            data_rows = list(csv.reader(csv_file))
            ports_imps_row = data_rows[6]
            R = float(ports_imps_row[0].split(' ')[1])
            data_row = data_rows[8]
            freq0 = float(data_row[0])

            s = [[0, 0], [0, 0]]  # s-matrix
            for i in range(0, 2):
                for j in range(0, 2):
                    s[i][j] = complex(float(data_row[1 + 2 * (i * 2 + j)]), float(data_row[1 + 2 * (i * 2 + j) + 1]))
            import math

            # formula to transform S-params to Y-params is taken from
            # https://en.wikipedia.org/wiki/Admittance_parameters#Two_port

            # calculations for C1 and C12 are taken from Sonnet equation plotting
            delta = (1 + s[0][0]) * (1 + s[1][1]) - s[0][1] * s[1][0]
            y21 = -2 * s[1][0] / delta * 1 / R
            C_qr = 1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y21).imag)

            y11 = ((1 - s[0][0])*(1 + s[1][1]) + s[0][1]*s[1][0]) / delta * 1/R
            C1 = -1e15 / (2 * math.pi * freq0 * 1e9 * (1 / y11).imag)

        print(C_qr, "  ", C1)

        # ### SAVING RESULT SECTION START ###
        # output_filepath = os.path.join(project_dir, "Xmon_resonator_Cg.csv")
        # if os.path.exists(output_filepath):
        #     # append data to file
        #     with open(output_filepath, "a") as csv_file:
        #         writer = csv.writer(csv_file)
        #         writer.writerow(
        #             [cross_width / 1e3, cross_len / 1e3, cross_gnd_gap / 1e3,
        #              xmon_fork_penetration / 1e3, xmon_fork_gnd_gap / 1e3, C1, C_qr / C1]
        #         )
        # else:
        #     # create file, add header, append data
        #     with open(output_filepath, "w") as csv_file:
        #         writer = csv.writer(csv_file)
        #         # create header of the file
        #         writer.writerow(
        #             ["cross_width, um", "cross_len, um", "cross_gnd_gap, um",
        #              "xmon_fork_penetration, um", "xmon_fork_gnd_gap, um", "C1, fF", "beta"])
        #         writer.writerow(
        #             [cross_width / 1e3, cross_len / 1e3, cross_gnd_gap / 1e3,
        #              xmon_fork_penetration / 1e3, xmon_fork_gnd_gap/1e3, C1, C_qr / C1]
        #         )
        # ### SAVING RESULT SECTION END ###