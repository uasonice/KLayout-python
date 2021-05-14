# $description: test_script  
# $version: 1.0.0  
# $show-in-menu
import pya
from math import cos, sin, atan2, pi
from pya import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

import classLib
from classLib import *
from classLib import revUtility

from classLib.coplanars import *

class CHIP:
	dx = 10.1e6
	dy = 5.1e6
	L1 = 2.5e6
	gap = 150.e3
	width = 260.e3
	b = 2*gap + width
	origin = DPoint( 0,0 )
	box = pya.DBox( origin, origin + DPoint( dx,dy ) )
	# only 4 connections programmed by now
	connections = [box.p1 + DPoint( L1 + b/2,0 ), box.p1 + DPoint( dx - (L1+b/2),0 ), box.p2 - DPoint( L1 + b/2,0 ),  box.p1 + DPoint( L1 + b/2, dy )]


if ( __name__ ==  "__main__" ):
	# find or create the desired by programmer cell and layer
	layout, lv = revUtility.klayout_init_layout()
	layout.dbu = 0.001
	if( layout.has_cell( "testScript") ):
		pass
	else:
		cell = layout.create_cell( "testScript" )

	layer_info_photo = pya.LayerInfo(10,0)
	layer_info_el = pya.LayerInfo(1,0)
	layer_photo = layout.layer( layer_info_photo )
	layer_el = layout.layer( layer_info_el )

	# clear this cell and layer
	cell.clear()

	# setting layout view
	if lv:
		lv.select_cell(cell.cell_index(), 0)
		lv.add_missing_layers()


	## DRAWING SECTION START ##
	cell.shapes( layer_photo ).insert( pya.Box( Point(0,0), Point( CHIP.dx, CHIP.dy ) ) )
	origin = DPoint( 0,0 )

	contact_L = 1e6

	# main drive line coplanar
	width = 24.1e3
	gap = 12.95e3
	p1 = DPoint( 0 + contact_L, CHIP.dy/2 )
	p2 = DPoint( CHIP.dx - contact_L, CHIP.dy/2 )
	Z0 = CPW( width, gap, p1, p2 )
	Z0.place( cell, layer_photo )

	# left contact pad
	width1 = CHIP.width*2
	gap1 = CHIP.gap*2
	p3 = DPoint( 0, CHIP.dy/2 )
	p4 = DPoint( contact_L/2, CHIP.dy/2 )
	Z1 = CPW( width1, gap1, p3, p4 )
	Z1.place( cell, layer_photo )

	adapter1 = CPW2CPW( Z1, Z0, p4, p1 )
	adapter1.place( cell, layer_photo )

	# right contact pad
	width1 = CHIP.width*2
	gap1 = CHIP.gap*2
	p5 = DPoint( CHIP.dx - contact_L/2, CHIP.dy/2 )
	p6 = DPoint( CHIP.dx, CHIP.dy/2 )
	Z2 = CPW( width1, gap1, p5, p6 )
	Z2.place( cell, layer_photo )

	adapter1 = CPW2CPW( Z2, Z0, p5, p2 )
	adapter1.place( cell, layer_photo )

	# resonator
	L_coupling = list(map( lambda x: x*1e3, [387.398, 367.161, 335.258, 298.067, 273.721, 257.209, 244.608, 234.3, \
											 225.5] ))
	L1 = list(map( lambda x: x*1e3, [290.033, 262.729, 238.458, 216.743, 197.199, 179.516, 163.441, \
									 148.763, 135.309] ))
	r = 50e3
	L2 = 270e3
	toLine = 21e3
	N = 4
	width_res = 9.6e3
	gap_res = 5.2e3
	Z_res = CPW( width_res, gap_res, origin, origin )

	# qBit
	origin = DPoint(1e6,0)
	a = list(map(lambda x: 4.66*x*1e3,[1., 1.03078, 1.0933, 1.1914, 1.33203, 1.52603, 1.78943, 2.14545, \
									   2.62762]))
	b = a
	alpha = 0.43
	jos1_b = 487*alpha
	jos1_a = 0.1e3
	f1 = 0.13e3
	d1 = 0.05e3
	jos2_b = jos1_b/alpha
	jos2_a = 0.09e3
	f2 = 0.25e3
	d2 = d1
	w = 0.2e3
	dCap = 0
	dSquares = 30e3
	gap = 25e3
	square_a = 150e3
	alum_over = 20e3
	B1_width = 3e3
	qbit_gnd_gap = 15e3
	qbit_coupling_gap = 15e3
	qbit_params = [a,b,
				   jos1_b,jos1_a,f1,d1,
				   jos2_b,jos2_a,f2,d2,
				   w, dCap,gap, square_a, dSquares, alum_over, B1_width]


	# place bottom array of resonators with qBits
	N_bottom = 5
	delta = 1e6
	step_bot = (CHIP.dx - 2*contact_L - 2*delta)/(N_bottom-1)
	resonators = []
	qbits = []
	empties = []
	for i in range( N_bottom ):
		point = DPoint( contact_L + delta + i*step_bot, CHIP.dy/2 - ((Z0.width + Z_res.width)/2 + toLine)  )

		worm = EMResonator_TL2Qbit_worm( Z_res, point, L_coupling[i], L1[i], r, L2, N )
		worm.place( cell, layer_photo )
		resonators.append( worm )

		dy_qbit = (worm.cop_tail.dr.abs() - square_a)/2
		qbit_params = [a[i],b[i],
					   jos1_b,jos1_a,f1,d1,
					   jos2_b,jos2_a,f2,d2,
					   w, dCap,gap, square_a, dSquares, alum_over, B1_width]
		dy_qbit = (worm.cop_tail.dr.abs() - square_a)/2
		qbit_start = worm.end + DPoint( B1_width/2,0 )
		qbit = QBit_Flux_Сshunted_3JJ( worm.end + DPoint( qbit_coupling_gap + Z_res.width/2, dy_qbit ) , qbit_params, DCplxTrans( 1,90,False,0,0 ) )
		qbit_bbox = pya.DBox().from_ibox( qbit.metal_regions["photo"].bbox() )
		p1 = qbit.origin + DPoint(-qbit_coupling_gap,qbit.square_a/2)
		p2 = p1 + DPoint( qbit_bbox.width() + qbit_coupling_gap + qbit_gnd_gap, 0 )
		empty = CPW( 0, (square_a + 2*qbit_gnd_gap)/2, p1, p2 )
		empty.place( cell, layer_photo )
		empties.append( empty )

		qbit.place( cell, layer_photo, layer_el )
		qbits.append(qbit)

	# place top array of resonators with qBits
	N_top = N_bottom - 1
	delta = 1e6
	step_top = step_bot
	resonators = []
	qbits = []
	empties = []
	for i in range( N_top ):
		point = DPoint( contact_L + delta + step_bot/2 + i*step_top, CHIP.dy/2 + (Z0.width + Z_res.width)/2 + toLine )

		worm = EMResonator_TL2Qbit_worm( Z_res, point, L_coupling[i + N_bottom], L1[i + N_bottom], r, L2, N, Trans( Trans.M0 ) )
		worm.place( cell, layer_photo )
		resonators.append( worm )

		qbit_params = [a[i+N_bottom],b[i+N_bottom],
					   jos1_b,jos1_a,f1,d1,
					   jos2_b,jos2_a,f2,d2,
					   w, dCap,gap, square_a, dSquares, alum_over, B1_width]
		dy_qbit = -(worm.cop_tail.dr.abs() - square_a)/2 - square_a
		qbit_start = worm.end + DPoint( B1_width/2,0 )
		qbit = QBit_Flux_Сshunted_3JJ( worm.end + DPoint( qbit_coupling_gap + Z_res.width/2, dy_qbit ) , qbit_params, DCplxTrans( 1,90,False,0,0 ) )
		qbit_bbox = pya.DBox().from_ibox( qbit.metal_regions["photo"].bbox() )
		p1 = qbit.origin + DPoint(-qbit_coupling_gap,qbit.square_a/2)
		p2 = p1 + DPoint( qbit_bbox.width() + qbit_coupling_gap + qbit_gnd_gap, 0 )
		empty = CPW( 0, (square_a + 2*qbit_gnd_gap)/2, p1, p2 )
		empty.place( cell, layer_photo )
		empties.append( empty )

		qbit.place( cell, layer_photo, layer_el )
		qbits.append(qbit)


	## DRAWING SECTION END ##
	if lv:
		lv.zoom_fit()
	else:
		layout.write("t.gds")