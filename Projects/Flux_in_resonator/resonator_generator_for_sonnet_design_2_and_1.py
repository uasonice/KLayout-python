# $description: test_script  
# $version: 1.0.0  
# $show-in-menu
import pya
from math import cos, sin, atan2, pi
from pya import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans
from classLib import *
from classLib import revUtility


def fill_holes( poly, dx=10e3, dy=8e3, width=5e3, height=5e3, d=0 ):
	bbox = poly.bbox()
	poly_reg = Region( poly )
	t_reg = Region()

	y = bbox.p1.y + height
	while( y < bbox.p2.y - height ):
		x = bbox.p1.x + width
		while( x < bbox.p2.x -  width ):
			box = pya.Box().from_dbox( pya.DBox(DPoint(x,y), DPoint(x + width,y + height)) )
			x += dx
			t_reg.clear()
			t_reg.insert( box )
			qwe = t_reg.select_inside( poly_reg )
			if( qwe.is_empty() ):
				continue
			edge_pairs = qwe.inside_check( poly_reg, d )
			if( edge_pairs.is_empty() ):
				poly.insert_hole( box )
		y += dy
	return poly


class CHIP:
	dx = 0.6e6
	dy = 1.8e6
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
	x = 0.25*CHIP.dx
	y = 0.9*CHIP.dy
	p1 = DPoint( 0, y )
	p2 = DPoint( CHIP.dx, y )
	Z0 = CPW( width, gap, p1, p2 )
	Z0.place( cell, layer_photo )


	# resonator
	L_coupling = 300e3
	# corresponding to resonanse freq is linspaced in interval [6,9) GHz
	L1_list = [261653., 239550., 219456., 201109., 184291., 168819., 154536., 141312., 129032., 117600.]
	r = 50e3
	L2 = 270e3
	gnd_width = 35e3
	N = 4
	width_res = 9.6e3
	gap_res = 5.2e3
	Z_res = CPW( width_res, gap_res, origin, origin )
	to_line = 20e3
	worm = EMResonator_TL2Qbit_worm( Z_res, DPoint(x, y - Z0.width/2 - Z_res.width/2 - to_line) , L_coupling, L1_list[0], r, L2, N )
	worm.place( cell, layer_photo )

	# placing holes around the resonator
	'''
	reg = worm.gnd_reg
	new_reg = Region()
	#reg.merged_semantics=False
	r_cell = Region( cell.begin_shapes_rec( layer_photo ) )    
	reg = reg & r_cell
	for poly in reg:
		poly_temp = fill_holes( poly )
		new_reg.insert( poly_temp )
		
	#r_cell.merged_semantics=False   
	r_cell = r_cell - reg
	r_cell = r_cell + new_reg
	temp_i = cell.layout().layer( pya.LayerInfo(PROGRAM.LAYER1_NUM,0) ) 
	cell.shapes( temp_i ).insert( r_cell )
	cell.layout().clear_layer( layer_photo )
	cell.layout().move_layer( temp_i, layer_photo )
	cell.layout().delete_layer( temp_i )
	'''
	## DRAWING SECTION END ##
	if lv:
		lv.zoom_fit()
	else:
		layout.write("t.gds")
