import pya
from math import sqrt, cos, sin, atan2, pi, copysign
from pya import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from pya import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans


class PROGRAM:
    LAYER1_NUM = 7
    LAYER2_NUM = 8
    
class CHIP:
    dx = 10.1e6
    dy = 5.1e6
    L1 = 2466721
    gap = 150e3
    width = 260e3
    b = 2*gap + width
    origin = DPoint( 0,0 )
    box = pya.DBox( origin, origin + DPoint( dx,dy ) )
    # only 4 connections programmed by now
    connections = [box.p1 + DPoint( L1 + b/2,0 ), box.p1 + DPoint( dx - (L1+b/2),0 ), box.p2 - DPoint( L1 + b/2,0 ),  box.p1 + DPoint( L1 + b/2, dy )]
    

class Element_Base():
    def __init__(self, origin, trans_in=None):
        ## MUST BE IMPLEMENTED ##
        self.connections = []       # DPoint list with possible connection points
        self.angle_connections = [] #list with angle of connecting elements
        self.connection_ptrs = [] # pointers to connected structures represented by their class instances
        ## MUST BE IMLPEMENTED END ##
        self.origin = origin
        self.metal_region = Region()
        self.empty_region = Region()
        self.metal_region.merged_semantics = False
        self.empty_region.merged_semantics = False
        self.DCplxTrans_init = None
        self.ICplxTrans_init = None
        
        if( trans_in is not None ):
        # if( isinstance( trans_in, ICplxTrans ) ): <==== FORBIDDEN
            if( isinstance( trans_in, DCplxTrans ) ):
                self.DCplxTrans_init = trans_in
                self.ICplxTrans_init = ICplxTrans().from_dtrans( trans_in )
            elif( isinstance( trans_in, CplxTrans ) ):
                self.DCplxTrans_init = DCplxTrans().from_itrans( trans_in )
                self.ICplxTrans_init = ICplxTrans().from_trans( trans_in )
            elif( isinstance( trans_in, DTrans ) ):
                self.DCplxTrans_init = DCplxTrans( trans_in, 1 )
                self.ICplxTrans_init = ICplxTrans( Trans().from_dtrans( trans_in ), 1 )
            elif( isinstance( trans_in, Trans ) ):
                self.DCplxTrans_init = DCplxTrans( DTrans().from_itrans( trans_in ), 1 )
                self.ICplxTrans_init = ICplxTrans( trans_in, 1 )
                
        self._init_regions_trans()
    
    def init_regions( self ):
        raise NotImplementedError
    
    # first it makes trans_init displacement
    # then the rest of the trans_init
    # then displacement of the current state to the origin
    # after all, origin should be updated
    def _init_regions_trans( self ):
        self.init_regions()         # must be implemented in every subclass
        dr_origin = DSimplePolygon( [DPoint(0,0)] )
        if( self.DCplxTrans_init is not None ):
            # constructor trans displacement
            dCplxTrans_temp = DCplxTrans( 1,0,False, self.DCplxTrans_init.disp )
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )
            
            # rest of the constructor trans functions
            dCplxTrans_temp = self.DCplxTrans_init.dup()
            dCplxTrans_temp.disp = DPoint(0,0)
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )                  
            
        # translation to the old origin (self.connections are alredy contain proper values)
        self.make_trans( DCplxTrans( 1,0,False, self.origin ) ) # move to the origin
        self.origin += dr_origin.point( 0 )
        
        
    def make_trans( self, dCplxTrans ):
        if( dCplxTrans is not None ):
            iCplxTrans = ICplxTrans().from_dtrans( dCplxTrans )
            self.metal_region.transform( iCplxTrans )
            self.empty_region.transform( iCplxTrans )
            self._update_connections( dCplxTrans )
            self._update_alpha( dCplxTrans )
    
    def _update_connections( self, dCplxTrans ):       
        if( dCplxTrans is not None ):
            # the problem is, if k construct polygon with multiple points
            # their order in poly_temp.each_point() doesn't coinside with the 
            # order of the list that was passed to the polygon constructor
            # so, when k perform transformation and try to read new values through poly_temp.each_point()
            # they values are rearranged
            # solution is: k need to create polygon for each point personally, and the initial order presists
            for i,pt in enumerate(self.connections):
                poly_temp = DSimplePolygon( [pt] )
                poly_temp.transform( dCplxTrans )
                self.connections[i] = poly_temp.point( 0 )
    
    def _update_alpha( self, dCplxTrans ):
        if( dCplxTrans is not None ):
            dCplxTrans_temp = dCplxTrans.dup()
            dCplxTrans_temp.disp = DPoint(0,0)
            
            for i,alpha in enumerate(self.angle_connections):
                poly_temp = DSimplePolygon( [DPoint( cos(alpha), sin(alpha) )] )
                poly_temp.transform( dCplxTrans_temp )
                pt = poly_temp.point( 0 )
                self.angle_connections[i] = atan2( pt.y, pt.x )
    
    def _update_origin( self, dCplxTrans ):
        if( dCplxTrans is not None ):     
            poly_temp = DSimplePolygon( [self.origin] )
            poly_temp.transform( dCplxTrans )
            self.origin = poly_temp.point( 0 )
    
    def place( self, dest, layer_i=-1 ):
        r_cell = None
        if( layer_i != -1 ):            
            r_cell = Region( dest.begin_shapes_rec( layer_i ) )        

        # how to interpret destination
        if( layer_i == -1 ):
            dest += self.metal_region
            dest -= self.empty_region
        else:
            temp_i = dest.layout().layer( pya.LayerInfo(PROGRAM.LAYER1_NUM,0) ) 
            dest.shapes( temp_i ).insert( r_cell + self.metal_region  - self.empty_region )
            dest.layout().clear_layer( layer_i )
            dest.layout().move_layer( temp_i, layer_i )
            dest.layout().delete_layer( temp_i )
class CPW( Element_Base ):
    """@brief: class represents single coplanar waveguide
      @params:   float width - represents width of the central conductor
                        float gap - spacing between central conductor and ground planes
                        float gndWidth - width of ground plane to be drawed 
                        DPoint start - center aligned point, determines the start point of the coplanar segment
                        DPoint end - center aligned point, determines the end point of the coplanar segment
    """
    def __init__(self, width, gap, start, end, gndWidth=-1, trans_in=None ):                    
        self.width = width
        self.gap = gap
        self.b = 2*gap + width
        self.gndWidth = gndWidth
        self.end = end
        self.start = start
        self.dr = end - start
        super( CPW, self ).__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_regions( self ):
        self.connections = [DPoint(0,0),self.dr]
        self.start = DPoint(0,0)
        self.end = self.start + self.dr
        alpha = atan2( self.dr.y, self.dr.x )
        self.angle_connections = [alpha,alpha]
        alpha_trans = ICplxTrans().from_dtrans( DCplxTrans( 1,alpha*180/pi,False, self.start ) )
        self.metal_region.insert( pya.Box( Point().from_dpoint(DPoint(0,-self.width/2)), Point().from_dpoint(DPoint( self.dr.abs(), self.width/2 )) ) )
        self.empty_region.insert( pya.Box( Point().from_dpoint(DPoint(0,self.width/2)), Point().from_dpoint(DPoint( self.dr.abs(), self.width/2 + self.gap )) ) )
        self.empty_region.insert( pya.Box( Point().from_dpoint(DPoint(0,-self.width/2-self.gap)), Point().from_dpoint(DPoint( self.dr.abs(), -self.width/2 )) ) )
        self.metal_region.transform( alpha_trans )
        self.empty_region.transform( alpha_trans )



# Enter your Python code here
if( __name__ == "__main__"):
# getting main references of the application
    app = pya.Application.instance()
    mw = app.main_window()
    lv = mw.current_view()
    cv = None
    
    #this insures that lv and cv are valid objects
    if( lv == None ):
        cv = mw.create_layout(1)
        lv = mw.current_view()
    else:
        cv = lv.active_cellview()

# find or create the desired by programmer cell and layer
    layout = cv.layout()
    layout.dbu = 0.001
    if( layout.has_cell( "testScript") ):
        pass
    else:
        cell = layout.create_cell( "testScript" )
    
    info = pya.LayerInfo(1,0)
    layer_i = layout.layer( info )
    
    lv.select_cell(cell.cell_index(), 0)
    lv.add_missing_layers()
    
    
    
    region = Region()
    Z0 = CPW( 1e3,1e3, DPoint(0,0), DPoint(1e6,0) )
    Z0.place( region ) 
    print( region.bbox() )  # output: ()
    cell.shapes( layer_i ).insert( region ) # output: nothing is drawn
    
    
    '''region = Region()
    Z0 = CPW( 1e3,1e3, DPoint(0,0), DPoint(1e6,0) )
    region = region + Z0.metal_region - Z0.empty_region 
    print( region.bbox() ) # output: (0,-500,1e6,500)
    cell.shapes( layer_i ).insert( region ) #output: line is drawn as intended
    '''
    lv.zoom_fit()