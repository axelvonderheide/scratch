
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from enthought.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from math  import \
     pow, fabs

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity, average, sqrt

from scipy.linalg import \
     inv, det

import time

from ibvpy.fets.fets_eval import FETSEval
from ibvpy.mats.mats_eval import MATSEval

from enthought.tvtk.api import tvtk

#------------------------------------------------------------------------------
# FETS2D9Q - 9 nodes isoparametric quadrilateral (2D, quadratic, Lagrange family) 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Element Information: 
#------------------------------------------------------------------------------
#
# Here an isoparametric element formulation is applied.
# The implemented shape functions are derived in femple 
# based on the following ordering of the nodes of the 
# parent element: 
#
#    _node_coord_map_dof = Array( Float, (9,2), 
#                                 [[ -1.,-1. ],
#                                  [  1.,-1. ],
#                                  [  1., 1. ],
#                                  [ -1., 1. ],
#                                  [  0.,-1. ],
#                                  [  1., 0. ],
#                                  [  0., 1. ],
#                                  [ -1., 0. ],
#                                  [  0., 0. ]])
#
#    _node_coord_map_geo = Array( Float, (9,2), 
#                                 [[ -1.,-1. ],
#                                  [  1.,-1. ],
#                                  [  1., 1. ],
#                                  [ -1., 1. ],
#                                  [  0.,-1. ],
#                                  [  1., 0. ],
#                                  [  0., 1. ],
#                                  [ -1., 0. ],
#                                  [  0., 0. ]])
#
#------------------------------------------------------------------------------

class FETS2D9Q(FETSEval):
    debug_on = True
    
    mats_eval = Instance(MATSEval)

    # Dimensional mapping
    dim_slice = slice(0,2)
    
    n_e_dofs = Int(9*2)
    t = Float( 1.0, label = 'thickness' )
    E = Float( 1.0, label = "Young's modulus" )
    nu = Float( 0., label = "Poisson's ratio" )

    # Integration parameters
    #
    ngp_r = 3
    ngp_s = 3

    field_entity_type = 'quad'
    # 
    vtk_r = [[-1.,-1.],
                        [ 0.,-1.],
                        [ 1.,-1.],
                        [-1., 0.],
                        [ 0., 0.],
                        [ 1., 0.],
                        [-1., 1.],
                        [ 0., 1.],
                        [ (1.-sqrt(0.5)), (1-sqrt(0.5))]]
    
    vtk_point_ip_map = [6,0,3,10,4,7,9,5,9]
    field_faces = [[0,1,4,3],
                   [1,2,5,4],
                   [3,4,7,6],
                   [4,5,8,7]]
    
    n_nodal_dofs = Int(2)
    
    gp_list = Property(Array,depends_on = 'ngp_r,ngp_s,ngp_t')
    @cached_property
    def _get_gp_list(self):
        points = tvtk.Points()
        points.insert_next_point(-1, -1, 0)
        points.insert_next_point(1, -1, 0)
        points.insert_next_point(1, 0, 0)
        points.insert_next_point((1-sqrt(0.5)),(1-sqrt(0.5)) , 0)
        points.insert_next_point(0, 1, 0)
        points.insert_next_point(-1, 1, 0)

        # Create a polydata with the points we just created.
        profile = tvtk.PolyData(points=points)


        # Perform a 2D Delaunay triangulation on them.
        delny = tvtk.Delaunay2D(input=profile)
        tri= delny.output
        tri.update()
        triangles = array(tri.polys.data,dtype=int_)    
        points = array(tri.points.data)
        #self.vtk_r = points[:,-1]
        ids = (triangles.reshape((triangles.shape[0]/4),4))[:,1:]
        #self.field_faces=ids
        int_order = 3
        gps=[]
        if int_order == 1:
            for id in ids:
                gp=[average(points[ix_(id)],0),1.]
                #print "gp ",gp
                gps.append(gp)
        elif int_order == 2:    
            raise NotImplementedError
        elif int_order == 3:    
            weigths = array([[0.6,0.2,0.2],[0.2,0.6,0.2],[0.2,0.2,0.6]])
            for id in ids:
                gps +=[average(points[ix_(id)],0),-0.5625],\
                    [average(points[ix_(id)],0, weigths[0]),0.52083333333333337],\
                    [average(points[ix_(id)],0, weigths[1]),0.52083333333333337],\
                    [average(points[ix_(id)],0, weigths[2]),0.52083333333333337]
    
                #gps.append(gp)
        elif int_order == 4:    
            raise NotImplementedError
        elif int_order == 5:    
            weigths = array([[0.0597158717, 0.4701420641, 0.4701420641],\
                             [0.4701420641, 0.0597158717, 0.4701420641],\
                             [0.4701420641, 0.4701420641, 0.0597158717],\
                             [0.7974269853, 0.1012865073, 0.1012865073],\
                             [0.1012865073, 0.7974269853, 0.1012865073],\
                             [0.1012865073, 0.1012865073, 0.7974269853]])
            for id in ids:
                weigts_sum = False#for debug
                gps +=[average(points[ix_(id)],0),0.225],\
                    [average(points[ix_(id)],0, weigths[0],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[1],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[2],weigts_sum),0.1323941527],\
                    [average(points[ix_(id)],0, weigths[3],weigts_sum),0.1259391805],\
                    [average(points[ix_(id)],0, weigths[4],weigts_sum),0.1259391805],\
                    [average(points[ix_(id)],0, weigths[5],weigts_sum),0.1259391805]
    
                #gps.append(gp)
        else:    
            raise NotImplementedError
        print "gps ",gps
        return gps

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r_pnt):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        N_geo_mtx = zeros((1,9), dtype = 'float_')
        N_geo_mtx[0,0] =   (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0
        N_geo_mtx[0,1] =   (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[1]) * ( 1 + r_pnt[0])) / 4.0
        N_geo_mtx[0,2] =   (r_pnt[0] * r_pnt[1] * ( 1 + r_pnt[1]) * ( 1 + r_pnt[0])) / 4.0
        N_geo_mtx[0,3] =   (r_pnt[0] * r_pnt[1] * ( 1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0
        N_geo_mtx[0,4] = - (r_pnt[1] * (-1 + r_pnt[0]) * (1 + r_pnt[0]) * (-1 + r_pnt[1])) / 2.0
        N_geo_mtx[0,5] = - (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[1]) * ( 1 + r_pnt[0])) / 2.0
        N_geo_mtx[0,6] = - (r_pnt[1] * (-1 + r_pnt[0]) * (1 + r_pnt[0]) * ( 1 + r_pnt[1])) / 2.0
        N_geo_mtx[0,7] = - (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0
        N_geo_mtx[0,8] =   (-1 + r_pnt[1]) * (1 + r_pnt[1]) * (-1 + r_pnt[0]) * (1 + r_pnt[0])
        return N_geo_mtx

    def get_dNr_geo_mtx(self, r_pnt):
        '''
        Return the matrix of shape function derivatives.
        Used for the construction of the Jacobi matrix.
        '''
        dNr_geo_mtx = zeros((2,9), dtype = 'float_')
        dNr_geo_mtx[0,0] =   (r_pnt[1] * (-1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[1])) / 4.0
        dNr_geo_mtx[0,1] =   (r_pnt[1] * (-1 + r_pnt[1]) * (1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[1])) / 4.0
        dNr_geo_mtx[0,2] =   (r_pnt[1] * (1 + r_pnt[1]) * (1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (1 + r_pnt[1])) / 4.0
        dNr_geo_mtx[0,3] =   (r_pnt[1] * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (1 + r_pnt[1])) / 4.0
        dNr_geo_mtx[0,4] = - (r_pnt[1] * (-1 + r_pnt[1]) * (1 + r_pnt[0])) / 0.2e1 -  (r_pnt[1] * (-1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[0,5] = - ((-1 + r_pnt[1]) * (1 + r_pnt[1]) * (1 + r_pnt[0])) / 2.0 -  (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[1])) / 2.0
        dNr_geo_mtx[0,6] = - (r_pnt[1] * (1 + r_pnt[1]) * (1 + r_pnt[0])) / 2.0 -  (r_pnt[1] * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[0,7] = - ((-1 + r_pnt[1]) * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0 -  (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[1])) / 2.0
        dNr_geo_mtx[0,8] =   (-1 + r_pnt[1]) * (1 + r_pnt[1]) * (1 + r_pnt[0]) + (-1 + r_pnt[1]) * (1 + r_pnt[1]) * (-1 + r_pnt[0])
        dNr_geo_mtx[1,0] =   (r_pnt[0] * (-1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[0])) / 4.0
        dNr_geo_mtx[1,1] =   (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (1 + r_pnt[0])) / 4.0
        dNr_geo_mtx[1,2] =   (r_pnt[0] * (1 + r_pnt[1]) * (1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (1 + r_pnt[0])) / 4.0
        dNr_geo_mtx[1,3] =   (r_pnt[0] * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 4.0 +  (r_pnt[0] * r_pnt[1] * (-1 + r_pnt[0])) / 4.0
        dNr_geo_mtx[1,4] = - ((-1 + r_pnt[0]) * (1 + r_pnt[0]) * (-1 + r_pnt[1])) / 2.0 -  (r_pnt[1] * (-1 + r_pnt[0]) * (1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[1,5] = - (r_pnt[0] * (1 + r_pnt[1]) * (1 + r_pnt[0])) / 2.0 -  (r_pnt[0] * (-1 + r_pnt[1]) * (1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[1,6] = - ((-1 + r_pnt[0]) * (1 + r_pnt[0]) * (1 + r_pnt[1])) / 2.0 -  (r_pnt[1] * (-1 + r_pnt[0]) * (1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[1,7] = - (r_pnt[0] * (1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0 -  (r_pnt[0] * (-1 + r_pnt[1]) * (-1 + r_pnt[0])) / 2.0
        dNr_geo_mtx[1,8] =   (-1 + r_pnt[0]) * (1 + r_pnt[0]) * (1 + r_pnt[1]) + (-1 + r_pnt[0]) * (1 + r_pnt[0]) * (-1 + r_pnt[1])
        return dNr_geo_mtx

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self,r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r_pnt.
        '''
        N = self.get_N_geo_mtx(r_pnt)    
        I_mtx = identity(self.n_nodal_dofs, float)
        N_mtx_list = [I_mtx*N[0,i] for i in range(0,N.shape[1])]
        N_mtx = hstack(N_mtx_list)
        return N_mtx

    def get_dNr_mtx(self, r_pnt):
        '''
        Return the derivatives of the shape functions used for the field approximation
        '''
        dNr_mtx = self.get_dNr_geo_mtx(r_pnt)    
        return dNr_mtx
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        Bx_mtx = zeros( (3,18 ), dtype = 'float_' )
        for i in range(0,9):
            Bx_mtx[0,i*2]   = dNx_mtx[0,i]
            Bx_mtx[1,i*2+1] = dNx_mtx[1,i]
            Bx_mtx[2,i*2]   = dNx_mtx[1,i]
            Bx_mtx[2,i*2+1] = dNx_mtx[0,i]
        return Bx_mtx


#----------------------- example --------------------

if __name__ == '__main__':
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, MGridDomain, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
        
    #from lib.mats.mats2D.mats_cmdm2D.mats_mdm2d import MACMDM
    #from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    #fets_eval = FETS2D9Q(mats_eval = MATS2DScalarDamage(strain_norm_type = 'Euclidean')) 
    #fets_eval = FETS2D9Q(mats_eval = MACMDM())  
    fets_eval = FETS2D9Q(mats_eval = MATS2DElastic()) 

    # Tseval for a discretized line domain
    tseval  = DOTSEval( fets_eval = fets_eval )

        # Define a mesh domain adaptor as a cached property to 
    # be constracted on demand
    from ibvpy.mesh.mgrid_domain import MeshGridAdaptor
    mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 2, 
                                     n_e_nodes_geo = (2,2,0), 
                                     n_e_nodes_dof = (2,2,0), 
                                     node_map_geo = [0,2,8,6,1,5,7,3,4], 
                                     node_map_dof = [0,2,8,6,1,5,7,3,4] )  

    # Discretization
    domain = MGridDomain( lengths = (3.,3.,0.), 
                             shape = (1,1,0),
                             adaptor = mgrid_adaptor )
                
    right_dof = 2
    ts = TS( tse = tseval,
         sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
         bcond_list =  [ BCDof(var='u', dof = i, value = 0.) for i in  domain.get_left_dofs()[:,0]  ] +
                    [ BCDof(var='u', dof = i, value = 0.) for i in [domain.get_left_dofs()[0,1]] ] +    
                    [ BCDof(var='u', dof = i, value = 0.002 ) for i in domain.get_right_dofs()[:,0] ],
         rtrace_list =  [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = right_dof,
                               var_x = 'U_k', idx_x = right_dof,
                               record_on = 'update'),
                         RTraceDomainField(name = 'Stress' ,
                         var = 'sig_app', idx = 0,
                         record_on = 'update'),
                         RTraceDomainField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
#                             RTraceDomainField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')
                      
                ]             
            )
    
    # Add the time-loop control
                #
    tl = TLoop( tstepper = ts,
                DT = 0.5,
                tline  = TLine( min = 0.0,  max = 1.0 ))
    
    tl.eval()
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tl )
    app.main()