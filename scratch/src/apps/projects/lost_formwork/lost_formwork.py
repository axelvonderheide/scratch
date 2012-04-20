'''
IBAC - lost formwork idealization for TRC

@todos 
1) geo_transform for the bottom bar [Done]

2) Include possibility to use slices for the get_dof_methods.
   (there is problem with combined application of constraints
    at the top of the slab)
3) avoid search for link_dofs
4) include the bottom reinforcement
5) provide the 2.5D model for CMDM - orient it in the slab to simulate distr cracking
6) calibrate the damage function from the tensile
7) apply the 3D mdm for the matrix
8) determine the bond law for the reinforcing bar?
9) prepare the response tracers that should be available 
   by default. this includes the setup of pipelines 
   to process the fracture energy and stresses
10) define a response tracer for load-deflection diagram 
   (requires accumulation of contributions 
    from all constraied points)
9) define an application plugin with the summarized 
   application parameters and an action menu with 
   the buttons for starting the calculation
   
12) Problem with the cyclic constraint - can be solved using open slices - cumbersome
'''

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup

from ibvpy.rtrace.rt_domain_list_field import RTraceDomainListField

from ibvpy.mats.mats1D5.mats1D5bond import MATS1D5Bond
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import MATS3DMicroplaneDamage

from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U

from ibvpy.mesh.fe_grid import FEGrid

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import array, tensordot, dot,zeros, c_
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos

fineness = 3 # elemes in all dims
fineness_x = 4
fineness_z = 2

mats_eval = MATS3DScalarDamage( E = 35000, nu = 0.2 )
mats_eval = MATS3DMicroplaneDamage( E = 35000, nu = 0.2, 
                                  model_version = 'stiffness' )
#mats_eval = MATS3DElastic(E = 35000, nu = 0.2)
fets_eval2 = FETS3D8H( mats_eval = mats_eval,
                       vtk_r = [[ -0.8, -0.8, -0.8],
                                     [  0.8, -0.8, -0.8],
                                     [ -0.8,  0.8, -0.8],
                                     [  0.8,  0.8, -0.8],
                                     [ -0.8, -0.8,  0.8],
                                     [  0.8, -0.8,  0.8],
                                     [ -0.8,  0.8,  0.8],
                                     [  0.8,  0.8,  0.8]])

span = 0.45
height_rib = 0.24
height_plate = 0.1
width_rib_bottom = 0.3
width_rib_top = 0.4
width_side_plate = 0.55
overspan = 0.3

def bar_z_widening( points ):
    x, y, z = points[:,0], points[:,1], points[:,2]
    h = height_rib
    w = width_rib_bottom
    dw = ( width_rib_top - width_rib_bottom ) / 2
    
    return c_[ x, y, ( z * ( w * h + 2 * dw * y + 2 * dw * h )  ) / ( w * h ) ]

fed_rib_left = FEGrid( coord_min = ( 0.0,  -height_rib, - width_rib_bottom / 2., 0.0 ),  
                           coord_max = ( span,         0.0,   width_rib_bottom / 2.  ),
                           geo_transform = bar_z_widening,
                           shape   = ( fineness_x, 1, 1 ),
                           fets_eval = fets_eval2 )

fed_slab_midleft = FEGrid( coord_min = ( 0.0,  0.0,  - width_rib_top / 2  ),  
                           coord_max = ( span, height_plate,  width_rib_top / 2. ), 
                           shape   = ( fineness_x, 1, 1 ),
                           fets_eval = fets_eval2 )

fed_slab_frontleft = FEGrid( coord_min = (  0.0, 0.0,  width_rib_top / 2.  ),  
                           coord_max = ( span, height_plate,  
                                         width_rib_top / 2. + width_side_plate ), 
                           shape   = ( fineness_x, 1, fineness_z ),
                           fets_eval = fets_eval2 )

fed_slab_backleft = FEGrid( coord_min = (  0.0,          0.0, - width_side_plate - width_rib_top / 2. ),  
                           coord_max = ( span, height_plate, - width_rib_top / 2. ), 
                           shape   = ( fineness_x, 1, fineness_z ),
                           fets_eval = fets_eval2 )

fed_rib_right = FEGrid( coord_min = ( span,  -height_rib, - width_rib_bottom / 2., 0.0 ),  
                           coord_max = ( span + overspan,         0.0,   width_rib_bottom / 2.  ),
                           geo_transform = bar_z_widening,
                           shape   = ( fineness_x, 1, 1 ),
                           fets_eval = fets_eval2 )

fed_slab_midright = FEGrid( coord_min = ( span,  0.0,  - width_rib_top / 2  ),  
                           coord_max = ( span + overspan, height_plate,  width_rib_top / 2. ), 
                           shape   = ( fineness_x, 1, 1 ),
                           fets_eval = fets_eval2 )

#time_fn = MFnLineArray( #xdata = arange(10),
#                       ydata = array([1.e-5]) )


ts = TS( sdomain = [fed_rib_left, 
                    # fed_rib_right, 
                    fed_slab_midleft, 
                    # fed_slab_midright,
                    fed_slab_frontleft, 
                    fed_slab_backleft 
                    ], 
         dof_resultants = True,
         bcond_list =  [    
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fed_rib_left.get_left_dofs ),

                        # right support - both in y and z direction
                        BCDofGroup( var='u', value = 0., dims = [1,2],
                                    get_dof_method = fed_rib_left.get_bottom_right_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fed_slab_midleft.get_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fed_slab_frontleft.get_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fed_slab_backleft.get_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
                                    get_link_dof_method = fed_slab_midleft.get_bottom_dofs,
                                    get_dof_method = fed_rib_left.get_top_dofs,
                                    link_coeffs = [1.] ),          
                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
                                    get_link_dof_method = fed_slab_frontleft.get_back_dofs,
                                    get_dof_method = fed_slab_midleft.get_front_dofs,
                                    link_coeffs = [1.] ),          
                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
                                    get_link_dof_method = fed_slab_backleft.get_front_dofs,
                                    get_dof_method = fed_slab_midleft.get_back_dofs,
                                    link_coeffs = [1.] ), 
                        # link the left 
#                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
#                                    get_link_dof_method = fed_rib_left.get_right_dofs,
#                                    get_dof_method = fed_rib_right.get_left_dofs,
#                                    link_coeffs = [1.] ),

#                        # glue the left and right ribs
#                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
#                                    get_link_dof_method = fed_slab_midright.get_bottom_dofs,
#                                    get_dof_method = fed_rib_right.get_top_dofs,
#                                    link_coeffs = [1.] ),
                                                            
                        # glue the left and right ribs
#                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
#                                    get_link_dof_method = fed_slab_midleft.get_right_dofs,
#                                    get_dof_method = fed_slab_midright.get_left_dofs,
#                                    link_coeffs = [1.] ),     
#                        BCDofGroup( var='u', value = -0.01, dims = [1],
#                                    get_dof_method = fed_slab_midleft.get_top_left_dofs ),
                        BCDofGroup( var='u', value = -0.0005, dims = [1],
                                    get_dof_method = fed_slab_frontleft.get_top_left_dofs ),
                        BCDofGroup( var='u', value = -0.0005, dims = [1],
                                    get_dof_method = fed_slab_backleft.get_top_left_dofs ),
#                        BCDofGroup( var='u', value = 0., dims = [2],
#                                    get_dof_method = fed_slab_backleft.get_back_dofs ),
                                     ],
         rtrace_list =  [ 
#                         RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = 0,
#                               var_x = 'U_k', idx_x = 1),
#                           RTraceDomainListField( name = 'G_f', 
#                                                     var = 'fracture_energy', 
#                                                     idx = 0, warp = False ),
#                        RTraceDomainListField(name = 'Damage' ,
#                                        var = 'microplane_damage', idx = 0, warp = True,
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Displacement' ,
#                                        var = 'u', idx = 0, warp = True,
#                                        record_on = 'update'),

#                        RTraceDomainListField(name = 'Main Stress' ,
#                                        var = 'msig_pm', idx = 0, warp = True,
#                                        position = 'int_pnts',
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Displacement' ,
#                                        var = 'u', idx = 0, warp = True,
#                                        record_on = 'update'),            
#                        RTraceDomainListField(name = 'Stress_x' ,
#                                        var = 'sig_app', idx = 0, warp = True,
#                                        #position = 'int_pnts',
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Stress_y' ,
#                                        var = 'sig_app', idx = 4, warp = True,
#                                        #position = 'int_pnts',
#                                        record_on = 'update'),
                        RTraceDomainListField(name = 'Stress_x' ,
                                        var = 'sig_app', idx = 0, warp = True,
#                                        position = 'int_pnts',
                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Sin' ,
#                                        var = 'sig_app', idx = 1, warp = True,
#                                        #position = 'int_pnts',
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Cos' ,
#                                        var = 'sig_app', idx = 2, warp = True,
#                                        #position = 'int_pnts',
#                                        record_on = 'update'),
                        
                ]             
            )


global tloop
# Add the time-loop control
tloop = TLoop( tstepper = ts,
               tolerance = 1e-3,
               tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))

use_profiling = True
start_ui = True
calculate = True

if calculate:
    if use_profiling:
        import cProfile
        cProfile.run('tloop.eval()', 'lost_formwork_tprof' )
        
        import pstats
        p = pstats.Stats('lost_formwork_tprof')
        p.strip_dirs()
        print 'cumulative'
        p.sort_stats('cumulative').print_stats(50)
        print 'time'
        p.sort_stats('time').print_stats(50)
    
    else:
        tloop.eval()
    
if start_ui:  
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
        
