from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVPSolve as IS, DOTSEval, DOTSListEval

from ibvpy.api import \
    RTraceDomainListField, RTraceDomainListInteg

from ibvpy.mats.mats1D5.mats1D5bond_elastic_frictional import MATS1D5Bond
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
    MATS3DMicroplaneDamage, PhiFnStrainSoftening

from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U

from ibvpy.mesh.fe_grid import FEGrid

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import array, tensordot, dot, zeros, c_
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos

fineness_width = 4 #  elemes in all dims
fineness_height = fineness_width * 1.5
fineness_bond = 5
fineness_no_bond = 2
bond_length = 0.003 # [m]
no_bond_length = 0.001 # [m]
width = 0.0035 # [m]
height = width * 1.5 # [m]

fets_reinf = FETS1D52B4ULRH(mats_eval = MATS1D5Bond(Ef = 17000.,
                                                    Af = 2.65e-6/4.,
                                                    Am = 0.,
                                                    Em = 0.,
                                                    tau_max = 8.23 * 2,
                                                    tau_fr = 8.23  * 2 ,
                                                    s_cr = 0.030e-3 * 10 ))        
fe_domain1 = FEGrid( coord_min = (0.,-0.0003, 0),
                     coord_max = (bond_length, 0 ,0.), 
                     shape   = ( fineness_bond *  1 ,1 ),
                     fets_eval = fets_reinf )

# characteristic element size

h =  bond_length / fineness_bond
print 'h',h

concrete = MATS3DMicroplaneDamage( model_version = 'stiffness',
                                   E = 34e3,
                                   nu = 0.2,
                                   phi_fn = PhiFnStrainSoftening( G_f = 0.001117,
                                                                  h = h,
                                                                  f_t = 2.8968 )
                                   )

fets_eval_mdm = FETS3D8H( mats_eval = concrete )

fets_eval2 = fets_eval_mdm
fets_eval2.vtk_r *= 0.8

# Geometric transformation
#
def gradual_mesh( points ):
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    
    L = x[-1] - x[0]
    B = y[-1] - y[0]
    H = z[-1] - z[0]

    X = x
    Y = y * (y/B)
    Z = z * (z/H)
    return c_[ X,Y,Z ]


fe_domain2 = FEGrid( coord_min = (0.,0.,0.),  
                     geo_transform = gradual_mesh,
                           coord_max = ( bond_length,height,width), 
                           shape   = ( fineness_bond, fineness_height, fineness_width ),
                           fets_eval = fets_eval2 )

fe_domain3 = FEGrid( coord_min = (bond_length,0.,0.),  
                     geo_transform = gradual_mesh,
                           coord_max = ( bond_length + no_bond_length,height,width), 
                           shape   = ( fineness_no_bond, fineness_height, fineness_width ),
                           fets_eval = fets_eval2 )

bc3_y = BCSlice( var='u', value = 0., dims = [2], 
                 slice = fe_domain3[:,1:,0,:,:,0] )
bc3_z = BCSlice( var='u', value = 0., dims = [1], 
                 slice = fe_domain3[:,0,1:,:,0,:] )

#time_fn = MFnLineArray( #xdata = arange(10),
#                       ydata = array([1.e-5]) )
loading_dof = fe_domain1.get_bottom_right_dofs()[0][0,0]

# original parameters from - CTRS 
# u_max = 12e-5 # 30.e-5 
# n_steps = 15
 
# parameters in order to enter the inelastic regime 
# - for the verification of the boundary conditions 
u_max = 16e-5
n_steps = 2

# parameters for spread damage -
# u_max = 32e-5
# n_steps = 10

loading_dof = fe_domain1.get_bottom_right_dofs()[0][0,0]

ts = TS( sdomain = [ fe_domain1, fe_domain2, fe_domain3 ], 
         dof_resultants = True,
         bcond_list =  [    
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain1.get_top_dofs,
                                    get_link_dof_method = fe_domain2.get_bottom_back_dofs,
                                    link_coeffs = [1.] ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain2.get_left_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [1],
                                    get_dof_method = fe_domain2.get_bottom_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [2],
                                    get_dof_method = fe_domain2.get_back_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
                                    get_link_dof_method = fe_domain3.get_left_dofs,
                                    get_dof_method = fe_domain2.get_right_dofs,
                                    link_coeffs = [1.] ),
                                    bc3_y,
                                    bc3_z,
#                        BCDofGroup( var='u', value = 0., dims = [1],
#                                    get_dof_method = fe_domain3.get_bottom_dofs ),
#                        BCDofGroup( var='u', value = 0., dims = [2],
#                                    get_dof_method = fe_domain3.get_back_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain1.get_bottom_left_dofs ),
                        BCDofGroup( var='u', value = u_max, dims = [0],
                                    #time_function = time_fn.get_value,
                                    #time_function = time_fn.get_value,
                                    get_dof_method = fe_domain1.get_bottom_right_dofs ),
                                     ],
         rtrace_list =  [ RTraceGraph( name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = loading_dof,
                               var_x = 'U_k', idx_x = loading_dof),                              
                        RTraceDomainListField(name = 'fracture_energy' ,
                                        var = 'fracture_energy', idx = 0, warp = True,
                                        record_on = 'update'), 
                        RTraceDomainListField(name = 'Debonding' ,
                                        var = 'debonding', idx = 0, warp = True,
                                        record_on = 'update'),            
                        RTraceDomainListInteg(name = 'Total fracture energy' ,
                                        var = 'fracture_energy', idx = 0, warp = False,
                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Strain' ,
#                                        var = 'eps_app', idx = 0, warp = True,
#                                        record_on = 'update'),   
#                        RTraceDomainListField(name = 'Stress' ,
#                                        var = 'sig_app', idx = 0, warp = True,
#                                        record_on = 'update'),                 
                ]             
            )

# Add the time-loop control
tloop = TLoop( tstepper = ts, debug = False,
               tolerance = 1e-5,
               tline  = TLine( min = 0.0,  step = 1.0/float(n_steps) , max = 1.0 ))

print tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    
