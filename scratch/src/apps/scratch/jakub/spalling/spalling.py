from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, BCDofGroup, IBVPSolve as IS, DOTSEval, DOTSListEval

from ibvpy.rtrace.rt_domain_list_field import RTraceDomainListField

#from ibvpy.mats.mats1D5.mats1D5bond import MATS1D5Bond
from ibvpy.mats.mats1D5.mats1D5bond_elastic_frictional import MATS1D5Bond
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import MATS3DMicroplaneDamage

from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U

from ibvpy.mesh.fe_grid import FEGrid

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import array, tensordot, dot,zeros
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos

class MATS3DElasticSigRotated( MATS3DElastic ):

    def get_sig_app(self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        sig_tensor = map3d_sig_eng_to_mtx( sig_eng )           
        # get the transformation matrix
        
        sctx.r_pnt = sctx.loc
        X_pnt = sctx.fets_eval.get_X_pnt( sctx )
        X, Y, Z = X_pnt
        
        L = sqrt( Y**2 + Z**2 )
        if L <= 1e-8:
            sig_another_one = 0.
            sa = 0.
            ca = 0.
        else:
            sa = Z / L
            ca = Y / L
            
            s_alpha = asin(Z / L)
            c_alpha = acos(Z / L)
    
            n = array([ 0, -sa, ca ],dtype='float_' )
            sig_another_one = tensordot( n, dot( n, sig_tensor ), [0,0] )

        T = array( [[1, 0, 0],
                    [0, ca, sa],
                    [0, -sa, ca]], dtype = 'float_' )
        sig_rotated = dot( T, dot( sig_tensor, T.T ) ) 
        return sig_rotated

fineness_width = 4 #  elemes in all dims
fineness_height = fineness_width * 1.5
fineness_bond = 5
fineness_no_bond = 2
bond_length = 0.003
no_bond_length = 0.001 
width = 0.0035
height = width * 1.5

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

fets_eval_elastic = FETS3D8H(mats_eval = MATS3DElasticSigRotated(E = 35000, nu = 0.2) )    

# characteristic element size
G_f = 0.001117 * 0.5

h =  bond_length / fineness_bond

concrete = MATS3DMicroplaneDamage( 
                                 model_version = 'stiffness',
                                 E = 34e3,
                                 nu = 0.2,
                                 f_t = 2.8968,
                                  )
#concrete.polar_fn.G_f = G_f / h
# Default settings for 'polar_fn_class':
concrete.polar_fn_class = 'Isotropic polar function'
# Default settings of 'PolarFnBase' for 'phi_fn_class':
concrete.polar_fn.phi_fn_class = 'QuasiBrittle'

fets_eval_mdm = FETS3D8H( mats_eval = concrete )

fets_eval2 = fets_eval_mdm
fets_eval2.vtk_r *= 0.8

fe_domain2 = FEGrid( coord_min = (0.,0.,0.),  
                           coord_max = ( bond_length,height,width), 
                           shape   = ( fineness_bond, fineness_height, fineness_width ),
                           fets_eval = fets_eval2 )

fe_domain3 = FEGrid( coord_min = (bond_length,0.,0.),  
                           coord_max = ( bond_length + no_bond_length,height,width), 
                           shape   = ( fineness_no_bond, fineness_height, fineness_width ),
                           fets_eval = fets_eval2 )

#time_fn = MFnLineArray( #xdata = arange(10),
#                       ydata = array([1.e-5]) )


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
                        BCDofGroup( var='u', value = 0., dims = [1],
                                    get_dof_method = fe_domain3.get_bottom_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [2],
                                    get_dof_method = fe_domain3.get_back_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain1.get_bottom_left_dofs ),
                        BCDofGroup( var='u', value = 15.e-5, dims = [0],
                                    #time_function = time_fn.get_value,
                                    #time_function = time_fn.get_value,
                                    get_dof_method = fe_domain1.get_bottom_right_dofs ),
                                     ],
         rtrace_list =  [ 
#                         RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                               var_y = 'F_int', idx_y = 0,
#                               var_x = 'U_k', idx_x = 1),
#                          RTraceDomainListField(name = 'Fracture Energy' ,
#                          var = 'fracture_energy', idx = 0,
#                          record_on = 'update',
#                          warp = True),
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
                        RTraceDomainListField(name = 'fracture_energy' ,
                                        var = 'fracture_energy', idx = 0, warp = True,
                                        record_on = 'update'),            
                        RTraceDomainListField(name = 'strain' ,
                                        var = 'eps_app', idx = 0, warp = False,
#                                        position = 'int_pnts',
                                        record_on = 'update'),
                        RTraceDomainListField(name = 'Stress_rotated' ,
                                        var = 'sig_app', idx = 0, warp = True,
#                                        position = 'int_pnts',
                                        record_on = 'update'),
                        
                ]             
            )



# Add the time-loop control
tloop = TLoop( tstepper = ts,
               tolerance = 1e-5,
               tline  = TLine( min = 0.0,  step = 1.0 / 15, max = 1.0 ))


print tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
app.main()    
