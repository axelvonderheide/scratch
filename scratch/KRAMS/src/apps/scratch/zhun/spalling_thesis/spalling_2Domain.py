# 
from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainField, TLoop, \
    TLine, BCDof, BCDofGroup, IBVPSolve as IS, DOTSEval, DOTSListEval

from ibvpy.rtrace.rt_domain_list_field import RTraceDomainListField

from ibvpy.mats.mats1D5.mats1D5bond import MATS1D5Bond
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import MATS3DScalarDamage
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage

from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import FETS3D8H27U

from ibvpy.mesh.fe_grid import FEGrid

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from numpy import array, tensordot, dot,zeros
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos

class MATS3DElasticSigRotated( MATS3DElastic ): # calculate the ring stress

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
        
        T = array( [[1, 0, 0],
                    [0, ca, sa],
                    [0, -sa, ca]], dtype = 'float_' )
        sig_rotated = dot( T, dot( sig_tensor, T.T ) ) 
        return sig_rotated
    

# the stiffness of the bond function
t_fr = 1.e-7
s_cr = 1.e-4     #stiffness=t_fr/s_cr
bond_fn = MFnLineArray(xdata = [0., s_cr, 1.5e-4,   2.e-4 ],
                       ydata = [0., t_fr, 1.e-20,   1.e-20])        # here must a number not zero.

fineness = 3 #  elemes in all dims
fineness_x_adhesive = 6 
#fineness_x_cohesive = 4    # used in 3.Domain  
fets_eval1 = FETS1D52B4ULRH(mats_eval = MATS1D5Bond(Ef = 17000.,
                                                    Af = 1.76e-6/4.,
                                                    Am = 0.,
                                                    Em = 0.,
                                                    bond_fn = bond_fn))        
fe_domain1 = FEGrid( coord_min = (0.,-0.0003, 0),
                     coord_max = (0.004, 0 ,0.), 
                     shape   = ( fineness_x_adhesive *  2 ,1 ),
                     fets_eval = fets_eval1 )


#fets_eval2 = FETS3D8H(mats_eval = MATS3DElasticSigRotated(E = 35000, nu = 0.2),
#                      vtk_r = [[ -.99, -.99, -.99],
#                        [  .99, -.99, -.99],
#                        [ -.99,  .99, -.99],
#                        [  .99,  .99, -.99],
#                        [ -.99, -.99,  .99],
#                        [  .99, -.99,  .99],
#                        [ -.99,  .99,  .99],
#                        [  .99,  .99,  .99]]
# )    
fets_eval2 = FETS3D8H20U(mats_eval = MATS3DElasticSigRotated(E = 35000, nu = 0.2),
                             # Used for Visualization 
          vtk_r = [[ -.99, -.99, -.99],
                        [  .99, -.99, -.99],
                        [  .99,  .99, -.99],
                        [ -.99,  .99, -.99],
                        [ -.99, -.99,  .99],
                        [  .99, -.99,  .99],
                        [  .99,  .99,  .99],
                        [ -.99,  .99,  .99],
                        [  0., -.99, -.99],
                        [  .99,  0., -.99],
                        [  0.,  .99, -.99],
                        [ -.99,  0., -.99],
                        [  0., -.99,  .99],
                        [  .99,  0.,  .99],
                        [  0.,  .99,  .99],
                        [ -.99,  0.,  .99],
                        [ -.99, -.99,  0.],
                        [  .99, -.99,  0.],
                        [  .99,  .99,  0.],
                        [ -.99,  .99,  0.]]
 )

fe_domain2 = FEGrid( coord_min = (0.,0.,0.),  
                     coord_max = ( 0.004,0.004,0.004), 
                     shape   = ( fineness_x_adhesive, fineness, fineness ),
                     fets_eval = fets_eval2 )

#fe_domain3 = FEGridDomain( coord_min = (0.002,0.,0.),  
#                           coord_max = ( 0.004,0.004,0.004), 
#                           shape   = ( fineness_x_cohesive, fineness, fineness ),
#                           fets_eval = fets_eval2 )

time_fn = MFnLineArray( ydata = array([0., 15., 15.08, 16., 16.08, 17.,17.08, 18.,18.08, 19.,19.08, 20.,20.08,21.,21.08,22.,22.08]) )

#time_fn = MFnLineArray( ydata = array([0., 5., 8., 15., 15.05]) )

# boundary conditions
ts = TS( sdomain = [ fe_domain1, fe_domain2 ], 
         dof_resultants = True,
         bcond_list =  [    
                        BCDofGroup( var='u', value = 0., dims = [0],        # in X-direction no displacement
                                    get_dof_method = fe_domain1.get_top_dofs, # top dofs in domain 1
                                    get_link_dof_method = fe_domain2.get_bottom_back_dofs, #bottom back dofs in domain 2
                                    link_coeffs = [1.] ),                   # combine the dofs in the two domains
                                    
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain2.get_left_dofs ), #the left dofs in domain 2 have no displacements in X direction.
                                    
                        BCDofGroup( var='u', value = 0., dims = [1],
                                    get_dof_method = fe_domain2.get_bottom_dofs ),#the bottom dofs in domain 2 have no displacements in Y direction
                                    
                        BCDofGroup( var='u', value = 0., dims = [2],
                                    get_dof_method = fe_domain2.get_back_dofs ),#the back dofs in domain 2 have no displacements in Z direction 
                                    
#                        BCDofGroup( var='u', value = 0., dims = [0,1,2],
#                                    get_link_dof_method = fe_domain3.get_left_dofs,
#                                    get_dof_method = fe_domain2.get_right_dofs,
#                                    link_coeffs = [1.] ),          
#                        BCDofGroup( var='u', value = 0., dims = [1],
#                                    get_dof_method = fe_domain3.get_bottom_dofs ),
#                        BCDofGroup( var='u', value = 0., dims = [2],
#                                    get_dof_method = fe_domain3.get_back_dofs ),
                        BCDofGroup( var='u', value = 0., dims = [0],
                                    get_dof_method = fe_domain1.get_bottom_left_dofs ), # the bottom left dofs in domain 1 have no displacements in 
                                    
                        BCDofGroup( var='u', value = 1.e-5, dims = [0],
                                    #time_function = time_fn.get_value,
                                    time_function = time_fn.get_value,
                                    get_dof_method = fe_domain1.get_bottom_right_dofs ), # the bottom right dof in domain 1 have displacements of 1.e5m 
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
#                        RTraceDomainListField(name = 'Stress_radial' ,
#                                        var = 'sig_app', idx = 4, warp = True,
##                                        position = 'int_pnts',
#                                        record_on = 'update'),
                        RTraceDomainListField(name = 'Stress_ring' ,
                                        var = 'sig_app', idx = 8, warp = True,
#                                        position = 'int_pnts',
                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'Stress_ring_ip' ,
#                                        var = 'sig_app', idx = 8, warp = True,
#                                        position = 'int_pnts',      #show the Gauss points
#                                        record_on = 'update'),
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



# Add the time-loop control
tloop = TLoop( tstepper = ts,
               #tolerance = 1e-3,
               tline  = TLine( min = 0.0,  step = 1.0, max = 20.0 ))


tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
#app.main()     #starts Mayavi2
