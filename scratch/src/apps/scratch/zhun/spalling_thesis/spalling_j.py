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

from ibvpy.mesh.fe_grid_domain import FEGridDomain
from ibvpy.mesh.fe_domain_list import FEDomainList
from ibvpy.mesh.fe_domain_tree import FEDomainTree

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
            
        T = array( [[1, 0, 0],
                    [0, ca, sa],
                    [0, -sa, ca]], dtype = 'float_' )
        sig_rotated = dot( T, dot( sig_tensor, T.T ) ) 
        return sig_rotated.flatten()
    
#t_fr = 13.e-6
#s_cr = t_fr/5.2e-3
#bond_fn = MFnLineArray(xdata = [0.,s_cr,1.01*s_cr, 10.*s_cr],
#                       ydata = [0.,t_fr,t_fr/200.,t_fr/300.])


t_fr = 12.5
s_cr = 0.028#
bond_fn = MFnLineArray(xdata = [0.,s_cr],
                       ydata = [0.,t_fr])

fineness = 2 #  elemes in all dims
fineness_x_adhesive = 5
fineness_x_cohesive = 2
fets_eval1 = FETS1D52B4ULRH(mats_eval = MATS1D5Bond(Ef = 17000.,
                                                    Af = 1.76e-6/4.,
                                                    Am = 0.,
                                                    Em = 0.,
                                                    bond_fn = bond_fn))        
fe_domain1 = FEGridDomain( coord_min = (0.,-0.0003, 0),
                           coord_max = (0.0035, 0 ,0.), 
                           shape   = ( fineness_x_adhesive *  2 ,1 ),
                           fets_eval = fets_eval1 )


#fets_eval2 = FETS3D8H(mats_eval = MATS3DElasticSigRotated(E = 35000, nu = 0.2),
#                      vtk_r = [[ -0.8, -0.8, -0.8],
#                        [  0.8, -0.8, -0.8],
#                        [ -0.8,  0.8, -0.8],
#                        [  0.8,  0.8, -0.8],
#                        [ -0.8, -0.8,  0.8],
#                        [  0.8, -0.8,  0.8],
#                        [ -0.8,  0.8,  0.8],
#                        [  0.8,  0.8,  0.8]]
# )    
fets_eval2 = FETS3D8H20U(mats_eval = MATS3DElasticSigRotated(E = 35000, nu = 0.2),
#                             # Used for Visualization 
#          vtk_r = [[ -.8, -.8, -.8],
#                        [  .8, -.8, -.8],
#                        [  .8,  .8, -.8],
#                        [ -.8,  .8, -.8],
#                        [ -.8, -.8,  .8],
#                        [  .8, -.8,  .8],
#                        [  .8,  .8,  .8],
#                        [ -.8,  .8,  .8],
#                        [  0., -.8, -.8],
#                        [  .8,  0., -.8],
#                        [  0.,  .8, -.8],
#                        [ -.8,  0., -.8],
#                        [  0., -.8,  .8],
#                        [  .8,  0.,  .8],
#                        [  0.,  .8,  .8],
#                        [ -.8,  0.,  .8],
#                        [ -.8, -.8,  0.],
#                        [  .8, -.8,  0.],
#                        [  .8,  .8,  0.],
#                        [ -.8,  .8,  0.]]
 )

#fets_eval2 = FETS3D8H(mats_eval = MATS2DMicroplaneDamage( dimensionality = '3D' ))
fe_domain2 = FEGridDomain( coord_min = (0.,0.,0.),  
                           coord_max = ( 0.0035,0.0035,0.0035), 
                           shape   = ( fineness_x_adhesive, fineness, fineness ),
                           fets_eval = fets_eval2 )

fe_domain3 = FEGridDomain( coord_min = (0.0035,0.,0.),  
                           coord_max = ( 0.0049,0.0035,0.0035), 
                           shape   = ( fineness_x_cohesive, fineness, fineness ),
                           fets_eval = fets_eval2 )

fe_domain  = FEDomainList( subdomains = [ fe_domain1, fe_domain2, fe_domain3 ] )
fe_domain_tree = FEDomainTree( domain_list = fe_domain )


#time_fn = MFnLineArray( #xdata = arange(10),
#                       ydata = array([1.e-5]) )


ts = TS( sdomain = fe_domain_tree, 
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
                                    
                        BCDofGroup( var='u', value = 1.e-5, dims = [0],
                                    #time_function = time_fn.get_value,
                                    #time_function = time_fn.get_value,
                                    get_dof_method = fe_domain1.get_bottom_right_dofs ),
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
                        RTraceDomainListField(name = 'Stress_radial' ,
                                        var = 'sig_app', idx = 4, warp = True,
#                                        position = 'int_pnts',
                                        record_on = 'update'),
                        RTraceDomainListField(name = 'Stress_ring' ,
                                        var = 'sig_app', idx = 8, warp = True,
                                        #position = 'int_pnts',
                                        record_on = 'update'),
                        RTraceDomainListField(name = 'Stress_ring_ip' ,
                                        var = 'sig_app', idx = 8, warp = True,
                                        position = 'int_pnts',
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



# Add the time-loop control
tloop = TLoop( tstepper = ts,
               #tolerance = 1e-3,
               tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))


tloop.eval()

from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp( ibv_resource = tloop )
#app.main()    
