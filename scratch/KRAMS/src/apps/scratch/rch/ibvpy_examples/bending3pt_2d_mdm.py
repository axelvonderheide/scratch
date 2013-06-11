
# Example - uniaxial tension with weakened cross section using Masar's isotropic damage model 

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, RTraceDomainListInteg, TLoop, \
    TLine, BCDof, IBVPSolve as IS, BCSlice
from ibvpy.api import \
    BCDofGroup
from ibvpy.mats.mats_proxy import \
    MATSProxy
from numpy import array, zeros
    
from ibvpy.mesh.fe_grid import FEGrid

def calculate_2D( length = 1., height = 0.3, shape = (21, 10) ):

    from ibvpy.fets.fets2D.fets2D4q import \
        FETS2D4Q
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
        MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import \
        MATS2DScalarDamage
    from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import \
        Energy, Euclidean, Mises, Rankine, Mazars
    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
        MATS2DMicroplaneDamage, PhiFnStrainSoftening
    
    h = length / float( shape[0] )

    mdm = MATS2DMicroplaneDamage( model_version = 'complience',
                                  E = 34000.,
                                  nu = 0.25,
                                  phi_fn = PhiFnStrainSoftening( G_f = 0.001117,
                                                                 f_t = 2.8968,
                                                                 h = h ) )
    
    fets_eval = FETS2D4Q( mats_eval = mdm ) 
    #fets_eval = FETS2D4Q( mats_eval = MATS2DElastic() ) 
    
    # Discretization
    domain = FEGrid( coord_max = (length,height,0.), 
                           shape   = shape,
                           fets_eval = fets_eval )
    
    middle_column = (shape[0]-1) / 2
    #domain.deactivate( ( middle_column, 0) )
    #
    bottom_dof = domain[middle_column,0,0,0].dofs[0,0,1]
    top_dof_slice = domain[middle_column,-1,:,-1]
    
    #bottom_middle_dofs, bottom_middle_dof_r = domain.get_bottom_middle_dofs()
    #bottom_middle_dof = bottom_middle_dofs[0,1]
    
    total_fracture_energy = RTraceDomainListInteg( name = 'Total fracture_energy' ,
                                                   var = 'fracture_energy',
                                                   record_on = 'update' )
    
    w_max = -.001
    
    tstepper = TS( dof_resultants = True,
                   sdomain = domain,
         bcond_list =  [ BCSlice( var = 'u', value = w_max, dims = [1], slice = top_dof_slice ),
                         BCDofGroup( var='u', value = 0., dims = [0,1],
                                  get_dof_method = domain.get_bottom_left_dofs ),
                         BCDofGroup( var='u', value = 0., dims = [1],
                                  get_dof_method = domain.get_bottom_right_dofs ) ],
         rtrace_list =  [
                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                 var_y = 'F_int', idx_y = top_dof_slice.dofs[0,0,1],
                                 var_x = 'U_k', idx_x = bottom_dof,
                                 record_on = 'update'),
                    total_fracture_energy,                                 
                  RTraceDomainListField(name = 'Displacement' ,
                  var = 'u', idx = 0,
                  record_on = 'update',
                  warp = True),
                  RTraceDomainListField(name = 'Fracture Energy' ,
                  var = 'fracture_energy', idx = 0,
                  record_on = 'update',
                  warp = True),
                  RTraceDomainListField(name = 'Strain' ,
                  var = 'eps_app', idx = 0,
                  record_on = 'update',
                  warp = False),
                  RTraceDomainListField(name = 'Stress' ,
                  var = 'sig_app', idx = 0,
                  record_on = 'update',
                  warp = False),                             
                ]
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = tstepper, tolerance = 1e-3,
                   tline  = TLine( min = 0.0,  step = 0.3, max = 1.001))                   
    
    tloop.eval()
    
    return tloop, total_fracture_energy.integ_val

def calculate_3D( length = 1., height = 0.3, shape = (21, 10) ):

    from ibvpy.fets.fets3D.fets3D8h import \
        FETS3D8H
    from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
        MATS3DMicroplaneDamage, PhiFnStrainSoftening
    
    h = length / float( shape[0] )

    mdm = MATS3DMicroplaneDamage( model_version = 'compliance',
                                  E = 34000.,
                                  nu = 0.25,
                                  phi_fn = PhiFnStrainSoftening( G_f = 0.001117,
                                                                 f_t = 2.8968,
                                                                 h = h ) )
    
    fets_eval = FETS3D8H( mats_eval = mdm ) 
    
    # Discretization
    domain = FEGrid( coord_max = (length,height,1.), 
                           shape   = shape + (1,) ,
                           fets_eval = fets_eval )
    
    middle_column = (shape[0]-1) / 2
    #domain.deactivate( ( middle_column, 0) )
    #
    bottom_dof = domain[middle_column,0,0,0,0,0].dofs[0,0,1]
    top_dof_slice = domain[middle_column,-1,:,:,-1,:]
    
    total_fracture_energy = RTraceDomainListInteg( name = 'Total fracture_energy' ,
                                                   var = 'fracture_energy',
                                                   record_on = 'update' )
    
    w_max = -.001
    
    tstepper = TS( dof_resultants = True,
                   sdomain = domain,
         bcond_list =  [ BCSlice( var = 'u', value = w_max, dims = [1], slice = top_dof_slice ),
                         BCDofGroup( var='u', value = 0., dims = [0,1,2],
                                  get_dof_method = domain.get_bottom_left_dofs ),
                         BCDofGroup( var='u', value = 0., dims = [1],
                                  get_dof_method = domain.get_bottom_right_dofs ) ],
         rtrace_list =  [
                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                 var_y = 'F_int', idx_y = top_dof_slice.dofs[0,0,1],
                                 var_x = 'U_k', idx_x = bottom_dof,
                                 record_on = 'update'),
                    total_fracture_energy,                                 
                 RTraceDomainListField(name = 'Fracture Energy' ,
                  var = 'fracture_energy', idx = 0,
                  record_on = 'update',
                  warp = True),
                  RTraceDomainListField(name = 'Strain' ,
                  var = 'eps_app', idx = 0,
                  record_on = 'update',
                  warp = False),
                  RTraceDomainListField(name = 'Stress' ,
                  var = 'sig_app', idx = 0,
                  record_on = 'update',
                  warp = False),                             
                ]
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = tstepper, tolerance = 1e-3,
#                   tline  = TLine( min = 0.0,  step = 0.5, max = 1.001))                   
                   tline  = TLine( min = 0.0,  step = 0.5, max = 1.001))                   
    
    tloop.eval()
    
    return tloop, total_fracture_energy.integ_val

if __name__ == '__main__':

    visual = True
    
    if visual:
        tloop, integ_val = calculate_3D( length = 1., height = 0.3, shape = (31, 10) )
        # Put the whole thing into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        #
        print integ_val
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = tloop )
        app.main()    
    else:
        elems = array( [7,15,21,27,33,39], dtype = int)        
        energies = zeros( elems.shape )
        for i, n_elem_x in enumerate( elems ):
            tloop, energy = create_model( length = 1., height = 0.3, shape = (n_elem_x, 10) )
            energies[i] = energy
            print 'shape' ,i , 'energy', energy
            
        print 'elemes'
        print elems
        print 'energies'
        print energies

        import matplotlib.pyplot as plt
        plt.plot(elems, energies)
        plt.title('fracture energy for different element sizes')
        plt.xlabel('number of elements in X direction')
        plt.ylabel('fracture energy')
        plt.show()
