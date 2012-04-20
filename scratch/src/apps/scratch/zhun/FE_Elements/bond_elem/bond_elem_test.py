def example_with_new_domain():    
    from ibvpy.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, IBVPSolve as IS, DOTSEval
    from ibvpy.api import BCDofGroup
    from ibvpy.mats.mats1D5.mats1D5bond import MATS1D5Bond
    from ibvpy.mesh.mgrid_domain import MeshGridAdaptor
    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.fets.fets1D5.fets1D52b4uLRH import FETS1D52B4ULRH
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
    
    fets_eval = FETS1D52B4ULRH(mats_eval = MATS1D5Bond(Ef = 0., Af =0.,
                                                       Em = 0., Am = 0.,
                                                       bond_fn = MFnLineArray(xdata = [0.,1.],
                                                                              ydata = [0.,1.])))        

    # Tseval for a discretized line domain
    tseval  = DOTSEval( fets_eval = fets_eval )
 
    domain = FEGrid( coord_max = (1., 0.1, 0.0), #new domain
                           shape   = (2,1),
                           n_nodal_dofs = 1,
                           dof_r =  [[-1,-1],[1,-1],[1,1],[-1,1]],           
                           geo_r =  [[-1,-1],[1,-1],[1,1],[-1,1]])
    
                            
    ts = TS( tse = tseval,
         sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
         bcond_list =  [BCDofGroup(var='u', value = 0.,dims = [0],
                               get_dof_method = domain.get_top_dofs),
#                        BCDofGroup(var='u', value = 0.,dims = [0],
#                               get_dof_method = domain.get_bottom_left_dofs),
                      # imposed displacement for all right dofs in y-direction:
#                        BCDofGroup(var='f', value = -1., dims = [0],
#                                get_dof_method = domain.get_bottom_left_dofs),
                        BCDofGroup(var='f', value = 1./3., dims = [0],
                                get_dof_method = domain.get_bottom_dofs )],
         
                            
         rtrace_list =  [ RTraceGraph(name = 'Internal Force - Displacement' ,
                                var_y = 'F_int', idx_y = domain.get_bottom_right_dofs()[0][0,0],
                                var_x = 'U_k', idx_x = domain.get_bottom_right_dofs()[0][0,0]),
                        RTraceDomainField(name = 'Stress' ,
                                #position = 'int_pnts',
                                var = 'sig_app', idx = 0),
                        RTraceDomainField(name = 'Displacement' ,
                                var = 'u', idx = 0),
                        RTraceDomainField(name = 'Strain' ,
                                var = 'eps_app', idx = 0),
                        RTraceDomainField(name = 'Slip' ,
                                var = 'slip', idx = 0),
                        RTraceDomainField(name = 'Shear' ,
                                var = 'shear', idx = 0)       
                                
#                             RTraceDomainField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')
                      
                ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
             DT = 1.,
             tline  = TLine( min = 0.0,  max = 1.0 ))
    
    tloop.eval()    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()


if __name__ == '__main__':
    example_with_new_domain()