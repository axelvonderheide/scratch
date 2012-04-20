'''
Example of a tensile test using a one element discretization
'''
if __name__ == '__main__':
    
    from ibvpy.api import \
        TStepper as TS, TLoop, TLine, \
        IBVPSolve as IS, DOTSEval, \
        RTraceGraph, RTraceDomainListField,  \
        BCDof, BCDofGroup,  BCSlice   
    
    
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q
    
    from ibvpy.mesh.fe_grid import FEGrid
    
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
        MATS2DElastic
        
    from mathkit.geo.geo_ndgrid import \
        GeoNDGrid
    from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
        MFnNDGrid, GridPoint
    
    from numpy import sin, cos, c_, arange, hstack, array, loadtxt, size
    
    from time import time
    from os.path import join
    
    from math import \
        pi as Pi, cos, sin, exp, sqrt as scalar_sqrt
    
    
    
    mats = MATS2DElastic( E = 35000., nu = 0.2, stress_state  = "plane_strain" )
    
    fets_eval = FETS2D4Q( mats_eval = mats ) 
    
    a=2.0
    
    
    
    domain = FEGrid( coord_max = (1., 1., 0.), 
                     shape   = (a, a),
                     fets_eval = fets_eval )
    
    loading_slice = domain[-1,:,-1,:]
    
    
    
    ts = TS(
            sdomain = domain,
    
    #        # simple shear test: clamped at left side loaded at right side
    #        bcond_list = [BCSlice( var = 'u', value = 0., dims = [0,1], slice = domain[0, 0, 0, :] ),
    #                      BCSlice( var = 'u', value = 0., dims = [0]  , slice = domain[0, 0,-1, :] ),
    #                      BCSlice( var = 'f', value = 1.0, dims = [1] , slice = domain[0, 0,-1, :] )],
    
            # shear test: fixed at 000, one support in x-direction at left top; load at top right
            bcond_list = [BCSlice( var = 'u', value = 0., dims = [0,1], slice = domain[0, :, 0, :] ),
                          BCSlice( var = 'f', value = 1/a, dims = [0], slice = loading_slice )], # for 2 MN force on right side udl
    
             
             
             rtrace_list = [ 
                            
                         RTraceDomainListField(name = 'Displacement' ,
                                        var = 'u', idx = 0, warp = True),
                         RTraceDomainListField(name = 'Stress' ,
                                        var = 'sig_app', idx = 0, warp = True, 
                                        record_on = 'update'),
                        ]             
            )
     
      
                 
    tloop = TLoop( tstepper = ts,
                   #tolerance = 1e-3,
                   tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))
    
    u = tloop.eval()
    
    #print 'loading dofs', loading_slice.dofs
    print u[ loading_slice.dofs ]
    print size(u[ loading_slice.dofs ])
    print 'x, y -displacement (right side top)', u[ loading_slice.dofs ][-1][-1]
    print 'x -displacement  (right side, middle)', u[ loading_slice.dofs ][len(loading_slice.dofs)/2][0][0]
      
    ui = False
    
    if ui:
            
            # Put the whole stuff into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = tloop )
        app.main()

