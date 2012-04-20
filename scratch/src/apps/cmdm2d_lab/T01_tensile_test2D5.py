'''
Example of a tensile test using a one element discretization
'''

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDofGroup, IBVPSolve as IS, DOTSEval, BCSlice

from ibvpy.mesh.fe_grid import FEGrid
    
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets2D5.fets2D58h import FETS2D58H
from ibvpy.util.simgrid import simgrid
from math import fabs        

fets_eval3D  = FETS3D8H(mats_eval  = MATS3DElastic(E = 34000, nu = 0.25))
fets_eval2D5 = FETS2D58H(mats_eval = MATS2DElastic(E = 34000, nu = 0.25, 
                                                        stress_state  = "plane_strain"))


support_slices = [
                  [ (-1   ,slice(None),slice(None),-1   ,slice(None),slice(None)), # yz plane  0
                    (slice(None) ,-1  ,slice(None),slice(None) ,-1  ,slice(None)), #  z-axis   1
                    (0   ,0   ,   0,0   ,0   ,0   )  #  origin   2
                  ]
                  ]


support_dirs = [[0],[1],[2]]

loading_slices = [ 
                  (-1  ,-1, -1, -1, -1, -1 ),  # loading in z dir
                ]

load_dirs = [2]
load = - 0.00025 # - 0.01
vars = ['eps_app']
for support_slice, loading_slice in zip( support_slices, loading_slices ): 

    for load_dir in load_dirs:
        tl, u3D, fields3D, integs, g  = simgrid( fets_eval3D, (0.5, 0.5, .03), (5,5,1),
                                   support_slice, support_dirs,
                                   loading_slice, load_dir, 
                                   load, 1, vars, var_type = 'f' )

        tl, u2D5, fields2D5, integs, g = simgrid( fets_eval2D5, (0.5, 0.5, .03), (5,5,1),
                                   support_slice, support_dirs,
                                   loading_slice, load_dir, 
                                   load, 1, vars, var_type = 'f' )
        
        
        print 'u215', u3D[215]
        
        for u1_, u2_ in zip( u3D.flatten(), u2D5.flatten() ):
            if fabs( u1_ - u2_ ) > 1e-10:
                print 'non matching', u1_, u2_

