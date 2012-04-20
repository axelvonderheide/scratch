
from ibvpy.mesh.cell_grid.cell_spec import CellSpec
from ibvpy.mesh.cell_grid.cell_grid import CellGrid
from numpy import frompyfunc, c_, array, cos, sin
from math import pi
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from ibvpy.fets.fets3D.fets3D8h27u import FETS3D8H27U
from ibvpy.mesh.fe_grid import FEGrid

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDofGroup, IBVPSolve as IS, DOTSEval

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic



def arch_3d( points ):
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    D_z = 3.
    R = 30.
    alpha = pi / 2.
    len_x = x[-1] - x[0]
    len_z = z[-1] - z[0]
    zz = z - len_z / 2.
    yy = y + 4 * D_z / len_z ** 2 * zz ** 2 + D_z

    phi_min, phi_max = pi / 2. - alpha / 2., pi / 2. + alpha / 2.
    delta_phi = phi_max - phi_min
    phi = phi_min + ( x / len_x ) * delta_phi
    r = R + yy
    #return c_[ x, yy, zz ]
    x, y, z = -r * cos( phi ), r * sin( phi ), zz * ( ( phi - delta_phi ) ** 2 + 2. )
    return c_[ x, y, z ]


if __name__ == '__main__':
    from ibvpy.plugins.ibvpy_app import IBVPyApp

    mats_eval = MATS3DElastic( E = 2.1e5, nu = 0.25 )

    fets_eval = FETS3D8H( mats_eval = mats_eval )
    fets_eval = FETS3D8H27U( mats_eval = mats_eval )

    fe_domain = FEGrid( fets_eval = fets_eval,
                           coord_max = ( 10., 2., 10. ),
#                           shape = ( 20, 3, 20 ),
#                           shape = ( 10, 3, 10 ),
                           shape = ( 4, 1, 4 ),
                           geo_transform = arch_3d )

    ts = TS( sdomain = fe_domain,
             bcond_list = [BCDofGroup( var = 'u', value = 0., dims = [0, 1, 2],
                                  get_dof_method = fe_domain.get_bottom_left_dofs ),
                           BCDofGroup( var = 'u', value = 0., dims = [1, 2],
                                  get_dof_method = fe_domain.get_bottom_right_dofs ),
                           BCDofGroup( var = 'f', value = -1., dims = [1],
                                  get_dof_method = fe_domain.get_top_dofs ),
#                           BCDofGroup( var='f', value = - 0.002, dims = [2],
#                                  get_dof_method = domain.get_left_dofs ) 
                    ],
             rtrace_list = [
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                  var_y = 'F_int', idx_y = right_dof,
#                                  var_x = 'U_k', idx_x = right_dof,
#                                  record_on = 'update'),
                        RTraceDomainListField( name = 'Deformation' ,
                                       var = 'eps_app', idx = 0,
                                       record_on = 'update' ),
                         RTraceDomainListField( name = 'Displacement' ,
                                        var = 'u', idx = 1,
                                        warp = False ),
                         RTraceDomainListField( name = 'Stress' ,
                                        var = 'sig_app', idx = 0, warp = True,
                                        record_on = 'update' ),
#                        RTraceDomainListField(name = 'N0' ,
#                                       var = 'N_mtx', idx = 0,
#                                       record_on = 'update')
                        ]
            )

    # Add the time-loop control
    tloop = TLoop( tstepper = ts, tolerance = 1e-5,
         tline = TLine( min = 0.0, step = 1.0, max = 1.0 ) )

    tloop.eval()

    ibvpy_app = IBVPyApp( ibv_resource = tloop )
    ibvpy_app.main()
