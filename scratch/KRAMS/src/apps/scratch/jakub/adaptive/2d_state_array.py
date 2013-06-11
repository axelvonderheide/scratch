'''
Created on Mar 30, 2011

@author: jakub
'''
if __name__ == '__main__':
    def example_2d():
        from ibvpy.api import FEDomain, FERefinementGrid, FEGrid, TStepper as TS, \
            BCDofGroup, RTraceDomainListField
        from ibvpy.core.tloop import TLoop, TLine
        from ibvpy.mesh.xfe_subdomain import XFESubDomain
        from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
        from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
        from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
        from ibvpy.fets.fets2D.fets2D4q9u import FETS2D4Q9U
        from ibvpy.fets.fets2D.fets2D4q12u import FETS2D4Q12U
        from ibvpy.fets.fets2D.fets2D4q16u import FETS2D4Q16U
        from ibvpy.fets.fets_ls.fets_crack import FETSCrack

        fets_eval = FETS2D4Q( mats_eval = MATS2DScalarDamage( E = 1., nu = 0. ) )

        xfets_eval = FETSCrack( parent_fets = fets_eval, int_order = 1 )

        # Discretization

        fe_domain = FEDomain()
        fe_level1 = FERefinementGrid( domain = fe_domain,
                                      fets_eval = fets_eval )
        fe_grid1 = FEGrid( coord_max = ( 2., 1., 0. ),
                           shape = ( 2, 1 ),
                           fets_eval = fets_eval,
                           level = fe_level1 )


        fe_xdomain = XFESubDomain( domain = fe_domain,
                                    rt_quad = False,
                                    fets_eval = xfets_eval,
                                    #fe_grid_slice = fe_grid1['Y - 0.5@ X < .5'] 
                                    fe_grid_slice = fe_grid1['X - 0.5']
                                    )

        fe_xdomain.deactivate_sliced_elems()

        ts = TS( dof_resultants = True,
                 sdomain = fe_domain
                    )



        ts.setup()
#        print 'parent elems ', fe_xdomain.fe_grid_slice.elems
#        print 'intersection points ', fe_xdomain.fe_grid_slice.r_i
#        print 'ip_coords ', fe_xdomain.dots.ip_coords
        print 'state array step 1 ', fe_xdomain.dots.state_array
        fe_xdomain.dots.state_array[:3] = [1, 2, 3]
        print 'state array write ', fe_xdomain.dots.state_array
        fe_xdomain.fe_grid_slice = fe_grid1['Y-0.5']
#        print 'parent elems ', fe_xdomain.fe_grid_slice.elems
#        print 'intersection points ', fe_xdomain.fe_grid_slice.r_i
#        print 'ip_coords ', fe_xdomain.dots.ip_coords
        print 'state array step 2 ', fe_xdomain.dots.state_array

    example_2d()
