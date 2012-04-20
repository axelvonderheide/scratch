

from numpy import tensordot, dot, zeros

from v_get_B_grid import get_B_mtx_grid

def get_D_mtx():
    E = 30e+5
    nu = 0.3
    D_stress = zeros([3,3])
    D_stress[0,0] = E/(1.0-nu*nu)
    D_stress[0,1] = E/(1.0-nu*nu)*nu
    D_stress[1,0] = E/(1.0-nu*nu)*nu
    D_stress[1,1] = E/(1.0-nu*nu)
    D_stress[2,2] = E/(1.0-nu*nu)*(1.0/2.0-nu/2.0)
    return D_stress

D_grid = get_D_mtx()#[None,None,:,:]

from ibvpy.mesh.fe_grid import FEGrid
# Discretization
domain = FEGrid( coord_min = (0., -1.,0.),
                       coord_max = (1.,-0.,0.), 
                       shape   = (2,1),
                       n_nodal_dofs = 2,
                       dof_r = [[-1,-1],[1,-1],[1,1],[-1,1]],
                       geo_r = [[-1,-1],[1,-1],[1,1],[-1,1]] )

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
fets_eval = FETS2D4Q()
B_grid = get_B_mtx_grid( domain, fets_eval )    

# DB_grid[ sig, el, ip, el_dof] = D_grid[ sig, eps ] * B_grid[ el, ip, eps, el_dof ]

DB_grid = dot( D_grid, B_grid )
print DB_grid.shape

# B_grid.t[ el, ip, el_dof, eps ] * DB_grid[ sig, el, ip, el_dof ] 

BtDB_grid = tensordot( B_grid, DB_grid, axes = (2,1) )

#print DtDB_grid