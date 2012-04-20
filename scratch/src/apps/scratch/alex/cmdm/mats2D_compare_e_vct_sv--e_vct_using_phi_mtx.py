from numpy import *

from scipy.linalg import \
    eigh, inv
    
from ibvpy.mats.mats2D.mats2D_tensor import \
    map2d_eps_eng_to_mtx, map2d_sig_eng_to_mtx, map2d_eps_mtx_to_eng, map2d_sig_mtx_to_eng, \
    map2d_ijkl2mn, map2d_tns2_to_tns4, map2d_tns4_to_tns2, compliance_mapping2d, \
    get_D_plane_stress, get_D_plane_strain, get_C_plane_stress, get_C_plane_strain

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

# LS 20, t_max = 0,0005, alpha_degree = 0, sum-type symmetrization


n_mp = 7
alpha_list = array([ Pi / n_mp * (i - 0.5) for i in range(1,n_mp+1) ])
_MPN = array([[ cos( alpha ), sin( alpha )] for alpha in alpha_list ])

print '# damage tensor and damage effect tensor of 2nd order:'
phi_mtx = array([[  2.70090707e-01,   6.07153217e-17],
                 [  6.07153217e-17,   8.11599595e-01]])
print 'phi_mtx', phi_mtx
print '\n'

inv_phi_mtx = inv(phi_mtx)
print 'inv_phi_mtx', inv_phi_mtx
print '\n'

phi_arr  = array([0.130308983633, 0.262649045704, 1.0, 1.0, 1.0, 0.262649045704, 0.130308983633])
print 'phi_arr', phi_arr
print '\n'

psi_mtx  = array([[  5.60550123e+00,  -3.88578059e-16],
                  [ -3.88578059e-16,   1.81245877e+00]])
print 'psi_mtx', psi_mtx
print '\n'

psi_arr  = array([7.6740679892, 3.80736201542, 1.0, 1.0, 1.0, 3.80736201542, 7.6740679892])
print 'psi_arr', psi_arr
print '\n'

print '# apparent strain tensor:'
eps_app_eng = array([  4.75000000e-04,  -3.83961167e-05,  -2.96738868e-20])
print 'eps_app_eng', eps_app_eng
print '\n'

print '# projection of the strain tensor:'
eps_app_mtx = map2d_eps_eng_to_mtx(eps_app_eng)
e_vct_arr = array( [ dot( eps_app_mtx, mpn ) for mpn in _MPN ] )
print 'e_vct_arr', e_vct_arr
print '\n'

eps_eff_mtx = dot( phi_mtx, eps_app_mtx )


print '# microplane strain vector based on stiffness version:'
e_app_vct_arr_sv = array([[  4.63090758e-04,  -8.54393975e-06],
                          [  3.71369954e-04,  -2.39395872e-05],
                          [  2.06094776e-04,  -3.45937058e-05],
                          [  1.42484181e-20,  -3.83961167e-05],
                          [ -2.06094776e-04,  -3.45937058e-05],
                          [ -3.71369954e-04,  -2.39395872e-05],
                          [ -4.63090758e-04,  -8.54393975e-06]])
print 'e_app_vct_arr_sv', e_app_vct_arr_sv
print '\n'

print '# microplane strain vector based on compliance version - using phi_mtx(!):'
e_app_vct_arr_phi_mtx  = array([[  4.63090758e-04,  -8.54393975e-06],
                                [  3.71369954e-04,  -2.39395872e-05],
                                [  2.06094776e-04,  -3.45937058e-05],
                                [  3.35552504e-20,  -3.83961167e-05],
                                [ -2.06094776e-04,  -3.45937058e-05],
                                [ -3.71369954e-04,  -2.39395872e-05],
                                [ -4.63090758e-04,  -8.54393975e-06]])
print 'e_app_vct_arr_phi_mtx', e_app_vct_arr_phi_mtx
print '\n'

e_eff_vct_arr = array([[  8.26136218e-05,  -4.71400502e-06],
                       [  6.62509808e-05,  -1.32083485e-05],
                       [  3.67665205e-05,  -1.90866167e-05],
                       [  4.51759790e-21,  -2.11845463e-05],
                       [ -3.67665205e-05,  -1.90866167e-05],
                       [ -6.62509808e-05,  -1.32083485e-05],
                       [ -8.26136218e-05,  -4.71400502e-06]])
print 'e_eff_vct_arr', e_eff_vct_arr
print '\n'

e_app_vct_arr = array( [ dot( psi_mtx, e_eff_vct ) for e_eff_vct in e_eff_vct_arr ] )
print 'e_app_vct_arr', e_app_vct_arr
print '\n'
     

print '# microplane strain vector based on compliance version:'
e_app_vct_arr_cv = array([[  6.33982550e-04,  -3.61755950e-05],
                          [  2.52241468e-04,  -5.02889645e-05],
                          [  3.67665205e-05,  -1.90866167e-05],
                          [  4.51759790e-21,  -2.11845463e-05],
                          [ -3.67665205e-05,  -1.90866167e-05],
                          [ -2.52241468e-04,  -5.02889645e-05],
                          [ -6.33982550e-04,  -3.61755950e-05]])
print 'e_app_vct_arr_cv', e_app_vct_arr_cv
print '\n'

