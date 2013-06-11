'''
Created on Sep 20, 2012

The CompositeCrackBridge class has a method for evaluating fibers and matrix
strain in the vicinity of a crack bridge.
Fiber diameter and bond coefficient can be set as random variables.
Reinforcement types can be combined by creating a list of Reinforcement
instances and defining it as the reinforcement_lst Trait in the
CompositeCrackBridge class.
The evaluation is array based.

@author: rostar
'''
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
import numpy as np
from stats.spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array
from types import FloatType
from reinforcement import Reinforcement, WeibullFibers
from scipy.optimize import fsolve, broyden2
import time as t
from scipy.integrate import cumtrapz
import time
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from math import pi
from scipy.interpolate import interp1d

def H( x ):
        return x > 0


class CompositeCrackBridge( HasTraits ):

    reinforcement_lst = List( Instance( Reinforcement ) )
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float

    V_f_tot = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_V_f_tot( self ):
        V_f_tot = 0.0
        for reinf in self.reinforcement_lst:
            V_f_tot += reinf.V_f
        return V_f_tot

    E_c = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_E_c( self ):
        E_fibers = 0.0
        for reinf in self.reinforcement_lst:
            E_fibers += reinf.V_f * reinf.E_f
        return self.E_m * ( 1. - self.V_f_tot ) + E_fibers
    
    grid_ones = Property()
    @cached_property
    def _get_grid_ones( self ):
        '''symmetric, quadratic matrix used in perc'''
        return np.ones( ( len( self.sorted_depsf ), len( self.sorted_depsf ) ) )
    

    sorted_theta = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_theta( self ):
        '''sorts the integral points by bond in descending order'''
        depsf_arr = np.array( [] )
        V_f_arr = np.array( [] )
        E_f_arr = np.array( [] )
        xi_arr = np.array( [] )
        stat_weights_arr = np.array( [] )
        nu_r_arr = np.array( [] )
        lf_arr = np.array( [] )
        phi_arr = np.array( [] )
        for reinf in self.reinforcement_lst:
            n_int = len( np.hstack( ( np.array( [] ), reinf.depsf_arr ) ) )
            depsf_arr = np.hstack( ( depsf_arr, reinf.depsf_arr ) )
            phi_arr = np.hstack( ( phi_arr, reinf.phi_arr ) )
            V_f_arr = np.hstack( ( V_f_arr, np.repeat( reinf.V_f, n_int ) ) )
            E_f_arr = np.hstack( ( E_f_arr, np.repeat( reinf.E_f, n_int ) ) )
            xi_arr = np.hstack( ( xi_arr, np.repeat( reinf.xi, n_int ) ) )
            stat_weights_arr = np.hstack( ( stat_weights_arr,
                                          np.repeat( reinf.stat_weights, n_int ) ) )
            nu_r_arr = np.hstack( ( nu_r_arr, reinf.nu_r ) )
            lf_arr = np.hstack( ( lf_arr, np.repeat( reinf.l_f, n_int ) ) )
        argsort = np.argsort( depsf_arr )[::-1]
        idxs = np.array( [] )
        for i, reinf in enumerate( self.reinforcement_lst ):
            idxs = np.hstack( ( idxs, i * np.ones_like( reinf.depsf_arr ) ) )
        masks = []
        for i, reinf in enumerate( self.reinforcement_lst ):
            masks.append( ( idxs == i )[argsort] )
        return depsf_arr[argsort], V_f_arr[argsort], E_f_arr[argsort], \
                xi_arr[argsort], stat_weights_arr[argsort], \
                nu_r_arr[argsort], masks, lf_arr[argsort], phi_arr[argsort]

    sorted_depsf = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_depsf( self ):
        return self.sorted_theta[0]

    sorted_V_f = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_V_f( self ):
        return self.sorted_theta[1]

    sorted_E_f = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_E_f( self ):
        return self.sorted_theta[2]

    sorted_xi = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_xi( self ):
        return self.sorted_theta[3]

    sorted_stats_weights = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_stats_weights( self ):
        return self.sorted_theta[4]

    sorted_nu_r = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_nu_r( self ):
        return self.sorted_theta[5]

    sorted_masks = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_masks( self ):
        return self.sorted_theta[6]
    
    sorted_lf = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_lf( self ):
        return self.sorted_theta[7]
    
    sorted_phi = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_phi( self ):
        return self.sorted_theta[8]

    sorted_xi_cdf = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_xi_cdf( self ):
        '''breaking strain: CDF for random and Heaviside for discrete values'''
        # TODO: does not work for reinforcement types with the same xi
        methods = []
        masks = []
        for reinf in self.reinforcement_lst:
            masks.append( self.sorted_xi == reinf.xi )
            if isinstance( reinf.xi, FloatType ):
                methods.append( lambda x: 1.0 * ( reinf.xi <= x ) )
            elif isinstance( reinf.xi, RV ):
                methods.append( reinf.xi._distr.cdf )
            elif isinstance( reinf.xi, WeibullFibers ):
                methods.append( reinf.xi.weibull_fibers_Pf )
        return methods, masks

    def vect_xi_cdf( self, epsy, x_short, x_long ):
        Pf = np.zeros_like( self.sorted_depsf )
        methods, masks = self.sorted_xi_cdf
        for i, method in enumerate( methods ):
            if method.__name__ == 'weibull_fibers_Pf':
                Pf += method( epsy * masks[i], self.sorted_depsf,
                             x_short = x_short, x_long = x_long )
            else:
                Pf += method( epsy * masks[i] )
        return Pf
    
    def perc( self, a , Kf, damage, depsf ):
        '''evaluates muT and Kf for every point in a '''
        #transforming Kf and damage into matrix array with length of a rows
        Kf = Kf * self.grid_ones
        damage = damage * self.grid_ones
        #getting lf of fibers
        lf = self.sorted_lf
        #reshaping a and phi to get an evaluation in every a
        a = a.reshape( len( a ), 1 )
        phi = self.sorted_phi.reshape( 1, len( self.sorted_phi ) )
        #geometrical condition evaluated in every a
        m_p = a * 2 / np.cos( phi ) / lf  
        # continious fibers are flagged with a negative value in m_p- The following mask is to filter the array.
        mask1 = m_p > 0
        m_p = m_p * mask1
        #m_p is a matrix array giving for every phi in a the percentage of fibers that are geometrically dropped out. 
        #p is the percentage of active fibers. Heaviside cuts all values exceeding 100% drop out.
        p = np.abs( H( 1 - m_p ) * ( 1 - m_p ) )
        #creating matrix that implements the mechanical condition. 
        #In every  step in 'a' there is one integral point of bond less active.
        #So it is a diagonal matrix with ones on the right side and zeroes on the other.
        stepm = np.eye( len( a ), len( phi[0] ) )
        eym = np.cumsum( stepm , -1 )
        #summing up all the matrix arrays in the axis of 'a' for muT and Kf
        muT = np.sum( depsf * ( 1 - damage ) * Kf * p * eym, 1 )
        Kf = np.sum( ( 1 - damage ) * Kf * p * eym, 1 )[::-1]
        return muT, Kf

    Emtrx = Array
    def dem_depsf_vect( self, depsf, damage, a ):
        '''evaluates the deps_m given deps_f
        at that point and the damage array'''
        Kf = self.sorted_V_f * self.sorted_nu_r * \
            self.sorted_stats_weights * self.sorted_E_f
        Kf_broken = np.sum( Kf * damage )
        Km = ( 1. - self.V_f_tot ) * self.E_m
        #in case of first iteration or only continous fibers, the first if clause evaluates depsm
        if a[0] == 1e-10 or self.sorted_lf.all() < 0:
            mu_T = np.cumsum( ( depsf * Kf * ( 1. - damage ) )[::-1] )[::-1]
            Kf_intact_bonded = np.hstack( ( 0.0, np.cumsum( ( Kf * ( 1. - damage ) ) ) ) )[:-1]
            Kf_add = Kf_intact_bonded + Kf_broken
            E_mtrx = Km + Kf_add
        #every further iteration of composites with short fibers is evaluated with perc function.
        else:
            mu_T, Kf_intact_bonded = self.perc( a[( len( a ) - 1 ) / 2 + 1:-1] , Kf, damage, depsf )
            Kf_add = Kf_intact_bonded + Kf_broken
            E_mtrx = Km + Kf_add
        #print E_mtrx, Kf_add
        self.Emtrx = E_mtrx
        return mu_T / E_mtrx

    def F( self, dems, amin ):
        '''
        
        #
        deps_min = []
        deps_max = []
        f_chi_list = []
        for mask in self.sorted_masks :
            min_max_mask = mask[::-1]
            deps_min.append( np.nonzero( min_max_mask )[0][0] )
            deps_max.append( np.nonzero( min_max_mask )[0][-1] )
            #####
            depsfi = self.sorted_depsf[mask]
            demsi = dems[mask]
            fi = 1. / ( depsfi + demsi )
            f_chi_list.append( interp1d( depsfi[::-1], demsi[::-1], bounds_error = False, fill_value = 0 ) )
            
        argsort = np.argsort( deps_min )
        deps_min_sorted = np.array( deps_min )[argsort]
        deps_max_sorted = np.array( deps_max )[argsort]
        turned_depsf = self.sorted_depsf[::-1]
        integ_deps_list = list( turned_depsf )
        cf = 0
        for i, max in enumerate( deps_max_sorted[:-1] ):
            if max < deps_min_sorted[i + 1]:
                integ_deps_list.insert( max + 1 + cf, integ_deps_list[ max + cf  ] + 1e-15 )
                integ_deps_list.insert( max + 2 + cf, integ_deps_list[ max + 2 + cf ] - 1e-15 )
                cf += 2

        integ_deps_arr = np.array( integ_deps_list )
        
        fi_integ = np.zeros_like( integ_deps_arr )
        for func in f_chi_list:
            fi_integ += func( integ_deps_arr )
            #print f_integ
        
        F_ix = cumtrapz( fi_integ[::-1], -integ_deps_arr[::-1], initial = 0 ) 
        F_func = interp1d( -integ_deps_arr[::-1], F_ix )
        F_test = F_func( -self.sorted_depsf )
        #plt.plot( -integ_deps_arr[::-1], fi_integ[::-1] )
        # plt.plot( -integ_deps_arr[::-1], F_ix )
        #plt.show()
        
        
        
        #######
        depsfi_rr = np.array( [] )
        fi_rr = np.array( [] )
        depsfi_list = []
        fi_list = []
        '''
        F = np.zeros_like( self.sorted_depsf )
        for i, mask in enumerate( self.sorted_masks ):
            depsfi = self.sorted_depsf[mask]
            demsi = dems[mask]
            fi = 1. / ( depsfi + demsi )
            #fi_list.append( fi )
            #fi_rr = np.hstack( ( fi_rr, fi[::-1] ) )
            #depsfi_rr = np.hstack( ( depsfi_rr, depsfi[::-1] ) )
            #depsfi_list.append( depsfi )
            F[mask] = np.hstack( ( cumtrapz( fi, -depsfi, initial = 0 ) ) )

            if i == 0:
                C = 0.0
            else:
                depsf0 = self.sorted_depsf[self.sorted_masks[i - 1]]
                depsf1 = depsfi[0]
                idx = np.sum( depsf0 > depsf1 ) - 1
                depsf2 = depsf0[idx]
                a1 = np.exp( F[self.sorted_masks[i - 1]][idx] / 2. + np.log( amin ) )
                p = depsf2 - depsf1
                q = depsf1 + demsi[0]
                amin_i = np.sqrt( a1 ** 2 + p / q * a1 ** 2 )
                C = np.log( amin_i / amin )
            F[mask] += 2 * C
        #print depsfi_rr
        
        #plt.plot( depsfi_rr, fi_rr )
        #plt.plot( self.sorted_depsf[::-1], F[::-1] , 'ro' )
        #zeros = np.zeros_like( self.sorted_depsf[::-1] )
        #plt.plot( self.sorted_depsf[::-1], zeros, 'ro' )
        #for i, x in enumerate( depsfi_list ):
        #    plt.plot( x, fi_list[i] )
        #fi_arr = np.array( fi_list )
        #plt.show()
        #dfb = MFnLineArray( xdata = self.sorted_depsf[::-1], ydata = F[::-1] )
        #plt.plot( self.sorted_depsf[::-1], dfb.get_diffs( self.sorted_depsf[::-1] ) )
        #plt.plot( -self.sorted_depsf, F )
       # print len( self.sorted_depsf ), len( fi_list )
        #plt.plot( self.sorted_depsf, fi_list )
        #plt.show()
        return F

    def profile( self, iter_damage, Lmin, Lmax ):
        # matrix strain derivative with resp. to z as a function of T
        dems = self.dem_depsf_vect( self.sorted_depsf, iter_damage, self._x_arr )
        # initial matrix strain derivative
        init_dem = dems[0]
        # debonded length of fibers with Tmax
        amin = ( self.w / ( np.abs( init_dem ) + np.abs( self.sorted_depsf[0] ) ) ) ** 0.5
        # integrated f(depsf) - see article
        F = self.F( dems, amin )
        # a(T) for double sided pullout
        a1 = np.exp( F / 2. + np.log( amin ) )
        if Lmin < a1[0] and Lmax < a1[0]:
            # all fibers debonded up to Lmin and Lmax
            a = np.hstack( ( -Lmin, 0.0, Lmax ) )
            em = np.hstack( ( init_dem * Lmin, 0.0, init_dem * Lmax ) )
            epsf0 = ( self.sorted_depsf / 2. * ( Lmin ** 2 + Lmax ** 2 ) + 
                     self.w + em[0] * Lmin / 2. + em[-1] * Lmax / 2. ) / ( Lmin + Lmax )
        elif Lmin < a1[0] and Lmax >= a1[0]:
            # all fibers debonded up to Lmin but not up to Lmax
            amin = -Lmin + np.sqrt( 2 * Lmin ** 2 + 2 * self.w / ( self.sorted_depsf[0] + init_dem ) )
            C = np.log( amin ** 2 + 2 * Lmin * amin - Lmin ** 2 )
            a2 = np.sqrt( 2 * Lmin ** 2 + np.exp( ( F + C ) ) ) - Lmin
            if Lmax <= a2[-1]:
                idx = np.sum( a2 < Lmax ) - 1
                a = np.hstack( ( -Lmin, 0.0, a2[:idx + 1], Lmax ) )
                em2 = np.cumsum( np.diff( np.hstack( ( 0.0, a2 ) ) ) * dems )
                em = np.hstack( ( init_dem * Lmin, 0.0, em2[:idx + 1], em2[idx] + ( Lmax - a2[idx] ) * dems[idx] ) )
                um = np.trapz( em, a )
                epsf01 = em2[:idx + 1] + a2[:idx + 1] * self.sorted_depsf[:idx + 1]
                epsf02 = ( self.w + um + self.sorted_depsf[idx + 1:] / 2. * ( Lmin ** 2 + Lmax ** 2 ) ) / ( Lmin + Lmax )
                epsf0 = np.hstack( ( epsf01, epsf02 ) )
            else:
                a = np.hstack( ( -Lmin, 0.0, a2, Lmax ) )
                em2 = np.cumsum( np.diff( np.hstack( ( 0.0, a2 ) ) ) * dems )
                em = np.hstack( ( init_dem * Lmin, 0.0, em2, em2[-1] ) )
                epsf0 = em2 + self.sorted_depsf * a2
        elif a1[0] < Lmin and a1[-1] > Lmin:
            # boundary condition position
            idx1 = np.sum( a1 <= Lmin )
            # a(T) for one sided pullout
            # first debonded length amin for one sided PO
            depsfLmin = self.sorted_depsf[idx1]
            p = ( depsfLmin + dems[idx1] )
            a_short = np.hstack( ( a1[:idx1], Lmin ) )
            em_short = np.cumsum( np.diff( np.hstack( ( 0.0, a_short ) ) ) * dems[:idx1 + 1] )
            emLmin = em_short[-1]
            umLmin = np.trapz( np.hstack( ( 0.0, em_short ) ), np.hstack( ( 0.0, a_short ) ) )
            amin = -Lmin + np.sqrt( 4 * Lmin ** 2 * p ** 2 - 4 * p * emLmin * Lmin + 4 * p * umLmin - 2 * p * Lmin ** 2 * depsfLmin + 2 * p * self.w ) / p
            C = np.log( amin ** 2 + 2 * amin * Lmin - Lmin ** 2 )
            a2 = ( np.sqrt( 2 * Lmin ** 2 + np.exp( F + C - F[idx1] ) ) - Lmin )[idx1:]
            # matrix strain profiles - shorter side
            a_short = np.hstack( ( -Lmin, -a1[:idx1][::-1], 0.0 ) )
            dems_short = np.hstack( ( dems[:idx1], dems[idx1] ) )
            em_short = np.hstack( ( 0.0, np.cumsum( np.diff( -a_short[::-1] ) * dems_short ) ) )[::-1]
            if a2[-1] > Lmax:
                idx2 = np.sum( a2 <= Lmax )
                # matrix strain profiles - longer side
                a_long = np.hstack( ( a1[:idx1], a2[:idx2] ) )
                em_long = np.cumsum( np.diff( np.hstack( ( 0.0, a_long ) ) ) * dems[:idx1 + idx2] )
                a = np.hstack( ( a_short, a_long, Lmax ) )
                em = np.hstack( ( em_short, em_long, em_long[-1] + ( Lmax - a_long[-1] ) * dems[idx1 + idx2] ) )
                um = np.trapz( em, a )
                epsf01 = em_long + a_long * self.sorted_depsf[:idx1 + idx2]
                epsf02 = ( self.w + um + self.sorted_depsf [idx1 + idx2:] / 2. * ( Lmin ** 2 + Lmax ** 2 ) ) / ( Lmin + Lmax )
                epsf0 = np.hstack( ( epsf01, epsf02 ) )
            else:
                a_long = np.hstack( ( 0.0, a1[:idx1], a2, Lmax ) )
                a = np.hstack( ( a_short, a_long ) )
                dems_long = np.hstack( ( dems, dems[-1] ) )
                em_long = np.hstack( ( 0.0, np.cumsum( np.diff( a_long ) * dems_long ) ) )
                em = np.hstack( ( em_short, em_long ) )
                epsf0 = em_long[1:-1] + self.sorted_depsf * a_long[1:-1]
        elif a1[-1] <= Lmin:
            #double sided pullout
            a = np.hstack( ( -Lmin, -a1[::-1], 0.0, a1, Lmin ) )
            em1 = np.cumsum( np.diff( np.hstack( ( 0.0, a1 ) ) ) * dems )
            em = np.hstack( ( em1[-1], em1[::-1], 0.0, em1, em1[-1] ) )
            epsf0 = em1 + self.sorted_depsf * a1
    
        self._x_arr = a
        self._epsm_arr = em
        self._epsf0_arr = epsf0
        a_short = -a[a < 0.0][1:][::-1]
        if len( a_short ) < len( self.sorted_depsf ):
            a_short = np.hstack( ( a_short, Lmin * np.ones( len( self.sorted_depsf ) - len( a_short ) ) ) )
        a_long = a[a > 0.0][:-1]
        if len( a_long ) < len( self.sorted_depsf ):
            a_long = np.hstack( ( a_long, Lmax * np.ones( len( self.sorted_depsf ) - len( a_long ) ) ) )
        return epsf0, a_short, a_long

    def damage_residuum( self, iter_damage ):
        Lmin = min( self.Ll, self.Lr )
        Lmax = max( self.Ll, self.Lr )
        epsf0, x_short, x_long = self.profile( iter_damage, Lmin, Lmax )
        residuum = self.vect_xi_cdf( epsf0, x_short = x_short, x_long = x_long ) - iter_damage
        #print residuum
        return np.abs( residuum ) * -1

    _x_arr = Array
    def __x_arr_default( self ):
        return np.repeat( 1e-10, len( self.sorted_depsf ) )

    _epsm_arr = Array
    def __epsm_arr_default( self ):
        return np.repeat( 1e-10, len( self.sorted_depsf ) )

    _epsf0_arr = Array
    def __epsf0_arr_default( self ):
        return np.repeat( 1e-10, len( self.sorted_depsf ) )

    damage = Property( depends_on = 'w, Ll, Lr, reinforcement+' )
    @cached_property
    def _get_damage( self ):
        ff = time.clock()
        if self.w == 0.:
            damage = np.zeros_like( self.sorted_depsf )
        else:
            ff = t.clock()
            try:
                damage = broyden2( self.damage_residuum, 0.2 * np.ones_like( self.sorted_depsf ) * ( self.sorted_lf < 0 ), iter = 20 )#maxiter = 20 , )
            except:
                print 'broyden2 does not converge fast enough: switched to fsolve for this step'
                damage = fsolve( self.damage_residuum, 0.2 * np.ones_like( self.sorted_depsf ) * ( self.sorted_lf < 0 ) )
                
            cont_fibers = self.sorted_lf < 0 
            dam_cont = np.sum( damage[cont_fibers] )
            num_cont = np.sum( cont_fibers )
            if num_cont > 0:
                print 'damage of continiuos fibers =', dam_cont / num_cont ,
            print'iteration time =', time.clock() - ff, 'sec'
        return damage

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    reinfglass = Reinforcement( r = 0.003, #RV('uniform', loc=0.001, scale=0.005),
                          tau = RV( 'uniform', loc = 4., scale = 2. ),
                          V_f = 0.2,
                          E_f = 70e3,
                          xi = RV( 'weibull_min', shape = 5., scale = 0.04 ),
                          n_int = 30,
                          label = 'AR glass' )

    reinf = Reinforcement( r = 0.001, #RV('uniform', loc=0.002, scale=0.002),
                          tau = RV( 'uniform', loc = .2, scale = .3 ),
                          V_f = 0.2,
                          E_f = 200e3,
                          xi = WeibullFibers( shape = 100., scale = 0.02 ),
                          n_int = 100,
                          label = 'carbon' )
    
    reinfimage = Reinforcement( r = 0.001, #RV('uniform', loc=0.002, scale=0.002),
                          tau = RV( 'uniform', loc = 3., scale = 1000 ),
                          V_f = 0.01,
                          E_f = 200e3,
                          xi = WeibullFibers( shape = 100., scale = 0.02 ),
                          n_int = 100,
                          label = 'image' )
    

    
    reinf2 = Reinforcement( r = 0.1,
                          tau = 7. ,
                          lf = 5.,
                          snub = .7,
                          phi = RV( 'sin2x', loc = 0., scale = 1. ),
                          V_f = 0.03,
                          E_f = 200e3,
                          xi = WeibullFibers( shape = 100., scale = 100 ),
                          n_int = 100,
                          label = 'Short Fibers' )
    

    ccb = CompositeCrackBridge( E_m = 25e3,
                                 reinforcement_lst = [   reinf2, reinf  ],
                                 Ll = 20.,
                                 Lr = 20.,
                                 w = 0.01 )

    ccb.damage
    plt.plot( ccb._x_arr, ccb._epsm_arr, lw = 2, color = 'red', ls = 'dashed', label = 'analytical' )
   
    print ccb._epsm_arr[-1]
    print 0.000506374804145
    #plt.plot( np.zeros_like( ccb._epsf0_arr ), ccb._epsf0_arr, 'ro' )
   
    em_func = interp1d( ccb._x_arr, ccb._epsm_arr )
    ef_discr = np.linspace( 0, np.max( ccb.sorted_lf ) / 2., 1000 )
    em_ef = em_func( ef_discr )
    phi = ccb.sorted_phi
    if 1 == 1:
        res_list = []
        sigma_list = []
        for i, depsf in enumerate( ccb.sorted_depsf ):
            if ccb.sorted_lf[i] < 0:
                if i % 2 == 0:
                    plt.plot( ccb._x_arr, np.maximum( ccb._epsf0_arr[i] - depsf * np.abs( ccb._x_arr ), ccb._epsm_arr ), 'r' )
                else:
                    plt.plot( ccb._x_arr, np.maximum( ccb._epsf0_arr[i] - depsf * np.abs( ccb._x_arr ), ccb._epsm_arr ), 'k' )
            else:  
                a = ccb._x_arr[( len( ccb._x_arr ) - 1 ) / 2 + 1:-1]
                l = np.cos( ccb.sorted_phi[i] ) * ccb.sorted_lf[i] / 2.
                p = ( l - a[i] ) / l
                if p < 0:
                    al_discr = np.linspace( 0, l, 100 )
                    ems = em_func( al_discr ).reshape( len( al_discr ) , 1 )
                    al_y = al_discr.reshape( 1, len( al_discr ) )
                    al_x = al_discr.reshape( len( al_discr ) , 1 )
                    mmm = np.maximum( depsf * al_y + ems - depsf * al_x, ems )
                    E_rest = np.sum( mmm, 1 ) / len( al_discr )
                    Eeps = E_rest
                    #plt.plot( al_x, Eeps, 'k' )
                   
                else:
                    al_discr = np.linspace( 0, a[i], 100 )
                    ems = em_func( al_discr ).reshape( len( al_discr ) , 1 )
                    al_y = al_discr.reshape( 1, len( al_discr ) )
                    al_x = al_discr.reshape( len( al_discr ) , 1 )
                    mmm = np.maximum( depsf * al_y + ems - depsf * al_x, ems )
                    E_rest = np.sum( mmm, 1 ) / len( al_discr )
                    normeps = ccb._epsf0_arr[i] - depsf * np.abs( al_discr )
                    Eeps = ( p * normeps + ( 1 - p ) * E_rest ) 
                    #plt.plot( al_x, Eeps, 'y' )
                ef_func = interp1d( al_discr, Eeps, bounds_error = False, fill_value = 0 )
                ef_array = np.maximum( ef_func( ef_discr ), em_ef )
                plt.plot( ef_discr, ef_array )
                res_list.append( ef_array )
    res_arr = np.array( res_list )
    mu_eps_f = np.sum( res_arr, 0 ) / res_arr.shape[0]
    Km = ( 1. - ccb.V_f_tot ) * ccb.E_m
    Kf = ccb.V_f_tot * np.max( ccb.sorted_E_f )
    #Km_func = interp1d( a, ccb.Emtrx, bounds_error = False, fill_value = ccb.Emtrx[-1] )
    control = mu_eps_f * Kf + em_ef * Km
    #print control
    #plt.plot( ef_discr, control )
    #print res_x * Vf * 
    
    plt.plot( ef_discr, mu_eps_f, linewidth = '2' , color = 'r' )
    mu_eps_f_c = ( control[-1] - Km * em_ef ) / Kf
    plt.plot( ef_discr, mu_eps_f_c, linewidth = '3', color = 'g' )
    #print mu_eps_f[0], mu_eps_f_c[0], mu_eps_f[0] / mu_eps_f_c[0]   
                    #plt.plot( ccb._x_arr, np.maximum( ccb._epsf0_arr[i] - depsf * np.abs( ccb._x_arr ), ccb._epsm_arr ), 'k' )
                #else:
                #pass
                    #plt.plot( ccb._x_arr, np.maximum( ccb._epsf0_arr[i] - depsf * np.abs( ccb._x_arr ), ccb._epsm_arr ), 'b' )
    plt.legend( loc = 'best' )
    plt.xlim( -5, 5 )
   
    plt.show()
