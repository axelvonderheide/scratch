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
import numpy as np
from spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array
from types import FloatType
from reinforcement import Reinforcement, ContinuousFibers, ShortFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from scipy.optimize import fsolve, broyden2, root
import time as t
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

def H( x ):
        return x > 0

class CompositeCrackBridge( HasTraits ):

    reinforcement_lst = List( Instance( Reinforcement ) )
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float
    c_mask = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_c_mask( self ):
        return ( self.sorted_lf == np.infty )
        
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
        E_c = self.E_m * ( 1. - self.V_f_tot ) + E_fibers
        return E_c * ( 1. + 1e-15 )

    sorted_theta = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_theta( self ):
        print 'initiation'
        '''sorts the integral points by bond in descending order'''
        depsf_arr = np.array( [] )
        V_f_arr = np.array( [] )
        E_f_arr = np.array( [] )
        xi_arr = np.array( [] )
        stat_weights_arr = np.array( [] )
        nu_r_arr = np.array( [] )
        r_arr = np.array( [] )
        # ##
        lf_arr = np.array( [] )
        phi_arr = np.array( [] )
        # ##
        for reinf in self.reinforcement_lst:
            n_int = len( np.hstack( ( np.array( [] ), reinf.depsf_arr ) ) )
            depsf_arr = np.hstack( ( depsf_arr, reinf.depsf_arr ) )
            V_f_arr = np.hstack( ( V_f_arr, np.repeat( reinf.V_f, n_int ) ) )
            E_f_arr = np.hstack( ( E_f_arr, np.repeat( reinf.E_f, n_int ) ) )
            xi_arr = np.hstack( ( xi_arr, np.repeat( reinf.xi, n_int ) ) )
            stat_weights_arr = np.hstack( ( stat_weights_arr,
                                          np.repeat( reinf.stat_weights, n_int ) ) )
            nu_r_arr = np.hstack( ( nu_r_arr, reinf.nu_r ) )
            r_arr = np.hstack( ( r_arr, reinf.r_arr ) )
            if hasattr( reinf, 'phi_arr' ):
                phi_arr = np.hstack( ( phi_arr, reinf.phi_arr ) )
                lf_arr = np.hstack( ( lf_arr, np.repeat( reinf.l_f, n_int ) ) )
            else:
                phi_arr = np.hstack( ( phi_arr, np.zeros_like( reinf.depsf_arr ) ) )
                lf_arr = np.hstack( ( lf_arr, np.repeat( np.Infinity, n_int ) ) )
        argsort = np.argsort( depsf_arr )[::-1]
        # sorting the masks for the evaluation of F
        idxs = np.array( [] )
        for i, reinf in enumerate( self.reinforcement_lst ):
            idxs = np.hstack( ( idxs, i * np.ones_like( reinf.depsf_arr ) ) )
        masks = []
        for i, reinf in enumerate( self.reinforcement_lst ):
            masks.append( ( idxs == i )[argsort] )
        max_depsf = [np.max( reinf.depsf_arr ) for reinf in self.reinforcement_lst]
        masks = [masks[i] for i in np.argsort( max_depsf )[::-1]]
        return depsf_arr[argsort], V_f_arr[argsort], E_f_arr[argsort], \
                xi_arr[argsort], stat_weights_arr[argsort], \
                nu_r_arr[argsort], masks, r_arr[argsort], lf_arr[argsort], phi_arr[argsort]

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

    sorted_r = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_r( self ):
        return self.sorted_theta[7]
    
    sorted_lf = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_lf( self ):
        return self.sorted_theta[8]
    
    sorted_phi = Property( depends_on = 'reinforcement_lst+' )
    @cached_property
    def _get_sorted_phi( self ):
        return self.sorted_theta[9]

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
                methods.append( reinf.xi.weibull_fibers_cdf )
        return methods, masks

    def vect_xi_cdf( self, epsy, x_short, x_long ):
        Pf = np.zeros_like( self.sorted_depsf )
        methods, masks = self.sorted_xi_cdf
        for i, method in enumerate( methods ):
            if method.__name__ == 'weibull_fibers_cdf':
                Pf += method( epsy * masks[i], self.sorted_depsf,
                             x_short, x_long, self.sorted_r )
            else:
                Pf += method( epsy * masks[i] )
                # print 'test', method( epsy * masks[i] ), epsy * masks[i]
        return Pf
    
    # ##
    amin_it = Array
    def _amin_it_default( self ):
        return np.array( [0, 50 ] )

    def geo_amin_os( self, damage, Lmin, Lmax ):
            phi = self.sorted_phi
            Vf = self.sorted_V_f
            lf = self.sorted_lf
            depsmax = self.sorted_depsf[0]
            Ef = self.sorted_E_f
            a = np.linspace( 0, self.amin_it[-1] * 1.2 , self.discr_amin )
            Kf = Vf * self.sorted_nu_r * self.sorted_stats_weights * Ef
            a_shaped = a.reshape( len( a ), 1 )
            phi = phi.reshape( 1, len( phi ) )
            m_p = a_shaped * 2 / np.cos( phi ) / lf
            mask1 = m_p > 0
            m_p = m_p * mask1
            p = np.abs( H( 1 - m_p ) * ( 1 - m_p ) )
            ####
            muT = np.sum( self.sorted_depsf * ( 1 - damage ) * Kf * p , 1 )
            Kf_intact = np.sum( ( 1 - damage ) * Kf * ( 1 - p ) , 1 ) 
            Kf_broken = np.sum( Kf * damage )
            Emtrx = ( 1. - self.V_f_tot ) * self.E_m + Kf_broken + Kf_intact
            depsm = muT / Emtrx
            emr = np.hstack( ( 0, cumtrapz( depsm, a ) ) )
            umr = np.hstack( ( 0, cumtrapz( emr , a ) ) )
            em_func = interp1d( a, emr )
            um_func = interp1d( a, umr )
            emLl = em_func( Lmin )
            umLl = um_func( Lmin )
            condI = self.w - ( Lmin * ( emr - emLl + ( a - Lmin ) * depsmax ) + depsmax * Lmin ** 2. / 2. + emLl * Lmin - umLl + ( depsmax * a ** 2. / 2. + emr * a - umr ) )
            ip_f = interp1d( condI[::-1], a[::-1], bounds_error = False, fill_value = a[-1] )
            cut_at = np.nonzero( H( ip_f( 0 ) - a ) )
            a_new = np.hstack( ( a[cut_at], ip_f( 0 ) ) )
            amin = a_new
            self.amin_it = a_new
            # interp depsm and em
            ind = len( amin )
            interp_data = ( a_new[-1] - a[ind - 2] ) / ( a[ind - 1 ] - a[ind - 2] )
            diff = depsm[ind - 1 ] - depsm[ind - 2]
            depsm[ind - 1] = depsm[ind - 2] + interp_data * diff
            return amin, depsm[:len( amin ) ]

        
    def geo_amin_lmin( self, damage, depsf, umLmin, emLmin, Lmin, idx ):
            phi = self.sorted_phi
            Vf = self.sorted_V_f
            lf = self.sorted_lf
            depsmax = depsf
            Ef = self.sorted_E_f
            a = np.linspace( Lmin, Lmin * 1.5, self.discr_amin )
            Kf = Vf * self.sorted_nu_r * self.sorted_stats_weights * Ef
            a_shaped = a.reshape( len( a ), 1 )
            phi = phi.reshape( 1, len( phi ) )
            m_p = a_shaped * 2 / np.cos( phi ) / lf
            mask1 = m_p > 0
            m_p = m_p * mask1
            p = np.abs( H( 1 - m_p ) * ( 1 - m_p ) )
            ####
            muT = np.sum( self.sorted_depsf[idx:] * ( 1 - damage[idx:] ) * Kf[idx:] * p[:, idx:] , 1 )
            Kf_intact = np.sum( ( ( 1 - damage ) * Kf * ( 1 - p ) )[:, idx:] , 1 )
            Kf_broken = np.sum( Kf * damage )
            Emtrx = ( 1. - self.V_f_tot ) * self.E_m + Kf_broken + Kf_intact
            depsm = muT / Emtrx
            emr = np.hstack( ( emLmin, cumtrapz( depsm, a ) ) )
            umr = np.hstack( ( umLmin, cumtrapz( emr , a ) ) )
            condI = self.w - ( Lmin * ( emr - emLmin + ( a - Lmin ) * depsmax ) + depsmax * Lmin ** 2. / 2. + emLmin * Lmin - umLmin + ( depsmax * a ** 2. / 2. + emr * a - umr ) )
            ip_f = interp1d( condI[::-1], a[::-1], bounds_error = False, fill_value = a[-1] )
            cut_at = np.nonzero( H( ip_f( 0 ) - a ) )
            amin = np.hstack( ( a[cut_at], ip_f( 0 ) ) )
            return amin[-1]


    def geo_amin( self, damage, Lmax ):
            phi = self.sorted_phi
            Vf = self.sorted_V_f
            lf = self.sorted_lf
            depsmax = self.sorted_depsf[0]
            Ef = self.sorted_E_f
            a = np.linspace( 0, np.min( ( self.amin_it[-1] * 1.2, Lmax + 1e-6 ) ) , self.discr_amin )
            Kf = Vf * self.sorted_nu_r * self.sorted_stats_weights * Ef
            a_shaped = a.reshape( len( a ), 1 )
            phi = phi.reshape( 1, len( phi ) )
            m_p = a_shaped * 2 / np.cos( phi ) / lf
            mask1 = m_p > 0
            m_p = m_p * mask1
            p = np.abs( H( 1 - m_p ) * ( 1 - m_p ) )
            ####
            muT = np.sum( self.sorted_depsf * ( 1 - damage ) * Kf * p , 1 )
            Kf_intact = np.sum( ( 1 - damage ) * Kf * ( 1 - p ) , 1 ) 
            Kf_broken = np.sum( Kf * damage )
            Emtrx = ( 1. - self.V_f_tot ) * self.E_m + Kf_broken + Kf_intact
            depsm = muT / Emtrx
            em = np.hstack( ( 0, cumtrapz( depsm, a ) ) )
            um = np.hstack( ( 0, cumtrapz( em , a ) ) )
            condI = self.w / 2. - ( depsmax * a ** 2. / 2. + em * a - um )
            ip_f = interp1d( condI[::-1], a[::-1], bounds_error = False, fill_value = a[-1] )
            cut_at = np.nonzero( H( ip_f( 0 ) - a ) )
            a_new = np.hstack( ( a[cut_at], ip_f( 0 ) ) )
            amin = a_new
            self.amin_it = a_new
            # interp depsm and em
            ind = len( amin )
            interp_data = ( a_new[-1] - a[ind - 2] ) / ( a[ind - 1 ] - a[ind - 2] )
            diff = depsm[ind - 1 ] - depsm[ind - 2]
            depsm[ind - 1] = depsm[ind - 2] + interp_data * diff
            return amin, depsm[:len( amin ) ]
    
    def dem_depsf_vect( self , damage, a_unshaped ):
        '''evaluates the deps_m given deps_f
        at that point and the damage array'''
        # print a_unshaped
        # ##
        phi = self.sorted_phi
        Vf = self.sorted_V_f
        lf = self.sorted_lf
        depsmax = self.sorted_depsf[0]
        Ef = self.sorted_E_f
        # ##
        Kf = self.sorted_V_f * self.sorted_nu_r * \
            self.sorted_stats_weights * self.sorted_E_f
        
        Kf_broken = np.sum( Kf * damage )
        
        # transforming Kf and damage into matrix array with length of a rows
        # getting lf of fibers
        lf = self.sorted_lf
        # reshaping a and phi to get an evaluation in every a
        
        a = a_unshaped.reshape( len( a_unshaped ), 1 )
        phi = self.sorted_phi.reshape( 1, len( self.sorted_phi ) )
        # geometrical condition evaluated in every a
        m_p = a * 2 / np.cos( phi ) / lf  
        # continious fibers are flagged with a negative value in m_p- The following mask is to filter the array.
        mask1 = m_p > 0
        m_p = m_p * mask1
        # m_p is a matrix array giving for every phi in a the percentage of fibers that are geometrically dropped out. 
        # p is the percentage of active fibers. Heaviside cuts all values exceeding 100% drop out.
        p = np.abs( H( 1 - m_p ) * ( 1 - m_p ) )
        # creating matrix that implements the mechanical condition. 
        # In every  step in 'a' there is one integral point of bond less active.
        # So it is a diagonal matrix with ones on the right side and zeroes on the other.
        stepm = np.eye( len( a ), len( phi[0] ) )
        eym = np.cumsum( stepm , -1 )
        # summing up all the matrix arrays in the axis of 'a' for muT and Kf
        # print p.shape, depsf.shape, damage.shape, Kf.shape, eym.shape
        eym_inv = eym == 0
        mask_2conditions = ( ( 1 - p ) + eym_inv ) < 1
        mu_T = np.sum( self.sorted_depsf * ( 1 - damage ) * Kf * p * eym, 1 )
        Kf_intact_bonded = np.sum( ( 1 - damage ) * Kf * ( ( ( 1 - p ) + eym_inv ) * mask_2conditions + ( 1 - mask_2conditions ) ) , 1 )  # [::-1]
        Kf_add = Kf_intact_bonded + Kf_broken
        Km = ( 1. - self.V_f_tot ) * self.E_m
        E_mtrx = Km + Kf_add
        return mu_T / E_mtrx

    def F( self, dems, amin ):
        F = np.zeros_like( self.sorted_depsf )
        for i, mask in enumerate( self.sorted_masks ):
            depsfi = self.sorted_depsf[mask]
            demsi = dems[mask]
            fi = 1. / ( depsfi + demsi )
            F[mask] = np.hstack( ( np.array( [0.0] ), cumtrapz( fi, -depsfi ) ) )
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
        # plt.plot(-self.sorted_depsf, F)
        # plt.show()
        return F
    
    if 0 == 0:
        def profile( self, iter_damage, Lmin, Lmax ):
            # matrix strain derivative with resp. to z as a function of T
            if np.any( iter_damage < 0.0 ) or np.any( iter_damage > 1.0 ):
                return np.ones_like( iter_damage ) * 0.5, np.ones_like( self.sorted_depsf ), np.ones_like( self.sorted_depsf )
            # old code
            # ind_a = np.argmin( np.abs( self._x_arr ) )
            # imag_z = len( self.sorted_V_f ) - len( self._a_long )
            # if imag_z > 0:
            #    print 'Randbedingung greift'
            #    a_unshaped = np.hstack( ( self._a_long, self._a_long[-1] + np.cumsum( np.repeat( self._a_long[-1] - self._a_long[-2], imag_z ) ) ) )
            # else:
            #    a_unshaped = self._a_long
            #    print len(self._a_long)
            #####
            dems = self.dem_depsf_vect( iter_damage, self._a_long )
            # print 'dems', dems
            # initial matrix strain derivative
            init_dem = dems[0]
            # debonded length of fibers with Tmax
            if  self._x_arr[0] == 1e-10:
                amin = ( self.w / ( np.abs( init_dem ) + np.abs( self.sorted_depsf[0] ) ) ) ** 0.5
                self.amin_it = np.array( [0, amin] )
            a_geo, depsm_geo = self.geo_amin( iter_damage, Lmax )
            amin = a_geo[-1]
            # integrated f(depsf) - see article
            F = self.F( dems, amin )
            # a(T) for double sided pullout
            a1 = np.exp( F / 2. + np.log( amin ) )
            if Lmin < a1[0] and Lmax < a1[0]:
                # all fibers debonded up to Lmin and Lmax
                print 'Lmin=Lmax<amin'
                min_a = np.sum( ( Lmin - a_geo ) >= 0 )
                # min_a = a_d[min_a_i]
                max_a = np.sum( ( Lmax - a_geo ) >= 0 )
                # max_a = a_d[max_a_i]
                l_a = np.hstack( ( -Lmin, -a_geo[:min_a][::-1] ) )
                r_a = np.hstack( ( a_geo[:max_a], Lmax ) )
                r_depsm = depsm_geo[:max_a + 1] 
                l_depsm = depsm_geo[:min_a + 1]
                r_em = cumtrapz( r_depsm, r_a )
                l_em = cumtrapz( l_depsm, -l_a[::-1] )[::-1]
                a = np.hstack( ( l_a[:-1], 0.0, r_a[1:] ) )
                em = np.hstack( ( l_em, 0.0, r_em ) )
                umL = np.trapz( np.hstack( ( 0, l_em [::-1] ) ), -l_a[::-1] )
                umR = np.trapz( np.hstack( ( 0, r_em ) ), r_a )
                vxL = ( self.w - Lmin * l_em[0] - Lmax * r_em[-1] + umL + umR - self.sorted_depsf * ( Lmin ** 2 + Lmax ** 2 ) / 2 - \
                        Lmax * ( l_em[0] + ( Lmin - Lmax ) * self.sorted_depsf - em[-1] ) ) / ( Lmin + Lmax )
                epsf0 = Lmin * self.sorted_depsf + vxL + l_em[0]
            elif Lmin < a1[0] and Lmax >= a1[0]:
                # all fibers debonded up to Lmax but not up to Lmin
                
                min_a = np.sum( ( Lmin - a_geo ) > 0 )
                l_a = np.hstack( ( -Lmin, -a_geo[:min_a][::-1] ) )
                l_depsm = depsm_geo[:min_a + 1]
                l_em = cumtrapz( l_depsm, -l_a[::-1] )[::-1]
                umL = np.trapz( np.hstack( ( 0, l_em [::-1] ) ), -l_a[::-1] )
                # amin numerisch neu bestimmen fuer SF (geo_amin one sided)
                almin_arr, depsm_lmin = self.geo_amin_os( iter_damage, Lmin, Lmax )
                a_geo = almin_arr
                depsm_geo = depsm_lmin
                amin = almin_arr[-1]
                C = np.log( amin ** 2 + 2 * Lmin * amin - Lmin ** 2 )
                a2 = np.sqrt( 2 * Lmin ** 2 + np.exp( ( F + C ) ) ) - Lmin
                if Lmax <= a2[-1]:
                    print 'Lmin <amin und Lmax <a2[-1]'
                    idx = np.sum( a2 < Lmax ) - 1
                    a = np.hstack( ( l_a[:-1], 0.0, a_geo[1:-1], a2[:idx + 1], Lmax ) )
                    em22 = cumtrapz( depsm_geo[:-1], a_geo[:-1] )
                    em2 = em22[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a2 ) ) ) * dems ) 
                    em = np.hstack( ( l_em, 0.0, em22, em2[:idx + 1], em2[idx] + ( Lmax - a2[idx] ) * dems[idx] ) )
                    um = np.trapz( em, a )
                    epsf01 = em2[:idx + 1] + a2[:idx + 1] * self.sorted_depsf[:idx + 1]
                    # if a2[0] > Lmax:
                    #    print 'dangerous'
                    #    idx = 0
                    vx02L = ( self.w - Lmin * l_em[0] - Lmax * em[-1] + um - self.sorted_depsf[idx + 1:] * ( Lmin ** 2 + Lmax ** 2 ) / 2 - \
                           Lmax * ( l_em[0] + ( Lmin - Lmax ) * self.sorted_depsf[idx + 1:] - em[-1] ) ) / ( Lmin + Lmax )
                    epsf02 = vx02L + self.sorted_depsf[idx + 1:] * Lmin + l_em[0]
                    epsf0 = np.hstack( ( epsf01, epsf02 ) )
                else:
                    print 'Lmin <amin und Lmax >a2[-1]'
                    a = np.hstack( ( l_a[:-1], 0.0, a_geo[1:-1], a2, Lmax ) )
                    em22 = cumtrapz( depsm_geo[:-1], a_geo[:-1] )
                    em2 = em22[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a2 ) ) ) * dems ) 
                    em = np.hstack( ( l_em, 0.0, em22, em2, em2[-1] ) )
                    epsf0 = em2 + self.sorted_depsf * a2
            
            elif a1[0] < Lmin and a1[-1] > Lmin:
                # boundary condition position
                print 'boundary'
                idx1 = np.sum( a1 <= Lmin )
                # a(T) for one sided pullout
                # amin for one sided PO
                depsfLmin = self.sorted_depsf[idx1]
                a_short = np.hstack( ( a1[:idx1], Lmin ) )
                # ##
                em11 = cumtrapz( depsm_geo[:-1], a_geo[:-1] )
                em1 = np.hstack( ( em11, em11[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a1 ) ) ) * dems ) ) )
                em = np.hstack( ( em1[-1], em1[::-1], 0.0, em1, em1[-1] ) )
                # ##
                em22 = cumtrapz( depsm_geo[:-1], a_geo[:-1] )
                em_short = em22[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a_short ) ) ) * dems[:idx1 + 1] )
                em_short_full = np.hstack( ( 0, em22, em_short ) )
                emLmin = em_short[-1]
                umLmin = np.trapz( em_short_full, np.hstack( ( a_geo[:-1], a_short ) ) )
                # amin numerisch
                amin = self.geo_amin_lmin( iter_damage, depsfLmin, umLmin, emLmin, Lmin, idx1 )
                
                # amin = -Lmin + np.sqrt( 4 * Lmin ** 2 * p ** 2 - 4 * p * emLmin * Lmin + 4 * p * umLmin - 2 * p * Lmin ** 2 * depsfLmin + 2 * p * self.w ) / p
                C = np.log( amin ** 2 + 2 * amin * Lmin - Lmin ** 2 )
                a2 = ( np.sqrt( 2 * Lmin ** 2 + np.exp( F + C - F[idx1] ) ) - Lmin )[idx1:]
                # matrix strain profiles - shorter side
                a_short = np.hstack( ( -a_short[::-1], 0.0 ) )
                dems_short = dems[:idx1]
                em_short = np.hstack( ( em_short[::-1] , 0 ) )
                if a2[-1] > Lmax:
                    print 'a2[-1]>Lmax'
                    idx2 = np.sum( a2 <= Lmax )
                    # matrix strain profiles - longer side
                    a_long = np.hstack( ( a1[:idx1], a2[:idx2] ) )
                    em_long = em22[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a_long ) ) ) * dems[:idx1 + idx2] )
                    a = np.hstack( ( a_short[:-1], -a_geo[::-1][1:-1], 0, a_geo[1:-1], a_long, Lmax ) )
                    em = np.hstack( ( em_short[:-1], em22[::-1], 0, em22, em_long, em_long[-1] + ( Lmax - a_long[-1] ) * dems[idx1 + idx2] ) )
                    um = np.trapz( em, a )
                    epsf01 = em_long + a_long * self.sorted_depsf[:idx1 + idx2]
                    vx02 = ( self.w - Lmin * em[0] - Lmax * em[-1] + um - self.sorted_depsf[idx1 + idx2:] * ( Lmin ** 2 + Lmax ** 2 ) / 2 - \
                           Lmax * ( em[0] + ( Lmin - Lmax ) * self.sorted_depsf[idx1 + idx2:] - em[-1] ) ) / ( Lmin + Lmax )
                    epsf02 = vx02 + self.sorted_depsf[idx1 + idx2:] * Lmin + em[0]
                    epsf0 = np.hstack( ( epsf01, epsf02 ) )
                else:
                    print 'a2[-1]<Lmax'
                    a_long = np.hstack( ( 0.0, a1[:idx1], a2, Lmax ) )
                    # print a1[:idx1], a2
                    a = np.hstack( ( a_short[:-1], -a_geo[::-1][1:-1], 0, a_geo[1:-1], a_long[1:] ) )
                    dems_long = dems
                    em_long = em22[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a_long[1:-1] ) ) ) * dems_long )
                    em_long = np.hstack( ( em_long, em_long[-1] ) )
                    em = np.hstack( ( em_short[:-1], em22[::-1], 0, em22, em_long ) )
                    epsf0 = em_long[:-1] + self.sorted_depsf * a_long[1:-1]
                    a_long = a_long[1:-1]
                    
            elif a1[-1] <= Lmin:
                # double sided pullout
                print 'double sided'
                a = np.hstack( ( -Lmin, -a1[::-1], -a_geo[1:-1][::-1], 0.0, a_geo[1:-1], a1, Lmin ) )
                em11 = cumtrapz( depsm_geo[:-1], a_geo[:-1] )
                em1 = np.hstack( ( em11, em11[-1] + np.cumsum( np.diff( np.hstack( ( a_geo[-2], a1 ) ) ) * dems ) ) )
                em = np.hstack( ( em1[-1], em1[::-1], 0.0, em1, em1[-1] ) )
                epsf0 = em1[len( em11 ):] + self.sorted_depsf * a1
            # print a[-5:],em[-5:]
            self._x_arr = a
            self._epsm_arr = em
            self._epsf0_arr = epsf0
            a_short = -a[a < 0.0][1:][::-1][len( a_geo[1:-1] ):]
            if len( a_short ) < len( self.sorted_depsf ):
                a_short = np.hstack( ( a_short, Lmin * np.ones( len( self.sorted_depsf ) - len( a_short ) ) ) )
            a_long = a[a > 0.0][:-1][len( a_geo[1:-1] ):]
            if len( a_long ) < len( self.sorted_depsf ):
                a_long = np.hstack( ( a_long, Lmax * np.ones( len( self.sorted_depsf ) - len( a_long ) ) ) )
            self._a_long = a_long
            return epsf0, a_short, a_long
    
    else: 
        '''
        def profile( self, iter_damage, Lmin, Lmax ):
            if np.any( iter_damage < 0.0 ) or np.any( iter_damage > 1.0 ):
                return np.ones_like( iter_damage ) * 0.5, np.ones_like( self.sorted_depsf ), np.ones_like( self.sorted_depsf )
            # matrix strain derivative with resp. to z as a function of T
            dems = self.dem_depsf_vect( iter_damage, np.array( [] ) )
            # initial matrix strain derivative
            init_dem = dems[0]
            # debonded length of fibers with Tmax
            amin = ( self.w / ( np.abs( init_dem ) + np.abs( self.sorted_depsf[0] ) ) ) ** 0.5
            # integrated f(depsf) - see article
            F = self.F( dems, amin )
            # a(T) for double sided pullout
            a1 = np.exp( F / 2. + np.log( amin ) )# * mask_po + Le * inv_mask_po
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
                    a = np.hstack( ( a_short, a_long[1:] ) )
                    dems_long = dems
                    em_long = np.hstack( ( np.cumsum( np.diff( a_long[:-1] ) * dems_long ) ) )
                    em_long = np.hstack( ( em_long, em_long[-1] ) )
                    em = np.hstack( ( em_short, em_long ) )
                    epsf0 = em_long[:-1] + self.sorted_depsf * a_long[1:-1]
            elif a1[-1] <= Lmin:
                #double sided pullout
                
                a = np.hstack( ( -Lmin, -a1[::-1], 0.0, a1, Lmax ) )
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
            self._a_long = a_long
            
            return epsf0, a_short, a_long
    '''

    def damage_residuum( self, iter_damage ):
        Lmin = min( self.Ll, self.Lr )
        Lmax = max( self.Ll, self.Lr )
        epsf0, x_short, x_long = self.profile( iter_damage, Lmin, Lmax )
        residuum = self.vect_xi_cdf( epsf0, x_short = x_short, x_long = x_long ) - iter_damage
        return residuum
    _x_arr = Array
    def __x_arr_default( self ):
        return np.hstack( ( np.repeat( 1e-10, len( self.sorted_depsf ) + 1 ), 0, np.repeat( 1e-10, len( self.sorted_depsf ) + 1 ) ) )
    
    _a_long = Array
    def __a_long_default( self ):
        return np.hstack( ( np.repeat( 1e-10, len( self.sorted_depsf ) ) ) )

    _epsm_arr = Array
    def __epsm_arr_default( self ):
        return np.repeat( 1e-10, len( self.sorted_depsf ) )

    _epsf0_arr = Array
    def __epsf0_arr_default( self ):
        return np.repeat( 1e-10, len( self.sorted_depsf ) )
    
    _first_guess = Array
    def __first_guess_default( self ):
        return  np.ones_like( self.sorted_depsf ) * 1e-10
    
    def set_sfarr_to_defaults( self ):
        self._a_long = self.__a_long_default()
        self._epsf0_arr = self.__epsf0_arr_default()
        self._x_arr = self.__x_arr_default()
        self._epsm_arr = self.__epsm_arr_default()

        
    damage = Property( depends_on = 'w, Ll, Lr, reinforcement+' )
    @cached_property
    def _get_damage( self ):
        if self.w == 0.:
            damage = np.zeros_like( self.sorted_depsf )
        else:
            ff = t.clock()
            try:
                self.set_sfarr_to_defaults()
                damage = root( self.damage_residuum, np.ones_like( self.sorted_depsf ) * 1e-10, method = 'excitingmixing' )
                if np.any( damage.x < 0.0 ) or np.any( damage.x > 1.0 ):
                    raise ValueError
                damage = damage.x
            except:
                print 'fast opt method does not converge: switched to a slower, robust method for this step'
                
                damage = root( self.damage_residuum, self._first_guess, method = 'krylov' )
                damage = damage.x
                self._first_guess = damage
            # print 'damage =', np.sum(damage) / len(damage), 'iteration time =', t.clock() - ff, 'sec'
        return damage

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    reinf1 = ContinuousFibers( r = 0.00345,
                          tau = RV( 'uniform', loc = 0.01, scale = 0.1 ),
                          V_f = 0.0111,
                          E_f = 180e3,
                          xi = WeibullFibers( shape = 4., sV0 = 0.0025 ),
                          n_int = 300,
                          label = 'carbon' )

    
    reinfSF = ShortFibers( r = 0.1,
                          tau = 7.,
                          lf = .3,
                          snub = 2.,
                          phi = RV( 'sin2x', loc = 0., scale = 1. ),
                          V_f = 0.01,
                          E_f = 200e3,
                          xi = 1.,  # WeibullFibers( shape = 1000., scale = 1000 ),
                          n_int = 300,
                          label = 'Short Fibers' )

    ccb = CompositeCrackBridge( E_m = 25e3,
                                 reinforcement_lst = [reinfSF, reinf1],
                                 Ll = 1.,
                                 Lr = 2.,
                                 w = 0.002,
                                 discr_amin = 100 )

#    reinf = ContinuousFibers(r=0.01,
#                          tau=RV('uniform', loc=0.01, scale=.5),
#                          V_f=0.05,
#                          E_f=200e3,
#                          xi=WeibullFibers(shape=5., sV0=0.00618983207723),
#                          n_int=10,
#                          label='carbon')
#
#    ccb = CompositeCrackBridge(E_m=25e3,
#                                 reinforcement_lst=[reinf],
#                                 Ll=4.,
#                                 Lr=87.,
#                                 w=0.004)

    ccb.damage
    plt.plot( ccb._x_arr, ccb._epsm_arr, lw = 2, color = 'red', label = 'analytical' )  # , ls = 'dashed' )
    plt.plot( np.zeros_like( ccb._epsf0_arr ), ccb._epsf0_arr, 'ro' )
    for i, depsf in enumerate( ccb.sorted_depsf ):
        # print ccb.sorted_lf[i]
        if ccb.sorted_lf[i] == np.infty:
           plt.plot( ccb._x_arr, np.maximum( ccb._epsf0_arr[i] - depsf * np.abs( ccb._x_arr ), ccb._epsm_arr ) )
    plt.legend( loc = 'best' )
    plt.xlim( -20, 20 )
    plt.show()
