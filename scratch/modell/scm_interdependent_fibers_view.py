'''
Created on Jul 26, 2012

@author: rostar
'''

from etsproxy.traits.api import \
    Instance, Array, List, cached_property, Property, Int
from etsproxy.traits.ui.api import ModelView
from spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
import copy
from math import pi
from scm_interdependent_fibers_model import SCM
from reinforcement import Reinforcement, ContinuousFibers, ShortFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from hom_CB_elastic_mtrx import CompositeCrackBridge
from hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
import pickle
import os
from matplotlib import pyplot as plt


class SCMView( ModelView ):

    model = Instance( SCM )
    def crack_widths( self, sigma_c ):
        # find the index of the nearest value in the load range
        idx = np.abs( self.model.load_sigma_c_arr - sigma_c ).argmin()
        # evaluate the relative strain e_rel between fibers
        # and matrix for the given load
        e_rel = self.mu_epsf_x[idx, :] - self.eps_m_x[idx, :]
        # pick the cracks that emerged at the given load
        cb_load = self.model.cb_list( sigma_c )
        if cb_load[0] is not None:
            # find the symmetry points between cracks as
            # the 0 element of their x range
            idxs = []
            for cb in cb_load:
                idxs.append( np.where( cb.position + 
                            cb.x[0] == self.model.x_arr )[0] )
            # add the index of the last point
            idxs.append( self.model.nx - 1 )
            # list of crack widths to be filled in a loop with integrated e_rel
            crack_widths = [np.trapz( e_rel[idx:idxs[i + 1]],
                            self.model.x_arr[idx:idxs[i + 1]] )
                            for i, idx in enumerate( idxs[:-1] )]
            return np.array( crack_widths, ndmin = 1 )
        else:
            return np.array( 0.0, ndmin = 1 )

    eval_w = Property( List, depends_on = 'model' )
    @cached_property
    def _get_eval_w( self ):
        return [self.crack_widths( load ) for load in self.model.load_sigma_c_arr]

    w_mean = Property( Array, depends_on = 'model' )
    @cached_property
    def _get_w_mean( self ):
        return np.array( [np.mean( w ) for w in self.eval_w] )

    w_median = Property( Array, depends_on = 'model' )
    @cached_property
    def _get_w_median( self ):
        return np.array( [np.median( w ) for w in self.eval_w] )

    w_stdev = Property( Array, depends_on = 'model' )
    @cached_property
    def _get_w_stdev( self ):
        return np.array( [np.std( w ) for w in self.eval_w] )

    w_max = Property( Array, depends_on = 'model' )
    @cached_property
    def _get_w_max( self ):
        return np.array( [np.max( w ) for w in self.eval_w] )

    x_area = Property( depends_on = 'model.' )
    @cached_property
    def _get_x_area( self ):
        return  np.ones_like( self.model.load_sigma_c_arr )[:, np.newaxis] \
            * self.model.x_arr[np.newaxis, :]

    sigma_m_x = Property( depends_on = 'model.' )
    @cached_property
    def _get_sigma_m_x( self ):
        sigma_m_x = np.zeros_like( self.model.load_sigma_c_arr[:, np.newaxis]
                                  * self.model.x_arr[np.newaxis, :] )
        for i, q in enumerate( self.model.load_sigma_c_arr ):
            sigma_m_x[i, :] = self.model.sigma_m( q )
        return sigma_m_x

    eps_m_x = Property( Array, depends_on = 'model.' )
    @cached_property
    def _get_eps_m_x( self ):
        return self.sigma_m_x / self.model.CB_model.E_m

    mu_epsf_x = Property( depends_on = 'model.' )
    @cached_property
    def _get_mu_epsf_x( self ):
        mu_epsf_x = np.zeros_like( self.model.load_sigma_c_arr[:, np.newaxis]
                                  * self.model.x_arr[np.newaxis, :] )
        for i, q in enumerate( self.model.load_sigma_c_arr ):
            mu_epsf_x[i, :] = self.model.epsf_x( q )
        return mu_epsf_x

    eps_sigma = Property( depends_on = 'model.' )
    @cached_property
    def _get_eps_sigma( self ):
        eps = np.trapz( self.mu_epsf_x, self.x_area, axis = 1 ) / self.model.length
        eps = eps[np.isnan( eps ) == False]
        if len( eps ) != len( self.model.load_sigma_c_arr ):
            eps = list( eps ) + [list( eps )[-1]]
            sigma = copy.copy( self.model.load_sigma_c_arr[:len( eps )] )
            sigma[-1] = 0.0
            return eps, sigma
        else:
            return eps, self.model.load_sigma_c_arr
        
    def save_to_file( self ):
        # 1 eps sigma
        eps, sigma = self.eps_sigma
        combined_es = open( 'combined_es.pkl', 'wb' )
        pickle.dump( [eps, sigma], combined_es, -1 )
        combined_es.close()
        # 2 w sigma
        sigma = self.model.load_sigma_c_arr
        combined_sw = open( 'combined_sw.pkl', 'wb' )
        pickle.dump( [sigma, self.w_mean, self.w_max, self.w_stdev], combined_sw, -1 )
        combined_sw.close()
        # 3 Histogram
        hist_list = []
        for load in sigma:
            hist_list.append( self.crack_widths( load ) )
        combined_hist = open( 'combined_hist.pkl', 'wb' )
        pickle.dump( np.array( hist_list ), combined_hist, -1 )
        combined_hist.close()
        return 0

if __name__ == '__main__':
    length = 550.
    nx = 2000
    

    def get_rf( VfTex, VfSf, Em, Ec ):
        print Ec
        c_scale = ( 1.4 + VfSf * 100.*1.6 ) * Em / Ec
        c_shape = 12  # + ( VfSf * 2 ) * 100
        return RandomField( seed = True,
                               lacor = .005,
                                xgrid = np.linspace( 0., length, 550 ),
                                nsim = 1,
                                loc = .0,
                                shape = c_shape,
                                scale = c_scale,
                                non_negative_check = True,
                                distribution = 'Weibull'
                               )

    

    
    # tex_carbonfibers_real = ContinuousFibers( r = 0.0035,
    #                      tau = RV( 'weibull_min', loc = 0.006, shape = .26, scale = .03 ),
    #                      V_f = 0.011,
    #                      E_f = 240e3,
    #                      xi = WeibullFibers( shape = 5.0, sV0 = 0.0026 ),
    #                      label = 'TEXCarbon_real' )
    
    tex_glas = ContinuousFibers( r = 0.0095,
                          tau = RV( 'weibull_min', loc = .0014, shape = 0.28, scale = 0.008 ),
                          V_f = 0.011,
                          E_f = 72e3,
                          xi = WeibullFibers( shape = 4.5, sV0 = 1.856e-3 ),
                          label = 'TEXGlas' )
    
    tex_carbonfibers = ContinuousFibers( r = 0.0035,
                          tau = RV( 'weibull_min', loc = 0.015, shape = 5., scale = .004 ) ,  # RV( 'weibull_min', loc = 0.006, shape = .23, scale = .03 ),  # RV( 'uniform', loc = 0.5, scale = 1.5 ),  # RV( 'weibull_min', loc = 0.006, shape = .23, scale = .03 ),
                          V_f = 0.011,
                          E_f = 180e3,
                          n_int = 300,
                          xi = WeibullFibers( shape = 5.00001, sV0 = 0.002501 ),
                          label = 'carbon' )

   
     
    sf_steel = ShortFibers( r = 0.3 ,
                          tau = 1.76,
                          lf = 17.,
                          snub = .03,
                          phi = RV( 'sin2x', loc = 0., scale = 1. ),  # RV( 'uniform', loc = 0., scale = 1e-12 ),
                          V_f = 0.01,
                          E_f = 200e3,
                          xi = np.infty,  # WeibullFibers( shape = 1000., scale = 1000 ),
                          label = 'SFSteel' )
    
    sf_polyethylen = ShortFibers( r = 3.8e-3 ,
                          tau = .102,
                          lf = 12.7 ,
                          snub = .7,
                          phi = RV( 'sin2x', loc = 0., scale = 1. ),  # RV( 'uniform', loc = 0., scale = 1e-12 ),
                          V_f = 0.01,
                          E_f = 120e3,
                          xi = np.infty,  # WeibullFibers( shape = 1000., scale = 1000 ),
                          label = 'SFPolyethylen' )
    
    sf_glas = ShortFibers( r = 9.5 * 1e-3,
                          tau = 1.5,
                          lf = 10. ,
                          snub = 0.7,
                          phi = RV( 'sin2x', loc = 0., scale = 1. ),
                          V_f = 0.01,
                          E_f = 72e3,
                          xi = np.infty,
                          label = 'SFGlas' )
    
    def open_CB( r1, r2 ):
        return CompositeCrackBridge( E_m = 25e3,
                                 reinforcement_lst = [r1, r2],
                                 )
    
    def open_ini( CB_model, rf ):
        scm = SCM( length = length,
                  nx = nx,
                  n_w_interp = 30,
                  n_BC_interp = 5,
                  n_x_interp = 300,
                  piees = False,
                  piers = False,
                  random_field = rf,
                  CB_model = CB_model,
                  load_sigma_c_arr = np.linspace( 0.01, 15., 100 ),
                  d_int = 2
                  )
        scm_view = SCMView( model = scm )
        return scm_view
    
    V_f_list = [ 0.01]  # , 0.01, .015, 0.02, 0.025, 0.03]
    reinf_tex_list = [tex_carbonfibers]
    reinf_sf_list = [sf_glas]
    if ( os.path.isdir( 'scm_data' ) == False ):os.mkdir( 'scm_data' )
    os.chdir( 'scm_data' )
    foldername = 'all' 
    if ( os.path.isdir( foldername ) == False ):os.mkdir( foldername ) 
    os.chdir( foldername )
    distr3D = RV( 'sin2x', loc = 0., scale = 1. )
    distr2D = RV( 'uniform', loc = 0., scale = pi / 2. )
    distr_ls = [distr3D]
    for distr in distr_ls: 
        if ( os.path.isdir( distr.type ) == False ): os.mkdir( distr.type )
        os.chdir( distr.type )
        for reinf_tex in reinf_tex_list:
            if ( os.path.isdir( reinf_tex.label ) == False ):os.mkdir( reinf_tex.label )
            os.chdir( reinf_tex.label )
            for reinf_sf in reinf_sf_list:
                if ( os.path.isdir( reinf_sf.label ) == False ):os.mkdir( reinf_sf.label )
                os.chdir( reinf_sf.label )
                reinfSF = reinf_sf
                reinfTEX = reinf_tex
                for i, V_fi in enumerate( V_f_list ):
                    reinfSF.V_f = V_fi
                    Ec = 25e3 * ( 1 - V_fi - reinf_tex.V_f ) + V_fi * reinf_sf.E_f + reinf_tex.V_f * reinf_tex.E_f
                    reinfSF.phi = distr
                    CB = open_CB( reinfTEX, reinfSF )
                    rfVf = get_rf( reinf_tex.V_f, V_fi, 25e3, Ec )
                    scm_view = open_ini( CB, rfVf )
                    if ( os.path.isdir( np.str( V_fi ) ) == False ):os.mkdir( np.str( V_fi ) )
                    os.chdir( np.str( V_fi ) )
                    if ( os.path.isdir( 'InterpolatorData' ) == False ):os.mkdir( 'InterpolatorData' )

                    scm_view.model.evaluate()
                    scm_view.save_to_file()
                    os.chdir( os.pardir )
                    del scm_view
                    del CB
                    del rfVf
                    print 'global status:', i + 1, 'of', len( V_f_list ), 'done' 
                os.chdir( os.pardir )
            os.chdir( os.pardir ) 
        os.chdir( os.pardir ) 
        
    
        
    def plot():
        # scm_view.save_to_file()
        eps, sigma = scm_view.eps_sigma
        plt.figure()
        plt.plot( eps, sigma, color = 'black', lw = 2, label = 'model' )
        plt.legend( loc = 'best' )
        plt.xlabel( 'composite strain [-]' )
        plt.ylabel( 'composite stress [MPa]' )
        plt.figure()
        plt.hist( scm_view.crack_widths( 4. ), bins = 20, label = 'load = 2 MPa' )
        plt.hist( scm_view.crack_widths( 5. ), bins = 20, label = 'load = 15 MPa' )
        plt.hist( scm_view.crack_widths( 7.5 ), bins = 20, label = 'load = 10 MPa' )
        plt.legend( loc = 'best' )
        plt.figure()
        plt.plot( scm_view.model.load_sigma_c_arr, scm_view.w_mean,
                 color = 'green', lw = 2, label = 'mean crack width' )
        plt.plot( scm_view.model.load_sigma_c_arr, scm_view.w_median,
                 color = 'blue', lw = 2, label = 'median crack width' )
        plt.plot( scm_view.model.load_sigma_c_arr, scm_view.w_mean + scm_view.w_stdev,
                 color = 'black', label = 'stdev' )
        plt.plot( scm_view.model.load_sigma_c_arr, scm_view.w_mean - scm_view.w_stdev,
                 color = 'black' )
        plt.plot( scm_view.model.load_sigma_c_arr, scm_view.w_max,
                 ls = 'dashed', color = 'red', label = 'max crack width' )
        plt.legend( loc = 'best' )
        plt.show()
    # plot()
