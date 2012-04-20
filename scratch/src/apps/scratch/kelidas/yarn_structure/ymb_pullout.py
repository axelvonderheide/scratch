#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 14, 2010 by: kelidas

from enthought.traits.api import \
    HasTraits, Int, Array, Str, implements, Range, Property, cached_property, \
     Float, Instance, Any, Interface

from math import pi, e

from numpy import \
    sign, linspace, array, cos, sqrt, argmax, hstack, max, zeros_like, argwhere

from matplotlib import pyplot as plt

from stats.spirrid.i_rf import \
    IRF

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from stats.spirrid.rf import \
    RF

from ymb_data import YMBData

from stats.spirrid import SPIRRID, RV

class RV_add( HasTraits ):
    '''Class representing the definition and discretization of a random variable.
    '''
    name = Str

    pd = Instance( Distrib )

    n_int = Int( 30 )

    # index within the randomization
    idx = Int( 0 )

    theta_arr = Property( Array( 'float_' ), depends_on = 'distr' )
    @cached_property
    def _get_theta_arr( self ):
        return linspace( 0, self.pd.par_b, self.n_int )

    pdf_arr = Property( Array( 'float_' ), depends_on = 'distr' )
    @cached_property
    def _get_pdf_arr( self ):
        return self.pd.pdf( self.theta_arr )

    pdf_theta_arr = Property( Array( 'float_' ), depends_on = 'distr' )
    @cached_property
    def _get_pdf_theta_arr( self ):
        pdf_theta_arr = self.pdf_arr * ( self.theta_arr[1] - self.theta_arr[0] )
        pdf_theta_arr[ 0] *= .5
        pdf_theta_arr[-1] *= .5
        return pdf_theta_arr

from quaducom.crackbridge.yarn_symmetrical import DoublePulloutSym

class YMBPullOut( HasTraits ):
    '''Idealization of the double sided pullout using the SPIRRID 
    statistical integration tool.
    '''
    data = Instance( YMBData )

    rf = Instance( DoublePulloutSym )
    def _rf_default( self ):
        return DoublePulloutSym( tau_fr = 2.5, l = 0.01, d = 26e-3, E_mod = 72.0e3,
                                 theta = 0.01, xi = 0.014, phi = 1. )


    figure = Instance( Figure )
    def _figure_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.08, 0.13, 0.85, 0.74] )
        return figure

    pdf_theta = Property( Instance( YMBDistrib ) )
    @cached_property
    def _get_pdf_theta( self ):
        return YMBDistrib( varname = 'slack', data = self.data )

    pdf_l = Property( Instance( YMBDistrib ) )
    @cached_property
    def _get_pdf_l( self ):
        return YMBDistrib( varname = 'cut_free_length', data = self.data )

    pdf_phi = Property( Instance( YMBDistrib ) )
    @cached_property
    def _get_pdf_phi( self ):
        return YMBDistrib( varname = 'contact_fraction', data = self.data )

#    pd1_pu = PUniform( par_a = 0.003, par_b = .03, par_m = 300 )
#    pd2_pu = PUniform( par_a = .1, par_b = 3., par_m = 6.5 )
#    pd3_pu = PUniform( par_a = .02, par_b = 1., par_m = 30 )

    data_changed = Event( True )
    @on_trait_change( '+modified,model.+modified' )
    def _redraw( self ):
        s = SPIRRID( rf = rf,
    #                 min_eps = 0.00, max_eps = 20.00, n_eps = 80,
                min_eps = 0.00, max_eps = 0.3, n_eps = 100,
                cached_qg = True,
                compiled_qg_loop = False,
                compiled_eps_loop = False
                )
        # construct the random variables

        n_int = 30

        #s.add_rv( 'tau1', distribution='uniform', loc=2. / 1000, scale=10.0 / 1000, n_int=n_int )
        #s.add_rv( 'l', distribution='uniform', loc=0.01, scale=.01, n_int=n_int )
        #s.add_rv( 'xi', distribution='weibull_min', shape=8.64, scale=0.0287, n_int=n_int )

        s.rv_dict['theta'] = RV( pd = self.pdf_theta, name = 'theta', n_int = n_int )

        s.rv_dict['l'] = RV( pd = self.pdf_l, name = 'l', n_int = n_int )

        s.rv_dict['phi'] = RV( pd = self.pdf_phi, name = 'phi', n_int = n_int )

        s.mean_curve.plot( self.figure, linewidth = 2, label = 'mean CB per filament response' )

        plt.xlabel( 'crack opening w[mm]' )
        plt.ylabel( 'force P[N]' )
        plt.legend( loc = 'best' )

        plt.figure( 10 )

        plt.show()
