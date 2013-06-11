'''target curve for the fitting procedure of the phi-function:
Data is obtained from the MATS2D_cmdm-model for quasibrittle, isotropic
behaviour using an exponential softening law 
'''
from numpy import array

# first principle strains:
xdata = array([\
                      0.00000000e+00,   3.33333333e-05,   6.66666667e-05,   1.00000000e-04,
                      1.33333333e-04,   1.66666667e-04,   2.00000000e-04,   2.33333333e-04,
                      2.66666667e-04,   3.00000000e-04,   3.33333333e-04,   3.66666667e-04,
                      4.00000000e-04,   4.33333333e-04,   4.66666667e-04,   5.00000000e-04,
                      5.33333333e-04,   5.66666667e-04,   6.00000000e-04,   6.33333333e-04,
                      6.66666667e-04,   7.00000000e-04,   7.33333333e-04,   7.66666667e-04,
                      8.00000000e-04,   8.33333333e-04,   8.66666667e-04,   9.00000000e-04,
                      9.33333333e-04,   9.66666667e-04,   1.00000000e-03,   1.03333333e-03,
                      1.06666667e-03,   1.10000000e-03,   1.13333333e-03,   1.16666667e-03,
                      1.20000000e-03,   1.23333333e-03,   1.26666667e-03,   1.30000000e-03,
                      1.33333333e-03,   1.36666667e-03,   1.40000000e-03,   1.43333333e-03,
                      1.46666667e-03,   1.50000000e-03,   1.53333333e-03,   1.56666667e-03,
                      1.60000000e-03,   1.63333333e-03,   1.66666667e-03,   1.70000000e-03,
                      1.73333333e-03,   1.76666667e-03,   1.80000000e-03,   1.83333333e-03,
                      1.86666667e-03,   1.90000000e-03,   1.93333333e-03,   1.96666667e-03,
                      2.00000000e-03])

# first principle stresses:
ydata = array([\
                      0.00000000e+00,   1.18055556e+00,   2.29430105e+00,   2.36326234e+00,
                      2.28213105e+00,   2.03030559e+00,   1.76395690e+00,   1.53034141e+00,
                      1.32610310e+00,   1.14796769e+00,   9.92883560e-01,   8.58069799e-01,
                      7.41025452e-01,   6.39520381e-01,   5.51577665e-01,   4.75452440e-01,
                      4.09609802e-01,   3.52703149e-01,   3.03553673e-01,   2.61131378e-01,
                      2.24537735e-01,   1.92990006e-01,   1.65807184e-01,   1.42397435e-01,
                      1.22246961e-01,   1.04910126e-01,   9.00007553e-02,   7.71844765e-02,
                      6.61719893e-02,   5.67131712e-02,   4.85919180e-02,   4.16216370e-02,
                      3.56413146e-02,   3.05120908e-02,   2.61142765e-02,   2.23447616e-02,
                      1.91147634e-02,   1.63478720e-02,   1.39783572e-02,   1.19497004e-02,
                      1.02133245e-02,   8.72749509e-03,   7.45637085e-03,   6.36918287e-03,
                      5.43952650e-03,   4.64475028e-03,   3.96542899e-03,   3.38490960e-03,
                      2.88891336e-03,   2.46523247e-03,   2.10338770e-03,   1.79439850e-03,
                      1.53058962e-03,   1.30539231e-03,   1.11318615e-03,   9.49163393e-04,
                      8.09212363e-04,   6.89817531e-04,   5.87973718e-04,   5.01112562e-04,
                      4.27039468e-04])


if __name__ == '__main__':
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
    from mathkit.mfn.mfn_line.mfn_matplotlib_editor import \
        MFnMatplotlibEditor
    from enthought.traits.api import Instance, HasTraits
    from enthought.traits.ui.api import View, Item
    class MFn(HasTraits):
        mfn = Instance( MFnLineArray )
        def _mfn_default(self):
            return MFnLineArray( xdata = xdata, ydata = ydata )
        view = View( Item('mfn', editor = MFnMatplotlibEditor() )
                     )

    mfn = MFn()
    mfn.configure_traits()