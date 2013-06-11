'''
Created on Jun 8, 2009

@author: alexander
'''

from scipy.optimize import newton, fsolve
from time import *


func = lambda x: (x-1)**2 - 1
#starting estimate:
x0 = 0.99

t_a = 0.0
t_b = 0.0

xtol=1.49012e-8

for n in range(1000):
    t_1 = time()
    a = newton( func,x0, tol=1.49012e-8)
    t_2 = time()
    t_a += t_2 - t_1

    t_1 = time()
    b = fsolve( func,x0, xtol=1.49012e-8 )
    t_2 = time()
    t_b += t_2 - t_1

print t_a
print t_b

print a
print b


#fsolve(func,x0,args=(),fprime=None,full_output=0,col_deriv=0,xtol=1.49012e-8,maxfev=0,band=None,epsfcn=0.0,factor=100,diag=None, warning=True)











#----------------
#def fsolve(func,x0,args=(),fprime=None,full_output=0,col_deriv=0,xtol=1.49012e-8,maxfev=0,band=None,epsfcn=0.0,factor=100,diag=None, warning=True):
#----------------
"""Find the roots of a function.

  Description:

    Return the roots of the (non-linear) equations defined by
    func(x)=0 given a starting estimate.

  Inputs:

    func -- A Python function or method which takes at least one
            (possibly vector) argument.
    x0 -- The starting estimate for the roots of func(x)=0.
    args -- Any extra arguments to func are placed in this tuple.
    fprime -- A function or method to compute the Jacobian of func with
            derivatives across the rows. If this is None, the
            Jacobian will be estimated.
    full_output -- non-zero to return the optional outputs.
    col_deriv -- non-zero to specify that the Jacobian function
                 computes derivatives down the columns (faster, because
                 there is no transpose operation).
    warning -- True to print a warning message when the call is
                unsuccessful; False to suppress the warning message.
  Outputs: (x, {infodict, ier, mesg})

    x -- the solution (or the result of the last iteration for an
         unsuccessful call.

    infodict -- a dictionary of optional outputs with the keys:
                'nfev' : the number of function calls
                'njev' : the number of jacobian calls
                'fvec' : the function evaluated at the output
                'fjac' : the orthogonal matrix, q, produced by the
                         QR factorization of the final approximate
                         Jacobian matrix, stored column wise.
                'r'    : upper triangular matrix produced by QR
                         factorization of same matrix.
                'qtf'  : the vector (transpose(q) * fvec).
    ier -- an integer flag.  If it is equal to 1 the solution was
           found.  If it is not equal to 1, the solution was not
           found and the following message gives more information.
    mesg -- a string message giving information about the cause of
            failure.

  Extended Inputs:

   xtol -- The calculation will terminate if the relative error
           between two consecutive iterates is at most xtol.
   maxfev -- The maximum number of calls to the function. If zero,
             then 100*(N+1) is the maximum where N is the number
             of elements in x0.
   band -- If set to a two-sequence containing the number of sub-
           and superdiagonals within the band of the Jacobi matrix,
           the Jacobi matrix is considered banded (only for fprime=None).
   epsfcn -- A suitable step length for the forward-difference
             approximation of the Jacobian (for fprime=None). If
             epsfcn is less than the machine precision, it is assumed
             that the relative errors in the functions are of
             the order of the machine precision.
   factor -- A parameter determining the initial step bound
             (factor * || diag * x||). Should be in interval (0.1,100).
   diag -- A sequency of N positive entries that serve as a
           scale factors for the variables.

  Remarks:

    "fsolve" is a wrapper around MINPACK's hybrd and hybrj algorithms.

  See also:

      scikits.openopt, which offers a unified syntax to call this and other solvers

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar and vector fixed-point finder

    """ 
# Netwon-Raphson method:
#----------------
#def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
#----------------
"""Given a function of a single variable and a starting point,
find a nearby zero using Newton-Raphson.

fprime is the derivative of the function.  If not given, the
Secant method is used.

See also:

  fmin, fmin_powell, fmin_cg,
         fmin_bfgs, fmin_ncg -- multivariate local optimizers
  leastsq -- nonlinear least squares minimizer

  fmin_l_bfgs_b, fmin_tnc,
         fmin_cobyla -- constrained multivariate optimizers

  anneal, brute -- global optimizers

  fminbound, brent, golden, bracket -- local scalar minimizers

  fsolve -- n-dimenstional root-finding

  brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

  fixed_point -- scalar and vector fixed-point finder

"""
