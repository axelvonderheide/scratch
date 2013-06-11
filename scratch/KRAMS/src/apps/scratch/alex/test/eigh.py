from numpy import *
from scipy.linalg import *
from time import *


print 'int', int(2.713)
print 'XXX', log(2.713)

'''def eigh(a, lower=True, eigvals_only=False, overwrite_a=False):
    Solve real symmetric or complex hermitian eigenvalue problem.

    Inputs:

      a            -- A hermitian N x N matrix.
      lower        -- values in a are read from lower triangle
                      [True: UPLO='L' (default) / False: UPLO='U']
      eigvals_only -- don't compute eigenvectors.
      overwrite_a  -- content of a may be destroyed

    Outputs:

      For eigvals_only == False (the default),
      w,v     -- w: eigenvalues, v: eigenvectors
      For eigvals_only == True,
      w       -- eigenvalues

    Definitions:

      a * v[:,i] = w[i] * vr[:,i]
      v.H * v = identity


def eig(a,b=None, left=False, right=True, overwrite_a=False, overwrite_b=False):
    Solve ordinary and generalized eigenvalue problem
    of a square matrix.

    Inputs:

      a     -- An N x N matrix.
      b     -- An N x N matrix [default is identity(N)].
      left  -- Return left eigenvectors [disabled].
      right -- Return right eigenvectors [enabled].
      overwrite_a, overwrite_b -- save space by overwriting the a and/or
                                  b matrices (both False by default)

    Outputs:

      w      -- eigenvalues [left==right==False].
      w,vr   -- w and right eigenvectors [left==False,right=True].
      w,vl   -- w and left eigenvectors [left==True,right==False].
      w,vl,vr  -- [left==right==True].

    Definitions:

      a * vr[:,i] = w[i] * b * vr[:,i]

      a^H * vl[:,i] = conjugate(w[i]) * b^H * vl[:,i]

    where a^H denotes transpose(conjugate(a)).
    '''

a = arange(0)

A12 = 1.e-15
A = array([[1.,A12],[A12,1.]])
#A = array([[1.,2],[2,3.]])
print 'A',A
nnn = 10000

t1 = time()
for i in range(nnn):
    ev, EV = eig(A)
t2 = time()
t_eig = t2-t1
print 't_eig',t_eig

print 'ev',ev
print 'EV',EV

t1 = time()
for i in range(nnn):
    evh, EVh = eigh(A)
t2 = time()
t_eigh = t2-t1
print 't_eigh',t_eigh

print 'evh',evh
print 'EVh',EVh

'''
'eigh' is twice as fast as 'eig'
'''

# ---
#print '\n'
#
#b = arange(0)
#
#B12 = 500000#1.e-15
#B21 = 1.e-2
#
## Inputs are read only from the lower triangle of B !!!
#
#B = array([[1.,B12],[B21,1.]])
##A = array([[1.,2],[2,3.]])
#
#
#
#print 'B',B
#ev, EV = eig(B)
#print 'ev',ev
#print 'EV',EV
#
#evh, EVh = eigh(B)
#print 'evh',evh
#print 'EVh',EVh
