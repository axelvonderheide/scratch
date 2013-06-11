'''
Created on Sep 3, 2012

@author: rch
'''

from ibvpy.api import IBVModel

from spirrid import SPIRRID

if __name__ == '__main__':
    q = lambda eps, xi: eps * xi
    s = SPIRRID(q = q,
                eps_vars = {'eps' : [0, 1]},
                theta_vars = {'xi' : 2.0 })
    print s.mu_q_arr

