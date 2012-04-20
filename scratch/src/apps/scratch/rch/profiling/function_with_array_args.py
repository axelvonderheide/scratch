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
# Created on Sep 28, 2011 by: rch

import numpy as np

def f(a, b, c, d):
    e = a * b
    f = e * c
    g = f * d
    return g

def run1():
    a = 0.
    b = 0.
    c = np.zeros((1000, 1), dtype = float)
    d = np.zeros((1, 1000), dtype = float)
    for i in range(1000):
        f(a, b, c, d)

def run2():
    a = np.zeros((1, 1, 1, 1), dtype = float)
    b = np.zeros((1, 1, 1, 1), dtype = float)
    c = np.zeros((1, 1, 1000, 1), dtype = float)
    d = np.zeros((1, 1, 1, 1000), dtype = float)
    for i in range(1000):
        f(a, b, c, d)

if __name__ == '__main__':

    from time import time as sysclock

    start_time = sysclock()
    run1()
    exec_time = sysclock() - start_time
    print '1', exec_time

    start_time = sysclock()
    run2()
    exec_time = sysclock() - start_time
    print '2', exec_time
