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
# Created on Jan 30, 2010 by: rch

from enthought.pyface.api import ProgressDialog
import time

def task_func(t):
    progress = ProgressDialog(title="progress", message="counting to %d"%t,
                              max=t, show_time=True, can_cancel=True)
    progress.open()

    for i in range(0,t+1):
        time.sleep(1)
        print i
        (cont, skip) = progress.update(i)
        if not cont or skip:
            break

    progress.update(t)
    
task_func(10)
