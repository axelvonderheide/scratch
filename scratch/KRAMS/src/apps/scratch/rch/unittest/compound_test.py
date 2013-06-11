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
# Created on Aug 18, 2009 by: rch

import unittest
import sys

def test_generator(func):
    def inner(self):
        failures = []
        errors = []
        for test, args in func(self):
            try:
                test(*args)
            except self.failureException, e:
                failures.append((test.__name__, args, e))
            except KeyboardInterrupt:
                raise
            except:
                # using sys.exc_info means we also catch string exceptions
                e = sys.exc_info()[1]
                errors.append((test.__name__, args, e))
        msg = '\n'.join('%s%s: %s: %s' % (name, args, e.__class__.__name__, e) for (name, args, e) in failures + errors)
        if errors:
            raise Exception(msg)
        raise self.failureException(msg)
    return inner


class Test2(unittest.TestCase):

    @test_generator
    def testSomething(self):
        for a, b in ((1, 2), (3, 3), (5, 4)):
            yield self.assertEqual, (a, b)

        def raises():
            raise Exception('phooey')
        yield raises, ()