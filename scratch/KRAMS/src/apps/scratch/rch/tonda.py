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
# Created on Nov 26, 2010 by: rch

from string import index

jair_sorted = sorted( ['jana', 'antonin', 'ifca', 'rosta' ] )

print jair_sorted

inicialy = ['j', 'a', 'i', 'r']

print '------------ spravne: -------------'
for ini in inicialy:
    for name in jair_sorted:
        if ini == name[0]:
            print name
            break


print 'kolik je 23% ze 4563.0?'
print 4563.0 * 0.23

from time import sleep

pocet_vterin = 2000
pocet_minut = pocet_vterin / 60
print '==================================='
print 'bode to trvat', pocet_minut, 'minut'
print '==================================='


for i in range( 1, 2000 ):
    sleep( 1.0 )
    print 'probudil jsem se', i
