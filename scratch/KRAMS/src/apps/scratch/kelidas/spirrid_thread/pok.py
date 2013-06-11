'''
Created on 19.10.2010

@author: Vasek
'''
from threading import Thread
from numpy import sqrt, exp


def thread( c, n ):
    ret = 0
    for i in range( 0, n ):
        ret = 3 * sqrt( i * n * n )
    print 'thread', c, 'finished'
    print ret

thr = []
for i in range( 1, 4 ):
    thr.append( Thread( target=thread,
                    args=( [i, i * 1000 - i * 100] ) ) )


for i in ( 3, 1, 2 ):
    thr[i - 1].start()
    

