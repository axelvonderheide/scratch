'''
Created on 19.10.2010

@author: Vasek
'''

from multiprocessing import Process
from numpy import sqrt, exp


def thread( c, n ):
    ret = 0
    for i in range( 0, n * 1000 ):
        ret = 3 * sqrt( i * n * n )
    print 'thread', c, 'finished'
    print ret

if __name__ == '__main__':

    thr = []
    for i in range( 1, 4 ):
        thr.append( Process( target=thread,
                        args=( [i, i * 1000 - i * 100] ) ) )
        
    for i in ( 3, 2, 1 ):
        thr[i - 1].start()


