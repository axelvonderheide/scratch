
import numpy as np
import scipy.weave as weave

import platform
if platform.system() == 'Linux':
    from time import time as sysclock
elif platform.system() == 'Windows':
    from time import clock as sysclock

l = 5000
x = np.linspace(0., 10., l)
y = np.linspace(0., 10., l)

def f(x, y):
    return np.sum(x[:, None] * y[None, :])

start = sysclock()
print f(x, y)
print sysclock() - start, 'numpy'

c_code = '''
double r = 0;
for( int i_x=0; i_x < l; i_x++){
    for( int i_y=0; i_y < l; i_y++){
        r += x(i_x) * y(i_y);// nefunguje *(x + i_x) * *(y + i_y);
    };
};
return_val = r;
'''

conv = weave.converters.blitz
#conv = weave.converters.default
compiler_verbose = 0
compiler = 'gcc'


arg_names = ['x', 'y', 'l']
arg_values = {'x':x, 'y':y, 'l':l}

start = sysclock()
print weave.inline(c_code, arg_names,
        local_dict = arg_values,
        type_converters = conv,
        compiler = compiler,
        verbose = compiler_verbose)
print sysclock() - start, 'c'





