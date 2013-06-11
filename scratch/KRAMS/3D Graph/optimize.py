from scipy.optimize import minimize
import numpy as np
def func(x, c):
    #get
    return 5 - (x - 2) ** (c[0] * c[1])

def getarg():
    return 1, 2, 3

c = [1, 2]

max_x = minimize(func, 2.0, args=([c]))
print max_x.fun
