'''
Created on Oct 14, 2011

@author: rostar
'''
import numpy as np


lengths = [50, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500]
test = ['Statimat', 'resin']
n_tests = 30

# get the unique test combinations
combinations = []
for t in test:
    for l in lengths:
        combinations.append((t, l))

# multiply them by the number of replications
ordered_tests = n_tests * combinations

# randomize the tests
np.random.shuffle(ordered_tests)

DOE = ordered_tests

print DOE
