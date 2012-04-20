#!/usr/bin/env python

# Mayor Library imports
from numpy import linspace, pi, sin, cos

from enthought.mayavi.core.source import Source

# Enthought library imports
from enthought.chaco.shell import plot, hold, title, show

# Create some data
x = linspace(-2*pi, 2*pi, 100)
y = sin(x)

# Create some line plots
plot(x, y, "b-", bgcolor="white")

#This command is only necessary if running from command line
#show()
