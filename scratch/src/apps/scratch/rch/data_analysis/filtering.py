
from numpy import cumsum
from pylab import *

def histOutline(dataIn, *args, **kwargs):
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

# Make some data to plot
data = randn(500)

figure(2, figsize=(10, 5))
clf()

##########
#
# First make a normal histogram
#
##########
subplot(1, 2, 1)
(n, bins, patches) = hist(data)

# Boundaries
xlo = -max(abs(bins))
xhi = max(abs(bins))
ylo = 0
yhi = max(n) * 1.1

axis([xlo, xhi, ylo, yhi])

##########
#
# Now make a histogram in outline format
#
##########
(bins, n) = histOutline(data)

subplot(1, 2, 2)
#plot(bins, n, 'k-')
plot(bins, cumsum(n))
#axis([xlo, xhi, ylo, yhi])
show()
