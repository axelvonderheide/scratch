'''
Created on Dec 4, 2009

@author: jakub
'''
import numpy as np
import matplotlib.pyplot as plt

N = 5
#jakub#rch#rostar#alex#MrsB
orderMeans = np.array([44./5.+8.30+11.+10.3, 44./5.*4., 14.3, 12.5+10.8, 0.])
orderStd =   (2, 3, 4, 1, 2)

ind = np.arange(N)  # the x locations for the groups
width = 0.3       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(ind, orderMeans, width, color='r', yerr=orderStd)

creditMeans = np.array([4., 4., 2., 2., 5.],dtype=float)*3.3 + np.array([45.40,0.,10.,0.,0.])
creditStd =   (3., 5., 2., 3., 3.)
rects2 = ax.bar(ind+width, creditMeans, width, color='y', yerr=creditStd)

resultMeans = orderMeans - creditMeans
print 'result'
print 'jakub\t rch\t rostar \talex \tMrsB'
print resultMeans
resultStd = orderStd
rects3 = ax.bar(ind+2*width, resultMeans, width, color='b', yerr=resultStd)

# add some
ax.set_ylabel(r'EURO')
ax.set_title('Bewerage Optimizing Problem')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('Jakub', 'RCH', 'ROSTAR', 'Alex', 'MrsB') )

ax.legend( (rects1[0], rects2[0], rects3[0]), ('Order', 'Credit','Result') )

def autolabel(rects, values):
    # attach some text labels
    for  rect, val in zip(rects, values):
        height = rect.get_height()
        if val > 0:
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%0.1f'%val,
                    ha='center', va='bottom')
        else:
            ax.text(rect.get_x()+rect.get_width()/2., 0.4, '%0.1f'%val,
                    ha='center', va='bottom')
autolabel(rects1, orderMeans)
autolabel(rects2, creditMeans)
autolabel(rects3, resultMeans)
plt.show()
