#!/usr/bin/python
#
###########################################################################
#                                                                         #
# Script to plot a CCall vs. CCall from SHELXD.                           #
#                                                                         #    
#                                                                         #   
#                                                                         #  
# Usage:   shelxd_cc.py                                                   #  
#                                                                         #  
# AW 01/05/2016                                                           #  
#                                                                         #  
###########################################################################

import sys, os, string, shutil
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import heapq

CCall = []
CCweak = []
CFOM = []
PATFOM = []

# get filename
filename = sys.argv[-1]
inputFile = open(filename)

for line in inputFile.readlines():
        bla = line.replace('U', ' ')
        bla2 = bla.replace('/', ' ')
        test = bla2.replace(',', ' ')
        column = test.split() 
	if column:
		if column[0] == 'Try':
			CCall.append(float(column[7]))
			CCweak.append(float(column[8]))
			CFOM.append(float(column[10]))
#			BEST.append(float(column[12]))

median = np.median(CCall)
mad = np.median(np.sqrt((CCall - median)**2))

CCall_max = heapq.nlargest(3, CCall)
CCall_mad = CCall - median

mad10 = sum(i > 10 * mad for i in CCall_mad)
mad9 = sum(i > 9 * mad for i in CCall_mad)
mad8 = sum(i > 8 * mad for i in CCall_mad)
mad7 = sum(i > 7 * mad for i in CCall_mad)
mad6 = sum(i > 6 * mad for i in CCall_mad)
mad5 = sum(i > 5 * mad for i in CCall_mad)


# print "Median=", median
# print "MAD=", mad
# print "3 largest CCall values", CCall_max

print "number of CCall with CCall - median > 10 * MAD", mad10
print "number of CCall with CCall - median > 9 * MAD", mad9
print "number of CCall with CCall - median > 8 * MAD", mad8
print "number of CCall with CCall - median > 7 * MAD", mad7
print "number of CCall with CCall - median > 6 * MAD", mad6
print "number of CCall with CCall - median > 5 * MAD", mad5

# next step is to create a histogram of bins of 0.5 * mad in size 

# print CCall_mad10
# print CCall_mad8

# print np.sqrt((CCall - median)**2)

# print median 
# print med_abs_deviation
# print max(CCall)

# next step is to count the number of CCall values > 5 MAD


	
plt.plot(CCweak, CCall, "o")
plt.title("SHELXD CCall/CCweak")
plt.xlabel("CCweak")
plt.ylabel("CCall")


fig = plt.gcf()
fig.savefig('CC.png', format='png')
# plt.show()
plt.close(fig)
