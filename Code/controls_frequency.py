#!/usr/bin/python

import argparse, random, os, re, sys
from operator import itemgetter
import numpy as np
import scipy.stats as stats

parser = argparse.ArgumentParser()
parser.add_argument("-x1","--x1", help = "Count in population 1")
parser.add_argument("-n1","--n1", help = "Total in population 1")
parser.add_argument("-x2","--x2", help = "Count in population 2")
parser.add_argument("-n2","--n2", help = "Total in population 2")
parser.add_argument("-pooled","--Pooled", action = 'store_true',help = "Pooled variance",default = False)
args = parser.parse_args()

n1 = float(args.n1)
n2 = float(args.n2)
x1 = float(args.x1)
x2 = float(args.x2)
pooled = args.Pooled
p1 = x1/n1
p2 = x2/n2

if not pooled:
    se2_1 = p1*(1-p1)/n1
    se2_2 = p2*(1-p2)/n2
    z = (p1-p2)/np.sqrt(se2_1+se2_2) # unpooled
else:
    p = (x1+x2)/(n1+n2)
    z = (p1-p2)/np.sqrt(p*(1-p)*((1/n1)+(1/n2))) # pooled
    
print("p = %0.6E" % (2*np.minimum(stats.norm.cdf(-z),stats.norm.cdf(z)))) # p-value for difference in proportions
