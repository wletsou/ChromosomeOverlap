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
args = parser.parse_args()

n1 = float(args.n1)
n2 = float(args.n2)
x1 = float(args.x1)
x2 = float(args.x2)
y1 = n1 - x1
y2 = n2 - x2

oddsratio,p_value=stats.fisher_exact([[x1,y1],[x2,y2]],'two-sided')

print("p = %0.6E" % p_value) # p-value for difference in proportions
