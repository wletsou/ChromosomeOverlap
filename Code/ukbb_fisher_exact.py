#!/usr/bin/python

import argparse, random, os, re, sys
from operator import itemgetter
import numpy as np
import scipy.stats as stats

parser = argparse.ArgumentParser()
parser.add_argument("-f","--File", help = "File to sample from")
parser.add_argument("-p","--Population", help = "Population names")
parser.add_argument("-d","--Directory", help = "Output directory",default = os.getcwd())
args = parser.parse_args()

file = args.File
directory = args.Directory
name_str = re.sub('(.*?)[.](.*)',"fisher_exact."r"\2",re.sub('.*/(.*)',r"\1",file))
f = open('{}/{}'.format(directory,name_str),"w") # initialize output file
f.close()
with open(file,'rt') as myfile:
    for line in myfile:
        elements=line.split()
        hap=elements[0]
        a=float(elements[1]) # cases counts
        b=float(elements[2]) # cases non-counts
        c=float(elements[3]) # controls counts
        d=float(elements[4]) # controls non-counts
        n=a+b # total cases
        m=c+d # total controls
        freq_cases=a/n
        freq_controls=c/m
        oddsratio,p_value=stats.fisher_exact([[a,b],[c,d]],'two-sided')
        str="{}\t{:.6e}\t{:.6e}\t{:.6f}\t{:.6e}\n".format(hap,freq_cases,freq_controls,oddsratio,p_value)
        with open('{}/{}'.format(directory,name_str),"a") as outfile:
            outfile.write(str)
