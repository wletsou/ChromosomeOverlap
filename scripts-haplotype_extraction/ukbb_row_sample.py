#!/usr/bin/python

#randomly select n lines from a combined genotype file and output a cases file (n lines) and controls file (total - n lines)

import argparse, random, os, re, sys
import pandas as pd
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("-n","--Sample", help = "Number of lines to draw",type=int)
parser.add_argument("-f","--File", help = "File to sample from")
parser.add_argument("-s","--Seed", help = "Random seed",type = int,default = 20200116)
parser.add_argument("-no_header","--No_Header", help = "Header line?",nargs='?',const = False,default = True)
parser.add_argument("-no_row_names","--No_Row_Names", help = "Is there a column of IDs?",nargs='?',const = False,default = True)
parser.add_argument("-d","--Directory", help = "Output directory",default = os.getcwd())
parser.add_argument("-r","--Repetition", help = "Sample lines with repetiton?",nargs='?',const = True,default = False)
args = parser.parse_args()

n = args.Sample # number of columns, not including index column
seed = args.Seed
directory = args.Directory
header = args.No_Header
if header:
    header = 0 # set first row as header
else:
    header = None
row_names = args.No_Row_Names
if row_names:
    row_names = 0 # set first column as index column
repetition = args.Repetition

name_head = re.sub('(.*)[.].*',r"\1",re.sub('.*/(.*)',r"\1",str(args.File))) # greedy match, get file name up to extension
name_ext = re.sub('.*[.](.*)',r"\1",str(args.File)) # get file extension

n_lines=len(open(args.File).readlines())
if (n>n_lines) and not repetition:
    n=n_lines-args.No_Header
df = pd.read_csv(args.File, delimiter = "\t", index_col = row_names, header = header)
df_sample = df.sample(n = n,axis = 'index',replace = repetition, random_state = seed) # sample n random columns, unordered

if header == 0:
    header = True
else:
    header = False

if row_names == 0:
    row_names = True
else:
    row_names = False

print("{0}/{1}.row_sample_{2}.{3}\n".format(directory,name_head,n,name_ext))
df_sample.to_csv("{0}/{1}.row_sample_{2}.{3}".format(directory,name_head,n,name_ext),sep = "\t",header = header,index = row_names)
