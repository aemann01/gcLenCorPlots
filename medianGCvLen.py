#!/usr/bin/python

'''Used in DVC project to plot median GC v length values for each individual sample'''

import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from scipy.interpolate import spline
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--type', help='file type, must either be fasta or fastq')
parser.add_argument('-o', '--out', help='output file name')

args = parser.parse_args()

#files = glob('*.human.merged.fastq')
fig, ax = plt.subplots()

files = ['F349Aroot.human.merged.fastq', 'S454root.human.merged.fastq', 'NF217root.human.merged.fastq', 'C214root.human.merged.fastq', 'S41broot.human.merged.fastq', 'H24broot.human.merged.fastq', 'H10broot.human.merged.fastq', 'F1948root.human.merged.fastq', 'S37root.human.merged.fastq', 'S108root.human.merged.fastq', 'cmol53root.human.merged.fastq', 'NF47root.human.merged.fastq', 'NF47calc.human.merged.fastq', 'S108calc.human.merged.fastq', 'H10bcalc.human.merged.fastq', 'S454calc.human.merged.fastq', 'C214calc.human.merged.fastq', 'S37calc.human.merged.fastq', 'H24bcalc.human.merged.fastq', 'F349Acalc.human.merged.fastq', 'F1948calc.human.merged.fastq', 'S41bcalc.human.merged.fastq', 'NF217calc.human.merged.fastq', 'cmol53calc.human.merged.fastq']

print(files)

colors = ('#000000','#000000','#000000','#000000','#000000','#000000','#000000','#000000','#000000','#000000','#000000','#000000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000',)
#colors = ('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#808080', '#b15928','#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#808080', '#b15928')
i = 0

for f in files:
	print("Reading in %s" %f)
	gcContent = [GC(rec.seq) for rec in SeqIO.parse(f, args.type)]
	lens = [len(rec) for rec in SeqIO.parse(f, args.type)]
	print("Number of reads: %i" % len(gcContent))
	
	df = pd.DataFrame({'length':lens, 'gcContent':gcContent})
	#print(df)
	dfGrouped = df.groupby(by='length').agg(['count','median']).reset_index()
	#print(dfGrouped)
	#print(dfGrouped['gcContent', 'median'])
#version below needs testing
#	plt.plot(pd.Series.rolling(center=False, window=25).mean(dfGrouped['gcContent', 'median']), pd.Series.rolling(center=False, window=25).mean(dfGrouped['length', '']), label=f.replace(args.type, ""), linewidth=1.0, color=colors[i])
	plt.plot(pd.rolling_mean(dfGrouped['gcContent', 'median'], window=25), pd.rolling_mean(dfGrouped['length', ''], window=25), label=f.replace(args.type, ""), linewidth=1.0, color=colors[i])
	i += 1
	print('--')

#plt.legend(loc='upper right')
plt.xlim(25, 75)
plt.ylim(20,200)
plt.xlabel('GC content (%)')
plt.ylabel('Read length (bp)')
plt.savefig(args.out)
