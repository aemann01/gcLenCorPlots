#!/usr/bin/python

'''Used in DVC project to plot median GC v length values for each individual sample'''

import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from scipy.interpolate import spline


files = glob('*calc.fasta')
fig, ax = plt.subplots()

colors = ('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#808080', '#b15928')
i = 0

for f in files:
	print("Reading in %s" %f)
	gcContent = [GC(rec.seq) for rec in SeqIO.parse(f, "fasta")]
	lens = [len(rec) for rec in SeqIO.parse(f, "fasta")]
	print("Number of reads: %i" % len(gcContent))
	
	df = pd.DataFrame({'length':lens, 'gcContent':gcContent})
	#print(df)
	dfGrouped = df.groupby(by='length').agg(['count','median']).reset_index()
	#print(dfGrouped)
	#print(dfGrouped['gcContent', 'median'])
	plt.plot(pd.rolling_mean(dfGrouped['gcContent', 'median'], window=25), pd.rolling_mean(dfGrouped['length', ''], window=25), label=f.replace(".fasta", ""), linewidth=1.0, color=colors[i])
	i += 1
	print(i)

plt.legend(loc='upper right')
plt.xlim(25, 75)
plt.ylim(20,200)
plt.xlabel('GC content (%)')
plt.ylabel('Read length (bp)')
plt.savefig('medianCalcGCvLen_bySample.pdf')
