#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='either a fastq or fasta file, must end with .fasta, .fna, .fa, .fastq, or .fq')
parser.add_argument('-t', '--plotTitle')
parser.add_argument('-m', '--method', help='options: mean, median')
parser.add_argument('-ec', '--errorbarColor', help='desired error bar color in hex color, default is grey', default='grey') 
parser.add_argument('-r', '--range', help='Range setting for the color bar. Accepted arguments: num, perc, max. Num colors by the absolute number of reads ranging from 1k to 10k, perc colors by percentage of total, max colors based on minimum and maximum read counts', default='num')

args = parser.parse_args()

fastaEnds = ('.fasta', '.fna', '.fa')
fastqEnds = ('.fastq', '.fq')

if args.input.endswith(fastaEnds):
        gcContent = [GC(rec.seq) for rec in SeqIO.parse(args.input, "fasta")]
        lens = [len(rec) for rec in SeqIO.parse(args.input, "fasta")]
elif args.input.endswith(fastqEnds):
        gcContent = [GC(rec.seq) for rec in SeqIO.parse(args.input, "fastq")]
        lens = [len(rec) for rec in SeqIO.parse(args.input, "fastq")]
else:
        print("File extension not recognized, see help file")
print("Number of reads: %i" % len(gcContent))

df = pd.DataFrame({'length': lens, 'gcContent': gcContent})

dfGrouped = df.groupby(by='length').agg(['count', 'mean', 'median', 'std']).reset_index()

dfGrouped['perc'] = dfGrouped['gcContent', 'count']/dfGrouped['gcContent', 'count'].sum()

#limit options
plt.xlim(15, 75)
plt.ylim(20, 200)
plt.suptitle(args.plotTitle + "\n" + "n= " + str(dfGrouped['gcContent', 'count'].sum()))
plt.xlabel('GC content (%)')
plt.ylabel('Read length (bp)')

#plot dat
cm = plt.cm.get_cmap('rainbow')
plt.errorbar(dfGrouped['gcContent', args.method], dfGrouped['length', ''], xerr=dfGrouped['gcContent', 'std'], linestyle="None", marker="None", color=args.errorbarColor)

#color range options
if args.range == 'num':
	#plot by number of reads, range from 1k to 10k
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=1000, vmax=10000, marker='o', edgecolors='None', s=25, zorder=2)
elif args.range == 'perc':
	#plot by percentage instead
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['perc']), cmap=cm, vmin=0.0, vmax=1.0, marker='o', edgecolors='None', s=25, zorder=2)
elif args.range == 'max':
	#plot by min to max count?
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=min(dfGrouped['gcContent', 'count']), vmax=max(dfGrouped['gcContent', 'count']), marker='o', edgecolors='None', s=25, zorder=2)

plt.colorbar()
plt.draw()
plt.show() 
