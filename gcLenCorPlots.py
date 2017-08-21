#!/usr/bin/python

'''Usage: python gcLenCorPlots.py -i <input fasta or fastq> [-m <method> -r <range for heatmap> -t <trim maximum length> -s <normalize to number> -ec <error bar color>]'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse
import scipy
from scipy.stats import mannwhitneyu
from scipy.signal import resample

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='either a fastq or fasta file, must end with .fasta, .fna, .fa, .fastq, or .fq')
parser.add_argument('-m', '--method', help='options: mean, median', default='mean')
parser.add_argument('-ec', '--errorbarColor', help='desired error bar color in hex color, default is grey', default='grey') 
parser.add_argument('-r', '--range', help='Range setting for the color bar. Accepted arguments: num, perc, max. Num colors by the absolute number of reads ranging from 100 to 2k, perc colors by percentage of total, max colors based on minimum and maximum read counts', default='max')
parser.add_argument('-t', '--trim', help='Maximum length trim, numeric')
parser.add_argument('-s', '--shuffle', help='Randomly shuffle results to a specific number')

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

#optional trimming/shuffle options
if args.trim is not None:
	df = df.drop(df[df.length > int(args.trim)].index).reset_index()
	print("Length trimmed to maximum %i" % int(args.trim))
	print("Number of trimmed reads: %i" % len(df.gcContent))

if args.shuffle is not None:
	df = df.sample(int(args.shuffle))
	print("Number of reads normalized to %i" % int(args.shuffle))

#stats
print("Mean GC content: %.2f" % np.mean(df['gcContent']))
print("Median GC content: %.2f" % np.median(df['gcContent']))
print("Mean fragment length: %i" % np.mean(df['length']))
print("Median fragment length: %i" % np.median(df['length']))
print("Fragment length range: %i : %i" % (min(df['length']), max(df['length'])))
print("GC content range: %.2f : %.2f" % (min(df['gcContent']), max(df['gcContent'])))	

#grouped data
dfGrouped = df.groupby(by='length').agg(['count', 'mean', 'median', 'std']).reset_index()
dfGrouped['perc'] = dfGrouped['gcContent', 'count']/dfGrouped['gcContent', 'count'].sum()

#print out data
with open('%s_data_out.txt' % args.input, 'w') as outfile:
	dfGrouped.to_csv(outfile, sep="\t", index=False)

#calculate significance of each grouping, compared to overall distribution
dfLensGroup = df.groupby(by='length')
overall = df['gcContent']
lines = []

#plotting

#limit options
plt.xlim(15, 90)
plt.ylim(20, 200)
plt.suptitle(args.input + "\n" + "n= " + str(dfGrouped['gcContent', 'count'].sum()))
plt.xlabel('GC content (%)')
plt.ylabel('Read length (bp)')

#plot dat
cm = plt.cm.get_cmap('YlOrRd')
plt.errorbar(dfGrouped['gcContent', args.method], dfGrouped['length', ''], xerr=dfGrouped['gcContent', 'std'], linestyle="None", marker="None", color=args.errorbarColor)

#color range options
if args.range == 'num':
	#plot by number of reads, range from 1k to 10k
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=100, vmax=2000, marker='o', edgecolors='None', s=25, zorder=2)
elif args.range == 'perc':
	#plot by percentage instead
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['perc']), cmap=cm, vmin=0.0, vmax=1.0, marker='o', edgecolors='None', s=25, zorder=2)
elif args.range == 'max':
	#plot by min to max count?
	plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=min(dfGrouped['gcContent', 'count']), vmax=max(dfGrouped['gcContent', 'count']), marker='o', edgecolors='None', s=25, zorder=2)

plt.colorbar()
plt.draw()
plt.savefig('%s_plot.pdf' % args.input)

