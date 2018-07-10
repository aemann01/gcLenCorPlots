#!/usr/bin/env python3

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''Usage: python gcLenCorPlots.py -i <input fasta or fastq> [-m <method> -r <range for heatmap> -t <trim maximum length> -s <normalize to number> -ec <error bar color>]'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import scipy
from scipy.stats import mannwhitneyu
from scipy.signal import resample
from Bio import SeqIO
from Bio.SeqUtils import GC

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='either a fastq or fasta file, must end with .fasta, .fna, .fa, .fastq, or .fq')
parser.add_argument('-m', '--method', help='options: mean, median', default='mean')
parser.add_argument('-ec', '--errorbarColor', help='desired error bar color in hex color, default is grey', default='grey') 
parser.add_argument('-r', '--range', help='Range setting for the color bar. Accepted arguments: num, perc, max. Num colors by the absolute number of reads ranging from 100 to 2k, perc colors by percentage of total, max colors based on minimum and maximum read counts', default='max')
parser.add_argument('-t', '--trim', help='Maximum length trim, numeric')
parser.add_argument('-s', '--shuffle', help='Randomly shuffle results to a specific number')

args = parser.parse_args()

def openFastx():
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
	printStats(df)

def printStats(df):
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
	plotDat(dfGrouped)

def plotDat(dfGrouped):
with open('%s_data_out.txt' % args.input, 'w') as outfile:
	dfGrouped.to_csv(outfile, sep="\t", index=False)
#plotting
	plt.xlim(15, 90)
	plt.ylim(20, 200)
	plt.suptitle(args.input + "\n" + "n= " + str(dfGrouped['gcContent', 'count'].sum()))
	plt.xlabel('GC content (%)')
	plt.ylabel('Read length (bp)')
	cm = plt.cm.get_cmap('YlOrRd')
	plt.errorbar(dfGrouped['gcContent', args.method], dfGrouped['length', ''], xerr=dfGrouped['gcContent', 'std'], linestyle="None", marker="None", color=args.errorbarColor)
#color range options
	if args.range == 'num':
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=100, vmax=2000, marker='o', edgecolors='None', s=25, zorder=2)
	elif args.range == 'perc':
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['perc']), cmap=cm, vmin=0.0, vmax=1.0, marker='o', edgecolors='None', s=25, zorder=2)
	elif args.range == 'max':
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=min(dfGrouped['gcContent', 'count']), vmax=max(dfGrouped['gcContent', 'count']), marker='o', edgecolors='None', s=25, zorder=2)
	plt.colorbar()
	plt.draw()
	plt.savefig('%s_plot.pdf' % args.input)

def main():
	print("Author: Allison E. Mann (allison.mann@botany.ubc.ca)")
	print("Copyright 2018: GPLv3.0")
	assert os.path.exists(args.input), 'Error! File does not exist: %s. Is the path correct?' % args.input
	openFastx()
main()

	

