#!/usr/bin/env python3

'''Usage: python gcLenCorPlots.py -i <input fasta or fastq> [-m <method> -r <range for heatmap> -t <trim maximum length> -s <normalize to number> -ec <error bar color>]'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
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
	if args.trim is not None: #optional length trimming and shuffle options
		df = df.drop(df[df.length > int(args.trim)].index).reset_index()
		print("Length trimmed to maximum %i" % int(args.trim))
		print("Number of trimmed reads: %i" % len(df.gcContent))
	if args.shuffle is not None:
		df = df.sample(int(args.shuffle))
		print("Number of reads normalized to %i" % int(args.shuffle))
	printStats(df)

def printStats(df): #generate statistics, print to screen, print out the statistics file
	print("Mean GC content: %.2f" % np.mean(df['gcContent']))
	print("Median GC content: %.2f" % np.median(df['gcContent']))
	print("Mean fragment length: %i" % np.mean(df['length']))
	print("Median fragment length: %i" % np.median(df['length']))
	print("Fragment length range: %i : %i" % (min(df['length']), max(df['length'])))
	print("GC content range: %.2f : %.2f" % (min(df['gcContent']), max(df['gcContent'])))	
	dfGrouped = df.groupby(by='length').agg(['count', 'mean', 'median', 'std']).reset_index()
	dfGrouped['perc'] = dfGrouped['gcContent', 'count']/dfGrouped['gcContent', 'count'].sum()
	with open('%s_data_out.txt' % args.input, 'w') as outfile:
		dfGrouped.to_csv(outfile, sep="\t", index=False)
	plotDat(dfGrouped)

def plotDat(dfGrouped):
	plt.xlim(15, 90) #change for new x axis range
	plt.ylim(20, 200) #change this for new y axis range
	plt.suptitle(args.input + "\n" + "n= " + str(dfGrouped['gcContent', 'count'].sum()))
	plt.xlabel('GC content (%)')
	plt.ylabel('Read length (bp)')
	cm = plt.cm.get_cmap('YlOrRd') #change this for color scheme other than yellow to red - for more options google matplotlib color maps
	plt.errorbar(dfGrouped['gcContent', args.method], dfGrouped['length', ''], xerr=dfGrouped['gcContent', 'std'], linestyle="None", marker="None", color=args.errorbarColor)
	if args.range == 'num': #color range options
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=100, vmax=2000, marker='o', edgecolors='None', s=25, zorder=2) #change vmin and vmax to alter range of absolute number of reads to custom number
	elif args.range == 'perc':
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['perc']), cmap=cm, vmin=0.0, vmax=1.0, marker='o', edgecolors='None', s=25, zorder=2)
	elif args.range == 'max':
		plt.scatter(dfGrouped['gcContent', args.method], dfGrouped['length', ''], c=list(dfGrouped['gcContent', 'count']), cmap=cm, vmin=min(dfGrouped['gcContent', 'count']), vmax=max(dfGrouped['gcContent', 'count']), marker='o', edgecolors='None', s=25, zorder=2)
	plt.colorbar()
	plt.draw()
	plt.savefig('%s_plot.pdf' % args.input)

def main():
	print("Author: Allison E. Mann")
	print("Cite: Mann AE et al. 2018. Differential preservation of endogenous human and microbial DNA in dental calculus and dentin. Scientific Reports 8:9822")
	print("Author: Allison E. Mann (allison.mann@unthsc.edu)")
	print("Copyright 2018: GPLv3.0\n")
	assert os.path.exists(args.input), 'Error! File does not exist: %s. Is the path correct?' % args.input
	openFastx()

main()

__author__ = "Allison E. Mann"
__license__ = "GPL"
__version__ = "1.0.1"
__email__="allison.e.mann@gmail.com"

	

