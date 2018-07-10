# gcLenCorPlots
This script calculates summary statistics for read length and GC content from a fasta or fastq file and generates a statistics table and plot. Written for the dental calculus versus dentin project (citation: Mann AE, Sabin S, Ziesemer KA, Vågene Å, Schroeder H, Ozga A, Sankaranarayanan K, Hofman CA, Fellows-Yates J, Salazar Garcia D, Frohlich B, Aldenderfer M, Hoogland M, Read C, Krause J, Hofman C, Bos K, Warinner C. (2018) Differential preservation of endogenous human and microbial DNA in dental calculus and dentin. Scientific Reports 8:9822.)

Usage: python3 gcLenCorPlots.py -i <fastx> [-m <method> -r <range for heatmap> -t <trim maximum length> -s <normalize to number> -ec <error bar color>]

For help or parameter description: python3 gcLenCorPlots.py -h

Requirements:
[python 3+](https://www.python.org/downloads/)
[scipy](https://www.scipy.org/)
[numpy](http://www.numpy.org/)
[BioPython](https://biopython.org/)
[Matplotlib](https://matplotlib.org/)
[pandas](https://pandas.pydata.org/)
