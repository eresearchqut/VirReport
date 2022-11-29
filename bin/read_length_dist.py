#!/usr/bin/env python

"""
This script will output read length distribution graph,
text file with counts table of unique read length and
fasta file of reads of wanted length specified with minlen and maxlen
# read_length_dist.py --input seqs.fasta --minlen 21 --maxlen 22
"""

###############################################################################
# Modules #
import argparse, sys, time, getpass, locale
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

################################################################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("--input", help="The fasta file to process", type=str)
args        = parser.parse_args()
input_path  = args.input
lengths = list(map(len, SeqIO.parse(input_path, 'fasta')))
sys.stderr.write("Read all lengths (%i sequences)\n" % len(lengths))
sys.stderr.write("Longest sequence: %i bp\n" % max(lengths))
sys.stderr.write("Shortest sequence: %i bp\n" % min(lengths))


sys.stderr.write("Deriving read length distribution\n")
sample = (args.input).replace(".rename.fa", "")
print(sample)
#build array to store all read lengths
arr = np.array(lengths)
#return counts for read lengths
u, c = np.unique(np.array(lengths), return_counts=True)
with open('%s_read_length_dist.txt' % sample, 'w') as f:
    np.savetxt(f, np.stack([u, c]).T,delimiter='\t', fmt='%12s')

sys.stderr.write("Making graph...\n")

fig = plt.figure(figsize=(20, 5))
labels, counts = np.unique(arr, return_counts=True)
plt.bar(labels, counts, color='green', align='center', width=0.5)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-6,9))
plt.gca().yaxis.set_major_formatter(formatter)
plt.gca().set_xticks(labels)
plt.xticks(rotation='vertical')
plt.title(sample)
plt.xlabel('read length (bp)')
plt.ylabel('# of reads')
plt.savefig('%s_read_length_dist.pdf' % sample, format='pdf')
plt.savefig('%s_read_length_dist.png' % sample, format='png')