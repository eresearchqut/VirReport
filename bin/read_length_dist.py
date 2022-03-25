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

# All the required arguments #
parser.add_argument("--input", help="The fasta file to process", type=str)
# parser.add_argument("--minlen", type=int)
# parser.add_argument("--maxlen", type=int)

args        = parser.parse_args()
input_path  = args.input
# minlen      = args.minlen
# maxlen      = args.maxlen

################################################################################
# Read #
#lengths = *map(len, SeqIO.parse(input_path, 'fasta')),
#test
lengths = list(map(len, SeqIO.parse(input_path, 'fasta')))
#lengths = map(len, SeqIO.parse(input_path, 'fasta')),

#print(lengths)

# Report #
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
#test
#print(np.stack([u, c]).T)
#with open('length_distribution.txt','w') as f:
with open('%s_read_length_dist.txt' % sample, 'w') as f:
    np.savetxt(f, np.stack([u, c]).T,delimiter='\t', fmt='%12s')

#extract reads of given length
# sys.stderr.write("Extract reads of given length\n")
# wanted_sequences = []
# for record in SeqIO.parse(input_path, 'fasta'):
#     if len(record.seq) == minlen or len(record.seq) == maxlen:
#         wanted_sequences.append(record)
# print("Found %i sequences of wanted lengths" % len(wanted_sequences))

# wanted_seq_fasta_file = open(sample + "_" + str(minlen) + "-" + str(maxlen) + "nt.fa","w")
# SeqIO.write(wanted_sequences, wanted_seq_fasta_file, "fasta")

sys.stderr.write("Making graph...\n")

fig = plt.figure(figsize=(20, 5))
labels, counts = np.unique(arr, return_counts=True)
#plt.figure(figsize=(20, 3))
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
