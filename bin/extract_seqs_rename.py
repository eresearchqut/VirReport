#!/usr/bin/env python

"""
Rename header of sequence and format fasta file so there is no line wrapping
and exactly 2 lines per record

"""

import sys, os.path

FROM = 0

def __print_contig__(name, contig):
    print(name.strip())
    print(contig.upper())

def scaff_split(lines, cutoff, FROM):
    current_name = lines[0].lstrip('>')
    contig = ""
    n = FROM
    for line in lines:
        if line[0] == ">":
            if len(contig) > cutoff:
                n = n + 1
                __print_contig__(">CONTIG"+str(n).zfill(6), contig) #+"|"+current_name,contig)
                current_name = line.lstrip('>')
            contig = ""
        else:
            contig += line.strip()
    if len(contig) >= cutoff:
        n = n + 1
        __print_contig__(">CONTIG"+str(n).zfill(6), contig) #+"|"+current_name,contig)

if __name__ == "__main__":
    import sys, os.path
    if len(sys.argv) >= 3 and os.path.exists(sys.argv[1]):
        FILE = open(sys.argv[1],"r")
        cutoff = int(sys.argv[2])
        lines = FILE.readlines()
        if len(sys.argv) == 4:
            FROM = int(sys.argv[3])
        scaff_split(lines, cutoff, FROM)
    else:
        print("i need to be given input!")
        print("usage: extract_seqs.py <input.fa> MIN_LENGTH")
