#!/usr/bin/env python
import argparse
import pandas as pd
import collections
from collections import OrderedDict
from operator import itemgetter

def main():
    parser = argparse.ArgumentParser(description="Load blast results")
    parser.add_argument("--virus_list", type=str)
    parser.add_argument("--contig_fasta", type=str)
    parser.add_argument("--out", type=str)
    args = parser.parse_args()
    viruslist = args.virus_list
    contigs_fasta = args.contig_fasta
    out = args.out


    with open(contigs_fasta, 'r') as file:
        counter = collections.Counter()
        header = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
                continue
            counter[header] += len(line)

    #print(counter)
    raw_data = pd.read_csv(viruslist, header=0, sep="\t",index_col=None)
    #print(raw_data)
    full_length_list = []
    full_unique_list = []
    sum_list = []
    min_list = []
    max_list = []
    contig_count_list = []
    longest_contig_list = []
    
    for index, row in raw_data.iterrows():
        #retrieve unique list of contig name for a given row
        unique_list = list(set(row['qseqids'].split(",")))
        contig_count = len(unique_list)
        length_list = []
        
        #extract length of each contig
        contig_len_dic = {}
        for i in unique_list:
            ind_length = (list(counter.values())[list(counter.keys()).index(i)])
            length_list.append(ind_length)
            contig_len_dic[i] = ind_length
        contig_string = ', '.join(map(str,unique_list))

        sort_orders = sorted(contig_len_dic.items(), key=lambda x: x[1])
        sorted_dict = OrderedDict(sorted(contig_len_dic.items(), key = itemgetter(1)))
        print(sorted_dict)
        longest_contig = list(sorted_dict.keys())[-1]

        contig_seq = ""
        string = ">" + longest_contig
        with open(contigs_fasta, 'r') as f:
            for line in f:
                if string in line:
                    contig_seq += line.strip()
                    contig_seq += ' '
                    contig_seq += next(f).strip()

            contig_seq = contig_seq.replace('"', '')
        longest_contig_list.append(contig_seq)

        sum_numbers = sum(length_list)
        max_length = max(length_list)
        min_length = min(length_list)
        sum_list.append(sum_numbers)
        full_length_list.append(sort_orders)
        full_unique_list.append(contig_string)
        max_list.append(max_length)
        min_list.append(min_length)
        contig_count_list.append(contig_count)

    raw_data['unique_contig_list'] = pd.Series(full_unique_list)
    raw_data['contig_ind_lengths'] = pd.Series(full_length_list)

    raw_data['cumulative_contig_len'] = pd.Series(sum_list)
    raw_data['contig_lenth_min'] = pd.Series(min_list)
    raw_data['contig_lenth_max'] = pd.Series(max_list)
    raw_data['contig_count'] = pd.Series(contig_count_list)
    raw_data['longest_contig_fasta'] = pd.Series(longest_contig_list)#, dtype=str)

    raw_data = raw_data.drop(["qseqids"], axis=1)
    raw_data = raw_data.rename(columns={"unique_contig_list": "qseqids"})
    raw_data = raw_data.drop(["naccs"], axis=1)
    raw_data = raw_data.rename(columns={"contig_count": "naccs"})

    raw_data = raw_data[['sacc', 'naccs', 'length', 'slen', 'cov', 'av-pident', 'stitle', 'qseqids', 'contig_ind_lengths', 'cumulative_contig_len', 'contig_lenth_min', 'contig_lenth_max', 'longest_contig_fasta']]
    dataTypeSeries = raw_data.dtypes
    print(dataTypeSeries)

    print(raw_data)
    raw_data.to_csv(out, index=None, sep="\t",float_format="%.2f")  

if __name__ == "__main__":
    main()