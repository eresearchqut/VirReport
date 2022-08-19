#!/usr/bin/env python
import argparse
import pandas as pd
import collections

# Using a stringio to emulate a file
def main():
    parser = argparse.ArgumentParser(description="Load blast results")
    # All the required arguments #
    parser.add_argument("--virus_list", type=str)
    parser.add_argument("--contig_fasta", type=str)
    parser.add_argument("--sample_name", type=str)
    parser.add_argument("--read_size", type=str)
    parser.add_argument("--out", type=str)
    args = parser.parse_args()

    viruslist = args.virus_list
    contigs_fasta = args.contig_fasta
    sample = args.sample_name
    read_size = args.read_size
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

    raw_data = raw_data.drop(["qseqids"], axis=1)
    raw_data = raw_data.rename(columns={"unique_contig_list": "qseqids"})
    raw_data = raw_data.drop(["naccs"], axis=1)
    raw_data = raw_data.rename(columns={"contig_count": "naccs"})

    raw_data = raw_data[['sacc', 'naccs', 'length', 'slen', 'cov', 'av-pident', 'stitle', 'qseqids', 'contig_ind_lengths', 'cumulative_contig_len', 'contig_lenth_min', 'contig_lenth_max']]
    dataTypeSeries = raw_data.dtypes
    print(dataTypeSeries)

    print(raw_data)
    #raw_data.to_csv("summary_" + sample + "_cap3_" + read_size + "_blastn_vs_NT_top5Hits_virus_viroids_final.txt", index=None, sep="\t",float_format="%.2f")
    raw_data.to_csv(out, index=None, sep="\t",float_format="%.2f")  

if __name__ == "__main__":
    main()