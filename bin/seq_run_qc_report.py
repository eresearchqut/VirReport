#!/usr/bin/env python
import pandas as pd
from functools import reduce
import glob
import re
import csv
import os
import numpy as np

def main():

    raw_read_counts_dict = {}
    for umitools_out in glob.glob("*_umi_tools.log"):
        raw_reads = ()
        umi_cleaned_reads = ()
        line_number = 0
        sample = (os.path.basename(umitools_out).replace('_umi_tools.log', ''))
        with open(umitools_out, 'r') as f:
            for line in f:
                string_to_search1 = ("Input Reads")
                string_to_search2 = ("Reads output")
                line_number += 1
                if string_to_search1 in line:
                    elements = line.split("Input Reads: ")
                    raw_reads = int(elements[1].strip())
                elif string_to_search2 in line:
                    elements = line.split("Reads output: ")
                    umi_cleaned_reads = int(elements[1].strip())
        raw_read_counts_dict[sample] = [raw_reads, umi_cleaned_reads]
        f.close()
    
    for cutadapt_qual_filt_out in glob.glob("*_qual_filtering_cutadapt.log"):
        qfiltered_reads = ()
        line_number = 0
        sample = (os.path.basename(cutadapt_qual_filt_out).replace('_qual_filtering_cutadapt.log', ''))
        with open(cutadapt_qual_filt_out, 'r') as f:
            for line in f:
                string_to_search1 = ("Reads written (passing filters):")
                line_number += 1
                if string_to_search1 in line:
                    qfiltered_reads = line.split("Reads written (passing filters): ")[1]
                    qfiltered_reads = re.sub(r' \(.*\)', '', qfiltered_reads).strip()
                    qfiltered_reads = int(re.sub(r',', '', qfiltered_reads).strip())
        raw_read_counts_dict[sample].append(qfiltered_reads)
        f.close()

    for fastp_out in glob.glob("*_fastp.json"):
        total_filtered_bases = ()
        q20_bases = ()
        q30_bases = ()
        gc_content = ()

        line_number = 0
        sample = (os.path.basename(fastp_out).replace('_fastp.json', ''))
        with open(fastp_out, 'r') as f:
            for line in f:
                string_to_search1 = ("total_bases")
                string_to_search2 = ("gc_content")
                line_number += 1
                if string_to_search1 in line:
                    total_filtered_bases = line.split(":")[1]
                    total_filtered_bases = int(re.sub(r',', '', total_filtered_bases).strip())
                    q20_bases = f.readline().split(":")[1]
                    q20_bases = int(re.sub(r',', '', q20_bases).strip())
                    q30_bases = f.readline().split(":")[1]
                    q30_bases = int(re.sub(r',', '', q30_bases).strip())
                elif string_to_search2 in line:
                    gc_content = float(line.split(":")[1].strip())
        f.close()
        raw_read_counts_dict[sample].append(total_filtered_bases)
        raw_read_counts_dict[sample].append(q20_bases)
        raw_read_counts_dict[sample].append(q30_bases)
        raw_read_counts_dict[sample].append(gc_content)


    for bowtie_blacklist_out in glob.glob("*_blacklist_filter.log"):
        usable_source_reads = ()
        line_number = 0
        sample = (os.path.basename(bowtie_blacklist_out).replace('_blacklist_filter.log', ''))
        with open(bowtie_blacklist_out, 'r') as f:
            for line in f:
                string_to_search1 = ("reads that failed to align")
                line_number += 1
                if string_to_search1 in line:
                    elements = line.split(": ")
                    usable_source_reads = elements[1].strip()
                    usable_source_reads = int(re.sub(r' \(.*\)', '', usable_source_reads).strip())
        f.close()
        raw_read_counts_dict[sample].append(usable_source_reads)

    for cutadapt_18_25_out in glob.glob("*_18-25nt_cutadapt.log"):
        usable_reads_18_25 = ()
        line_number = 0
        sample = (os.path.basename(cutadapt_18_25_out).replace('_18-25nt_cutadapt.log', ''))
        with open(cutadapt_18_25_out, 'r') as f:
            for line in f:
                string_to_search1 = ("Reads written (passing filters):")
                line_number += 1
                if string_to_search1 in line:
                    usable_reads_18_25 = line.split("Reads written (passing filters): ")[1]
                    usable_reads_18_25 = re.sub(r' \(.*\)', '', usable_reads_18_25).strip()
                    usable_reads_18_25 = int(re.sub(r',', '', usable_reads_18_25).strip())
        f.close()
        raw_read_counts_dict[sample].append(usable_reads_18_25)

    for cutadapt_21_22_out in glob.glob("*_21-22nt_cutadapt.log"):
        usable_reads_21_22 = ()
        line_number = 0
        sample = (os.path.basename(cutadapt_21_22_out).replace('_21-22nt_cutadapt.log', ''))
        with open(cutadapt_21_22_out, 'r') as f:
            for line in f:
                string_to_search1 = ("Reads written (passing filters):")
                line_number += 1
                if string_to_search1 in line:
                    usable_reads_21_22 = line.split("Reads written (passing filters): ")[1]
                    usable_reads_21_22 = re.sub(r' \(.*\)', '', usable_reads_21_22).strip()
                    usable_reads_21_22 = int(re.sub(r',', '', usable_reads_21_22).strip())
        f.close()
        raw_read_counts_dict[sample].append(usable_reads_21_22)

    for cutadapt_24_out in glob.glob("*_24nt_cutadapt.log"):
        usable_reads_24 = ()
        line_number = 0
        sample = (os.path.basename(cutadapt_24_out).replace('_24nt_cutadapt.log', ''))
        with open(cutadapt_24_out, 'r') as f:
            for line in f:
                string_to_search1 = ("Reads written (passing filters):")
                line_number += 1
                if string_to_search1 in line:
                    usable_reads_24 = line.split("Reads written (passing filters): ")[1]
                    usable_reads_24 = re.sub(r' \(.*\)', '', usable_reads_24).strip()
                    usable_reads_24 = int(re.sub(r',', '', usable_reads_24).strip())
        f.close()
        raw_read_counts_dict[sample].append(usable_reads_24) 
    
    run_data_df = pd.DataFrame([([k] + v) for k, v in raw_read_counts_dict.items()], columns=['sample','raw_reads','umi_cleaned_reads', 'quality filtered reads > 18bp', 'total filtered bases', 'q20 bases', 'q30 bases', 'gc content', 'usable reads type', 'usable reads 18-25 nt', 'usable reads 21-22nt', 'usable reads 24nt'])
    print(run_data_df.dtypes)

    run_data_df.set_index('sample')
    #For all columns in the dataframe that are of dtype int64, add commas
    run_data_df.update(run_data_df.select_dtypes(include=['int64']).applymap('{:,}'.format))
    #Retain 2 decimal point format for GC content column
    run_data_df.loc[:, 'gc content'] = run_data_df['gc content'].map('{:.2f}'.format)
    run_data_df = run_data_df.sort_values("sample")
    run_data_df.to_csv('run_qc_report.txt', index = None, sep="\t")

if __name__ == '__main__':
    main()
