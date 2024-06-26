#!/usr/bin/env python
import pandas as pd
from functools import reduce
import glob
import re
import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import time


def main():
    timestr = time.strftime("%Y%m%d-%H%M%S")
    read_origin_dict = {}
    for umitools_out in glob.glob("*bowtie.log"):
        sample = ()
        total_reads = ()
        rRNA = ()
        mt_pt_other = ()
        miRNA = ()
        plant_tRNA = ()
        plant_nc = ()
        artefacts = ()
        viral = ()
        leftover = ()

        #sample = umitools_out.replace('_bowtie.log', '')
        with open(umitools_out, 'r') as f:
            sample = f.readline().strip('\n')
            print(sample)
            for line in f:
                
                if ("rRNA alignment:") in line:
                    line = next(f)
                    elements = line.split("# reads processed: ")
                    total_reads = int(elements[1].strip())
                    line = next(f)
                    elements2 = line.split("# reads with at least one alignment: ")
                    rRNA = elements2[1].strip()
                    rRNA = int(re.sub(r' \(.*\)', '', rRNA).strip())

                elif ("miRNA alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    miRNA = elements[1].strip()
                    miRNA = int(re.sub(r' \(.*\)', '', miRNA).strip())

                elif ("plant_tRNA alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_tRNA = elements[1].strip()
                    plant_tRNA = int(re.sub(r' \(.*\)', '', plant_tRNA).strip())

                elif ("plant_pt_mt_other_genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    mt_pt_other = elements[1].strip()
                    mt_pt_other = int(re.sub(r' \(.*\)', '', mt_pt_other).strip())
        
                elif ("miRNA alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    miRNA = elements[1].strip()
                    miRNA = int(re.sub(r' \(.*\)', '', miRNA).strip())
                
                elif ("plant_noncoding alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_nc = elements[1].strip()
                    plant_nc = int(re.sub(r' \(.*\)', '', plant_nc).strip())
                
                
                elif ("artefacts alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    artefacts = elements[1].strip()
                    artefacts = int(re.sub(r' \(.*\)', '', artefacts).strip())
                
                elif ("plant_virus_viroid alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    viral = elements[1].strip()
                    viral = int(re.sub(r' \(.*\)', '', viral).strip())
                    line = next(f)
                    elements2 = line.split("# reads that failed to align: ")
                    leftover = elements2[1].strip()
                    leftover = int(re.sub(r' \(.*\)', '', leftover).strip())
                    
        f.close()
        read_origin_dict[sample] = [rRNA, plant_tRNA, mt_pt_other, plant_nc, artefacts,\
                                    miRNA, viral, leftover, total_reads]
        #sort dictionary by key (ie sample name)
        read_origin_dict = collections.OrderedDict(sorted(read_origin_dict.items()))
    
    #convert dictionary into pandas df
    read_origin_df = pd.DataFrame([([k] + v) for k, v in read_origin_dict.items()], columns=['sample','rRNA_total', 
                                'plant_tRNA', 'mitochondrial_and_plastid_genes', 'plant_non_coding_RNA',
                                'PhiX_and_artefacts', 'miRNA', 'virus_and_viroids', 'leftover', 'total_reads'])              
    
    read_origin_df = read_origin_df.set_index(read_origin_df.columns[0])
    read_origin_df = read_origin_df.sort_index(ascending=True)
    read_origin_df.to_csv('read_origin_counts.' + timestr + '.txt', sep="\t", float_format="%.2f")
    
    read_origin_df = read_origin_df.iloc[:, :-1]
    
    pc_df = read_origin_df.apply(lambda x: 100 * x / float(x.sum()), axis=1)

    pc_df.to_csv('read_origin_detailed_pc.' + timestr + '.txt', sep="\t", float_format="%.2f")

    pc_df.plot.barh(stacked=True, color=['#000000', '#C5C9C7', '#808080', 'purple', 'yellow', '#069AF3', '#15B01A', '#E6E6FA'], figsize=(8,15)).legend(loc='lower center',bbox_to_anchor=(0.5, -0.3))
    plt.tight_layout()
    plt.savefig('read_RNA_source.' + timestr + '.png', format="png")
    plt.savefig('read_RNA_source.' + timestr + '.pdf', format="pdf")
    plt.close()

    column_names = ['rRNA_total',  'plant_tRNA']
    pc_df['rRNA_and_tRNA']= pc_df[column_names].sum(axis=1)
    column_names = ['miRNA',  'virus_and_viroids']
    pc_df['miRNA/vsiRNA']= pc_df[column_names].sum(axis=1)
    pc_df['rRNA/tRNA_flag'] = pc_df['rRNA_and_tRNA'].apply(lambda x: 'High % of rRNA/tRNA' if x >= 50 else '')
    pc_df['miRNA/vsiRNA_flag'] = pc_df['miRNA/vsiRNA'].apply(lambda x: 'Low % of miRNA/vsiRNA' if x <= 10 else '')
    print(pc_df)
    pc_df.to_csv('read_origin_pc_summary.' + timestr + '.txt', sep="\t", float_format="%.2f")

if __name__ == '__main__':
    main()
