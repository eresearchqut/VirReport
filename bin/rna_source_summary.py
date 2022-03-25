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


def main():

    read_origin_dict = {}
    for umitools_out in glob.glob("*_bowtie.log"):
        sample = ()
        total_reads = ()
        rRNA_5S = ()
        nc_SSU_and_LSU_rRNA = ()
        mt_rRNA = ()
        pt_rRNA = ()
        mt_other = ()
        pt_other = ()
        plant_miRNA = ()
        other_miRNA = ()
        plant_tRNA = ()
        plant_nc = ()
        plant_transposons = ()
        phix = ()
        vector = ()
        viral = ()
        leftover = ()

        sample = umitools_out.replace('_bowtie.log', '')
        #print(sample_file)
        with open(umitools_out, 'r') as f:
            for line in f:
                #sample = line.rstrip()
                #sanity check
                #if sample == sample_file:
                #    print ("sample name matches")
                if ("5S rRNA genes alignment:") in line:
                    line = next(f)
                    elements = line.split("# reads processed: ")
                    total_reads = int(elements[1].strip())
                    line = next(f)
                    elements2 = line.split("# reads with at least one alignment: ")
                    rRNA_5S = elements2[1].strip()
                    rRNA_5S = int(re.sub(r' \(.*\)', '', rRNA_5S).strip())
                elif ("nc SSU and LSU rRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    nc_SSU_and_LSU_rRNA = elements[1].strip()
                    nc_SSU_and_LSU_rRNA = int(re.sub(r' \(.*\)', '', nc_SSU_and_LSU_rRNA).strip())
                elif ("mt rRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    mt_rRNA = elements[1].strip()
                    mt_rRNA = int(re.sub(r' \(.*\)', '', mt_rRNA).strip())
                elif ("pt rRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    pt_rRNA = elements[1].strip()
                    pt_rRNA = int(re.sub(r' \(.*\)', '', pt_rRNA).strip())
                elif ("mt other genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    mt_other = elements[1].strip()
                    mt_other = int(re.sub(r' \(.*\)', '', mt_other).strip())
                elif ("pt other genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    pt_other = elements[1].strip()
                    pt_other = int(re.sub(r' \(.*\)', '', pt_other).strip())
                elif ("plant miRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_miRNA = elements[1].strip()
                    plant_miRNA = int(re.sub(r' \(.*\)', '', plant_miRNA).strip())
                elif ("other miRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    other_miRNA = elements[1].strip()
                    other_miRNA = int(re.sub(r' \(.*\)', '', other_miRNA).strip())
                elif ("tRNA genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_tRNA = elements[1].strip()
                    plant_tRNA = int(re.sub(r' \(.*\)', '', plant_tRNA).strip())
                elif ("plant nc genes alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_nc = elements[1].strip()
                    plant_nc = int(re.sub(r' \(.*\)', '', plant_nc).strip())
                elif ("plant transposons alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    plant_transposons = elements[1].strip()
                    plant_transposons = int(re.sub(r' \(.*\)', '', plant_transposons).strip())
                elif ("PhiX alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    phix = elements[1].strip()
                    phix = int(re.sub(r' \(.*\)', '', phix).strip())
                elif ("Vector alignment:") in line:
                    line = next(f)
                    line = next(f)
                    elements = line.split("# reads with at least one alignment: ")
                    vector = elements[1].strip()
                    vector = int(re.sub(r' \(.*\)', '', vector).strip())
                elif ("Virus and viroid alignment:") in line:
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
        read_origin_dict[sample] = [rRNA_5S, nc_SSU_and_LSU_rRNA, mt_rRNA, pt_rRNA, mt_other,\
                                    pt_other, plant_tRNA, plant_nc, plant_transposons, phix, vector,\
                                    plant_miRNA, other_miRNA, viral, leftover, total_reads]
        #sort dictionary by key (ie sample name)
        read_origin_dict = collections.OrderedDict(sorted(read_origin_dict.items()))
    
    #convert dictionary into pandas df
    read_origin_df = pd.DataFrame([([k] + v) for k, v in read_origin_dict.items()], columns=['sample','5S_rRNA', 
                                'nuclear_SSU_and_LSU_rRNA', 'mitochondrial_rRNA', 'plastid_rRNA', 'mitochondrial_other', 
                                'plastid_other', 'plant_tRNA', 'plant_non_coding_RNA', 'plant_transposons', 'PhiX', 
                                'vector_and_artefacts', 'plant_miRNA', 'other_miRNA',  'virus_and_viroids', 'leftover', 'total_reads'])              
    
    read_origin_df = read_origin_df.set_index(read_origin_df.columns[0])
    read_origin_df = read_origin_df.sort_index(ascending=True)
    read_origin_df.to_csv('read_origin_counts.txt', sep="\t", float_format="%.2f")
    
    read_origin_df = read_origin_df.iloc[:, :-1]
    
    pc_df = read_origin_df.apply(lambda x: 100 * x / float(x.sum()), axis=1)

    pc_df.to_csv('read_origin_detailed_pc.txt', sep="\t")

    #print(pc_df)

    #transpose table
    #read_origin_df_T = read_origin_df.transpose()
    #read_origin_df_T.columns = read_origin_df_T.iloc[0] 
    #read_origin_df_T = read_origin_df_T[1:]
    #read_origin_df_T.to_csv('read_origin_counts.txt', sep="\t", float_format="%.2f")

    #derive percent
    #read_origin_df_T = read_origin_df_T[:-1]
    #print(read_origin_df_T)
    #pc_df = read_origin_df_T.apply(lambda x: 100 * x / float(x.sum()))
    #pc_df.to_csv('read_origin_pc.txt', sep="\t")

    #pc_df_T = pc_df.transpose()
    #pc_df_T = pc_df_T.sort_index(ascending=False)
    column_names = ['5S_rRNA', 'nuclear_SSU_and_LSU_rRNA', 'mitochondrial_rRNA', 'plastid_rRNA']
    pc_df['rRNA_total']= pc_df[column_names].sum(axis=1)
    
    
    column_names = ['plant_miRNA', 'other_miRNA']
    pc_df['miRNA_total']= pc_df[column_names].sum(axis=1)
 
    column_names = ['PhiX', 'vector_and_artefacts']
    pc_df['PhiX_and_artefacts']= pc_df[column_names].sum(axis=1)

    column_names = ['mitochondrial_other', 'plastid_other']
    pc_df['mitochondrial_and_plastid_genes']= pc_df[column_names].sum(axis=1)

    pc_df = pc_df.drop(columns=['5S_rRNA', 'nuclear_SSU_and_LSU_rRNA', 'mitochondrial_rRNA', 'plastid_rRNA', 'plant_miRNA', 'other_miRNA', 'PhiX', 'vector_and_artefacts', 'mitochondrial_other', 'plastid_other'])
    pc_df = pc_df[['rRNA_total',  'plant_tRNA', 'mitochondrial_and_plastid_genes', 'plant_non_coding_RNA', 'plant_transposons',
                        'PhiX_and_artefacts', 'miRNA_total',  'virus_and_viroids', 'leftover']]
    

    pc_df.plot.barh(stacked=True, color=['#000000', '#C5C9C7', '#808080', 'purple', 'red', 'yellow',  '#069AF3', '#15B01A', '#E6E6FA'], figsize=(8,15)).legend(loc='lower center',bbox_to_anchor=(0.5, -0.3))
    plt.tight_layout()
    plt.savefig("read_RNA_source.pdf", format='pdf')
    plt.close()

    column_names = ['rRNA_total',  'plant_tRNA']
    pc_df['rRNA_and_tRNA']= pc_df[column_names].sum(axis=1)
    column_names = ['miRNA_total',  'virus_and_viroids']
    pc_df['miRNA/vsiRNA']= pc_df[column_names].sum(axis=1)
    pc_df['rRNA/tRNA_flag'] = pc_df['rRNA_and_tRNA'].apply(lambda x: 'High % of rRNA/tRNA' if x >= 50 else '')
    pc_df['miRNA/vsiRNA_flag'] = pc_df['miRNA_total'].apply(lambda x: 'Low % of miRNA/vsiRNA' if x <= 10 else '')
    print(pc_df)
    pc_df.to_csv('read_origin_pc_summary.txt', sep="\t")

if __name__ == '__main__':
    main()