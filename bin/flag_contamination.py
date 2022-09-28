#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
import glob
import time


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load VSD pipeline results")
    # All the required arguments #
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--read_size", type=str)
    parser.add_argument("--method", type=str)
    parser.add_argument("--viral_db", type=str)
    parser.add_argument("--dedup", type=str)
    parser.add_argument("--diagno", type=str)
    parser.add_argument("--sampleinfopath", type=str)
    parser.add_argument("--targets", type=str)

    args = parser.parse_args()
    threshold = args.threshold
    readsize = args.read_size
    method = args.method
    viral_db = args.viral_db
    dedup = args.dedup
    diagno = args.diagno
    sampleinfo = args.sampleinfopath
    targets = args.targets

    timestr = time.strftime("%Y%m%d-%H%M%S")

    run_data = pd.DataFrame()
    for fl in glob.glob("*_top_scoring_targets_with_cov_stats*.txt"):
        sample_data = pd.read_csv(fl, header=0, sep="\t",index_col=None)
        run_data = run_data.append(sample_data)
    print (run_data)
    run_data["read size"] = readsize
    
    if viral_db == "true":
        if dedup == "true":
            run_data = run_data[["Sample","Species","sacc","naccs","length","slen","cov","av-pident","stitle", "qseqids", "contig_ind_lengths", "cumulative_contig_len", "contig_lenth_min", "contig_lenth_max", "ICTV_information", "Mean read depth","Read count","Dedup read count","Dup %","FPKM","PCT_5X","PCT_10X",]]
        else:
            run_data = run_data[["Sample","Species","sacc","naccs","length","slen","cov","av-pident","stitle", "qseqids", "contig_ind_lengths", "cumulative_contig_len", "contig_lenth_min", "contig_lenth_max", "ICTV_information", "Mean read depth","Read count","RPM","FPKM","PCT_5X","PCT_10X",]]
        
        run_data["FPKM"] = run_data["FPKM"].astype(float)
        run_data["FPKM_max"] = run_data.groupby(["Species"])["FPKM"].transform(max)
        run_data["threshold_value"]=run_data["FPKM_max"]*threshold
        run_data["contamination_flag"] = np.where(run_data["FPKM"] <= run_data["threshold_value"], True, False)
        run_data["contamination_flag"] = np.where(run_data["FPKM_max"] <= 10, "NA", run_data["contamination_flag"])
        run_data = run_data.sort_values(["Sample", "stitle"], ascending = (True, True))
        

        if diagno == "true":
            run_data['Evidence_category'] = np.where((run_data['av-pident'] >= 90) & (run_data['PCT_10X'] >= 0.7) & (run_data['length'] >= 45), "KNOWN",
                                        np.where((run_data['av-pident'] >= 90) & (run_data['PCT_10X'] >= 0.1) & (run_data['length'] >= 45), "KNOWN_FRAGMENT",
                                        np.where((run_data['av-pident'] < 90) & (run_data['av-pident'] > 70)  & (run_data['length'] >= 45) &  (run_data['PCT_10X'] >= 0.1), "PUTATIVE_NOVEL","EXCLUDE")))
            run_data = run_data.sort_values(["Sample", "Species"], ascending = (True, True))
            run_data.drop_duplicates()

            grouped_summary=run_data[['Sample', 'Species']]
            grouped_summary = grouped_summary.groupby('Sample', as_index=False).agg(','.join)
            grouped_summary["Species"] = grouped_summary["Species"].str.replace(",",", ")
            print(grouped_summary)

            if sampleinfo is not None:
                sampleinfo_data = pd.read_csv(sampleinfo, header=0, sep="\t",index_col=None)
                run_data = pd.merge(sampleinfo_data, run_data, on="Sample", how='outer').fillna('NA')
                grouped_summary = pd.merge(sampleinfo_data, grouped_summary, on="Sample", how='outer').fillna('NA')
            
            run_data.to_csv("VirReport_detection_summary_" +  "_" + readsize + "_viral_db_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")
            grouped_summary.to_csv("VirReport_detection_summary_collapsed" + readsize + "_viral_db_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")  
        
        else:
            run_data.to_csv("VirReport_detection_summary_" + readsize + "_viral_db_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")
    
    else:
        if dedup == "true":
            run_data = run_data[["Sample","Species","sacc","naccs","length","slen","cov","av-pident","stitle","qseqids","contig_ind_lengths","cumulative_contig_len","contig_lenth_min","contig_lenth_max","Mean read depth","Read count","Dedup read count","Dup %","FPKM","PCT_5X","PCT_10X"]]
        else:
            run_data = run_data[["Sample","Species","sacc","naccs","length","slen","cov","av-pident","stitle","qseqids","contig_ind_lengths","cumulative_contig_len","contig_lenth_min","contig_lenth_max","Mean read depth","Read count","FPKM","PCT_5X","PCT_10X"]]
        #testing both FPKM and RPM
        if method == "FPKM":
            run_data["FPKM_max"] = run_data.groupby(["Species"])["FPKM"].transform(max)
            run_data["threshold_value"]=run_data["FPKM_max"]*threshold
            run_data["contamination_flag"] = np.where(run_data["FPKM"] <= run_data["threshold_value"], True, False)
            run_data["contamination_flag"] = np.where(run_data["FPKM_max"] <= 10, "NA", run_data["contamination_flag"])
        
        elif method == "read_counts_normalised":
            run_data["FPKM_max"] = run_data.groupby(["Species"])["read_counts_normalised"].transform(max)
            run_data["threshold_value"]=run_data["FPKM_max"]*threshold
            run_data["contamination_flag"] = np.where(run_data["read_counts_normalised"] <= run_data["threshold_value"], True, False)
            run_data["contamination_flag"] = np.where(run_data["FPKM_max"] <= 10, "NA", run_data["contamination_flag"])
        
        run_data = run_data.sort_values(["Sample", "Species"], ascending = (True, True))

        if diagno == "true":
            run_data['Evidence_category'] = np.where((run_data['av-pident'] >= 90) & (run_data['PCT_10X'] >= 0.7) & (run_data['length'] >= 45), "KNOWN",
                                        np.where((run_data['av-pident'] >= 90) & (run_data['PCT_10X'] >= 0.1) & (run_data['length'] >= 45), "KNOWN_FRAGMENT",
                                        np.where((run_data['av-pident'] < 90) & (run_data['av-pident'] > 70)  & (run_data['length'] >= 45), "PUTATIVE_NOVEL","EXCLUDE")))
            
            
            targets_df = pd.read_csv(targets, header=0, sep="\t", index_col=None)
            targets_df["Species"] = targets_df["Species"].astype(str)
            #targets_df["Species"] = targets_df["Species"].str.lower()
            targets_list = targets_df["Species"].tolist()
            run_data["SSG_category"] = run_data['Species'].isin(targets_list)
            run_data['SSG_category'] = run_data['SSG_category'].map({True: 'Quarantinable', False: 'Higher_plant_viruses'})
            
            run_data = run_data.sort_values(["Sample", "Species"], ascending = (True, True))
            run_data.drop_duplicates()

            grouped_summary=run_data[['Sample', 'Species']]
            grouped_summary = grouped_summary.groupby('Sample', as_index=False).agg(','.join)
            grouped_summary["Species"] = grouped_summary["Species"].str.replace(",",", ")
            print(grouped_summary)

            if sampleinfo is not None:
                sampleinfo_data = pd.read_csv(sampleinfo, header=0, sep="\t",index_col=None)
                run_data = pd.merge(sampleinfo_data, run_data, on="Sample", how='outer').fillna('NA')
                grouped_summary = pd.merge(sampleinfo_data, grouped_summary, on="Sample", how='outer').fillna('NA')
            
            run_data.to_csv("VirReport_detection_summary_"  + readsize + "_ncbi_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")
            grouped_summary.to_csv("VirReport_detection_summary_"  + readsize + "_ncbi_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")
            
        else:
            run_data.to_csv("VirReport_detection_summary_" + readsize + "_ncbi_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")

if __name__ == "__main__":
    main()