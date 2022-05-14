#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from functools import reduce
import glob

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load VSD pipeline results")
    # All the required arguments #
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--read_size", type=str)
    parser.add_argument("--method", type=str)

    args = parser.parse_args()
    threshold = args.threshold
    readsize = args.read_size
    method = args.method

    run_data = pd.DataFrame()
    for fl in glob.glob("*_top_scoring_targets_with_cov_stats.txt"):
        sample_data = pd.read_csv(fl, header=0, sep="\t",index_col=None)
        run_data = run_data.append(sample_data)

    #run_data = run_data[["Sample","sacc","naccs","length","slen","cov","av-pident","qseqids","Targetted_sp_generic_name","Mean coverage","Read count","Dedup read count","Dup %","RPM","FPKM","PCT_1X","PCT_5X","PCT_10X","PCT_20X"]]
    run_data = run_data[["Sample","sacc","naccs","length","slen","cov","av-pident","qseqids","Targetted_sp_generic_name","Mean coverage","Read count","RPM","FPKM","PCT_1X","PCT_5X","PCT_10X","PCT_20X"]]
    run_data["read size"] = readsize
    if method == "FPKM":
        run_data["count_max"] = run_data.groupby(["Targetted_sp_generic_name"])["FPKM"].transform(max)
        run_data["threshold_value"]=run_data["count_max"]*threshold
        run_data["contamination_flag"] = np.where(run_data["FPKM"] <= run_data["threshold_value"], True, False)
        run_data["contamination_flag"] = np.where(run_data["count_max"] <= 10, "NA", run_data["contamination_flag"])
        
    elif method == "read_counts_normalised":
        run_data["count_max"] = run_data.groupby(["Targetted_sp_generic_name"])["read_counts_normalised"].transform(max)
        run_data["threshold_value"]=run_data["count_max"]*threshold
        run_data["contamination_flag"] = np.where(run_data["read_counts_normalised"] <= run_data["threshold_value"], True, False)
        run_data["contamination_flag"] = np.where(run_data["count_max"] <= 10, "NA", run_data["contamination_flag"])
    
    run_data = run_data.sort_values(["Sample", "Targetted_sp_generic_name"], ascending = (True, True))
    run_data.to_csv("run_top_scoring_targets_with_cov_stats_with_cont_flag" +  "_" + str(method) + "_" + str(threshold) + '_'  + readsize + ".txt", index=None, sep="\t",float_format="%.2f")
if __name__ == "__main__":
    main()
