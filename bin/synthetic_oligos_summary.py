#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
import glob
import re
import os
import time


def main():
    parser = argparse.ArgumentParser(description="Derive a summary of the synthetic oligos count")
    parser.add_argument("--sampleinfopath", type=str)
    parser.add_argument("--samplesheetpath", type=str)
    args = parser.parse_args()
    sampleinfo = args.sampleinfopath

    timestr = time.strftime("%Y%m%d-%H%M%S")
    
    synthetic_df = pd.DataFrame(columns=['Sample', 'Synthetic oligos', 'Read count', 'Dedup read count', 'FPKM', 'Dup %'])
    
    for synthetic_oligo_count in glob.glob("*synthetic_oligos_stats.txt"):
        individual_df = pd.read_csv(synthetic_oligo_count, header=0,index_col=None,sep="\t")
        synthetic_df = synthetic_df.append(individual_df, ignore_index=True)

    synthetic_flag(synthetic_df,'0.01')
    print(synthetic_df)

    if sampleinfo is not None:
        sampleinfo_data = pd.read_csv(sampleinfo, header=0, sep="\t",index_col=None)
        synthetic_df = pd.merge(sampleinfo_data, synthetic_df, on="Sample", how='outer').fillna('NA')
    
    synthetic_df.to_csv("synthetic_oligo_summary_" + timestr + ".txt", index = None, sep="\t")

def synthetic_flag(df, threshold):
    df["FPKM"] = df["FPKM"].astype(float)
    max_fpkm = df["FPKM"].max().astype(float)
    print(max_fpkm)
    df['1pcflag'] = np.where((df['FPKM'] <= float(threshold) * max_fpkm), "FLAG", "")

if __name__ == '__main__':
    main()