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
    
    df = pd.DataFrame(columns=['Sample', 'Synthetic oligos', 'Read count', 'Dedup read count', 'FPKM', 'Dup %'])
    
    for synthetic_oligo_count in glob.glob("*synthetic_oligos_stats.txt"):
        individual_df = pd.read_csv(synthetic_oligo_count, header=0,index_col=None,sep="\t")
        df = df.append(individual_df, ignore_index=True)
    
    

    if sampleinfo is not None:
        sampleinfo_data = pd.read_csv(sampleinfo, header=0, sep="\t",index_col=None)
        df = pd.merge(sampleinfo_data, df, on="Sample", how='outer').fillna('NA')
    print(df)
    df.to_csv("synthetic_oligo_summary_" + timestr + ".txt", index = None, sep="\t")
if __name__ == '__main__':
    main()