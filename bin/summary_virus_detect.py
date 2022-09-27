#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import time


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load VSD pipeline results")
    # All the required arguments #
    parser.add_argument("--read_size", type=str)

    args = parser.parse_args()
    readsize = args.read_size

    timestr = time.strftime("%Y%m%d-%H%M%S")

    run_data = pd.DataFrame()
    for fl in glob.glob("*nt.blastn.summary.spp.txt"):
        sample_data = pd.read_csv(fl, header=0, sep="\t", index_col=None, delim_whitespace=False)
        run_data = run_data.append(sample_data)
    run_data = run_data[["Sample","Reference","Length","%Coverage","#contig","Depth","Depth_Norm","%Identity","%Identity_max", "%Identity_min","Genus", "Description","Species"]]
    run_data = run_data.astype({'Sample': 'str', 'Reference': 'str','Length': 'int', '%Coverage': 'str' ,'#contig': 'int', 'Depth': 'float', 'Depth_Norm': 'float', '%Identity': 'float', '%Identity_max': 'float', '%Identity_min': 'float', 'Genus': 'str', 'Description': 'str', 'Species': 'str'})
    run_data = run_data.sort_values(["Sample", "Reference"], ascending = (True, True))
    run_data.to_csv("run_summary_top_scoring_targets_virusdetect_"  + readsize + '_' + timestr + ".txt", index=None, sep="\t",float_format="%.2f")
    
    run_data_filtered = pd.DataFrame()
    for flf in glob.glob("*nt.blastn.summary.filtered.txt"):
        sample_data_filtered = pd.read_csv(flf, header=0, sep="\t",index_col=None)
        run_data_filtered = run_data.append(sample_data_filtered)
    print (run_data_filtered)
    run_data_filtered = run_data_filtered[["Sample","Reference","Length","%Coverage","#contig","Depth","Depth_Norm","%Identity","%Identity_max", "%Identity_min","Genus", "Description","Species"]]
    run_data_filtered = run_data_filtered.astype({'Sample': 'str', 'Reference': 'str','Length': 'int', '%Coverage': 'str' ,'#contig': 'int', 'Depth': 'float', 'Depth_Norm': 'float', '%Identity': 'float', '%Identity_max': 'float', '%Identity_min': 'float', 'Genus': 'str', 'Description': 'str', 'Species': 'str'})
    run_data_filtered = run_data_filtered[~run_data_filtered["Species"].str.contains("pararetrovirus")]
    run_data_filtered.drop_duplicates(inplace=True)
    idx = run_data_filtered.groupby(["Species"])["Length"].transform(max) == run_data_filtered["Length"]
    run_data_filtered = run_data_filtered[idx]
    run_data_filtered = run_data_filtered.sort_values(["Sample", "Reference"], ascending = (True, True))
    
    run_data_filtered.to_csv("run_summary_top_scoring_targets_virusdetect_filtered_"  + readsize + "_" + timestr + ".txt", index=None, sep="\t",float_format="%.2f")

if __name__ == "__main__":
    main()