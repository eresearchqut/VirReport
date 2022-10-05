#!/usr/bin/env python
import argparse
import pandas as pd
import subprocess
from functools import reduce
from subprocess import run, PIPE


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--rawfastq", type=str)
    parser.add_argument("--fastqfiltbysize", type=str)
    parser.add_argument("--sample", type=str)
    parser.add_argument("--read_size", type=str)
    args = parser.parse_args()
    
    sample = args.sample
    rawfastq = args.rawfastq
    fastqfiltbysize = args.fastqfiltbysize
    read_size = args.read_size

    rawfastq_read_counts = (len(open(rawfastq).readlines(  ))/4)
    read_counts_dict = {}
    dedup_read_counts_dict = {}
    dup_pc_dict = {}
    fpkm_dict = {}

    #oligo = str("celmiR39")
    #fastafile = (sample + "_" + read_size + "_" + oligo + ".fa").replace(" ","_")
    #single_fasta_entry = open(fastafile, "w")
    #single_fasta_entry.write(">celmir39\nTCACCGGGTGTAAATCAGCTTG")
    #single_fasta_entry.close()

    fastafile = ("celmiR39.fa")
    single_fasta_entry = open(fastafile, "w")
    single_fasta_entry.write(">celmir39\nTCACCGGGTGTAAATCAGCTTG")
    single_fasta_entry.close()

    fastafile = ("celmiR54.fa")
    single_fasta_entry = open(fastafile, "w")
    single_fasta_entry.write(">celmir39\nTACCCGTAATCTTCATAATCCGAG")
    single_fasta_entry.close()
    
    fastafile = ("celmiR238.fa")
    single_fasta_entry = open(fastafile, "w")
    single_fasta_entry.write(">celmir238\nTTTGTACTCCGATGCCATTCAGA")
    single_fasta_entry.close()

    gene_list = ('celmiR39', 'celmiR54', 'celmiR238')
    for index in gene_list:
        #print("Building a bowtie index")
        #index = (sample + "_" + read_size + "_" + oligo)
        #print(index)
        #buildindex = ["bowtie-build", "-f", fastafile, index]
        buildindex = ["bowtie-build", "-f", index + ".fa", index]
        subprocess.call(buildindex)
        
        print("Aligning original reads")
        samoutput = str(index + ".sam")
        bowtie_output = str(index + "_bowtie_log.txt")
        aligning = ["bowtie", "-q", "-v", "1", "-k", "1", "-p", "4", "-x", index, fastqfiltbysize, "-S", samoutput]
        subprocess.call(aligning, stderr=open(bowtie_output,"w"))

        print("Derive a bam file")
        bamoutput = str(index + ".bam")
        derivebam = ["samtools", "view", "-@", "4", "-bS", samoutput]
        subprocess.call(derivebam, stdout=open(bamoutput,"w"))

        print("Sorting bam file")
        sortedbamoutput = str(index + ".sorted.bam")
        sorting = ["samtools", "sort", "-@", "4", bamoutput, "-o", sortedbamoutput]
        subprocess.call(sorting)

        print("Indexing bam file")
        bamindex = str(index + ".sorted.bam.bai")
        indexing = ["samtools", "index", sortedbamoutput]
        subprocess.call(indexing, stdout=open(bamindex,"w"))

        print("Deduping bam file")
        dedupbamoutput = str(index + ".dedup.bam")
        umi_dedup_log = str(index + "_umi_tools.log")
        dedup = ["umi_tools", "dedup", "-I", sortedbamoutput, "-L", umi_dedup_log]
        subprocess.call(dedup, stdout=open(dedupbamoutput,"w"))

        print("Indexing dedup bam file")
        dedupbamindex = str(index + ".dedup.bam.bai")
        dedup_indexing = ["samtools", "index", dedupbamoutput]
        subprocess.call(dedup_indexing, stdout=open(dedupbamindex,"w"))

        subprocess.call(["rm","-r", samoutput])
        subprocess.call(["rm","-r", bamoutput])
        subprocess.call(["rm","-r", sortedbamoutput])
        subprocess.call(["rm","-r", bamindex])

        reflen = ()
        print("Deriving synthetic oligo length")
        bamheaderout = str(index + "_bam_header.txt")
        header = ["samtools", "view", "-H", dedupbamoutput]
        subprocess.call(header, stdout=open(bamheaderout,"w"))

        with open(bamheaderout, 'r') as f:
            for line in f:
                string_to_search1 = ("@HD")
                if string_to_search1 in line:
                    line = next(f)
                    elements = line.split("LN:")
                    reflen = int(elements[1].strip())
        f.close()

        dup_pc = ()
        read_counts = ()
        dedup_read_counts = ()
        fpkm = ()


        dup_pc_dict[index] = dup_pc
        p = run(["samtools", "view", "-c", "-F", "260", dedupbamoutput], stdout=PIPE, encoding='ascii')
        dedup_read_counts = p.stdout.replace("\n","")
        dedup_read_counts_dict[index] = dedup_read_counts
        print(dedup_read_counts_dict)
        
        with open(bowtie_output) as bo:
            a = " "
            while(a):
                a = bo.readline()
                l = a.find("# reads with at least one alignment:") #Gives a non-negative value when there is a match
                if ( l >= 0 ):
                    print(a)
                    read_counts = a.split(" ")[7]

        read_counts_dict[index] = read_counts
        fpkm = round(int(dedup_read_counts)/(int(reflen)/1000*int(rawfastq_read_counts)/1000000))
        fpkm_dict[index] = fpkm
        dup_pc = round(100-(int(dedup_read_counts)*100/(int(read_counts)+0.1)))
        dup_pc_dict[index] = dup_pc

        read_counts_df = pd.DataFrame(read_counts_dict.items(),columns=["Synthetic oligos", "Read count"])
        read_counts_dedup_df = pd.DataFrame(dedup_read_counts_dict.items(),columns=["Synthetic oligos", "Dedup read count"]) 
        dup_pc_df = pd.DataFrame(dup_pc_dict.items(),columns=["Synthetic oligos", "Dup %"]) 
        fpkm_df = pd.DataFrame(fpkm_dict.items(),columns=["Synthetic oligos", "FPKM"])

        dfs = [read_counts_df, read_counts_dedup_df, fpkm_df, dup_pc_df]
        full_table = reduce(lambda left,right: pd.merge(left,right,on="Synthetic oligos"), dfs)

        full_table["Dup %"] = full_table["Dup %"].astype(float)
        
        full_table.insert(0, "Sample", sample)
        print(full_table)
        full_table.to_csv(sample + "_" + read_size + "_synthetic_oligos_stats.txt", index=None, sep="\t",float_format="%.2f")




if __name__ == "__main__":
    main()
