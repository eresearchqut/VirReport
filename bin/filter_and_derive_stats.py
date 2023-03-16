#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from functools import reduce
from glob import glob
from subprocess import run, PIPE

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--results", type=str)
    parser.add_argument("--rawfastq", type=str)
    parser.add_argument("--fastqfiltbysize", type=str)
    parser.add_argument("--sample", type=str)
    parser.add_argument("--read_size", type=str)
    parser.add_argument("--taxonomy", type=str)
    parser.add_argument("--blastdbpath", type=str)
    parser.add_argument("--dedup", type=str)
    parser.add_argument("--cpu", type=str)
    parser.add_argument("--mode", type=str)
    args = parser.parse_args()
    
    results_path = args.results
    sample = args.sample
    rawfastq = args.rawfastq
    fastqfiltbysize = args.fastqfiltbysize
    read_size = args.read_size
    taxonomy = args.taxonomy
    blastdbpath = args.blastdbpath
    dedup = args.dedup
    cpus = args.cpu
    mode = args.mode

    
    if mode == "ncbi":
        raw_data = pd.read_csv(results_path, header=0, sep="\t",index_col=None)
        if len(raw_data) == 0:
            print("DataFrame is empty!")
            csv_file1 = open(sample + "_" + read_size + "_all_targets_with_scores.txt", "w")
            csv_file1.write("sacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\nSpecies\tRNA_type\tSpecies_updated\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score")
            csv_file1.close()
            csv_file2 = open(sample + "_" + read_size + "_top_scoring_targets.txt", "w")
            csv_file2.write("sacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tSpecies\tRNA_type\tSpecies_updated\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score") 
            csv_file2.close()
            csv_file3 = open(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats.txt", "w")
            if dedup == "true": 
                csv_file3.write("Sample\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tSpecies\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score\tmean_read_depth\tread_count\tdedup_read_count\tduplication_rate\tRPM\tFPKM\tPCT_1X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            else:
                csv_file3.write("Sample\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tSpecies\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score\tmean_read_depth\tread_count\tRPM\tFPKM\tPCT_1X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            csv_file3.close()
            exit ()

        #load list of target viruses and viroids and matching official ICTV name

        if os.stat(taxonomy).st_size == 0:
            print('Taxonomy description file is empty!')
            exit ()
        else:
            taxonomy_df = pd.read_csv(taxonomy, header=None, sep="\t")
            taxonomy_df.columns =["sacc", "Species"]
            taxonomy_df["Species"] = taxonomy_df["Species"].str.replace("Hop_stunt_viroid_-_citrus","Hop_stunt_viroid")
            taxonomy_df["Species"] = taxonomy_df["Species"].str.replace("Hop_stunt_viroid;Hop_stunt_viroid","Hop_stunt_viroid")
        #print(taxonomy_df)

        #print(raw_data).head(10)
        print("Cleaning up the data")
        print("Remove double spacing")
        raw_data = raw_data.replace("\s+", " ", regex=True)

        print("Remove hyphens")
        raw_data["stitle"] = raw_data["stitle"].str.replace("-", " ")

        print("Remove underscores")
        raw_data["stitle"] = raw_data["stitle"].str.replace(",_"," ")
        raw_data["stitle"] = raw_data["stitle"].str.replace("_"," ")

        print("Remove commas")
        raw_data["stitle"] = raw_data["stitle"].str.replace(","," ")

        print("Remove problematic text in virus description names")
        #this fixes virus names like "Prunus_necrotic_ringspot_virus_Acot_genomic_RNA,_segment_RNA2,_complete_sequence"
        raw_data["stitle"] = raw_data["stitle"].str.replace(" genomic RNA segment","segment")
        
        raw_data["stitle"] = raw_data["stitle"].str.replace("{complete viroid sequence}","", regex=False)
        #This is a complete genome
        raw_data["stitle"] =  raw_data["stitle"].str.replace("Rubus yellow net virus isolate Canadian 2 hypothetical protein genes  partial cds; hypothetical proteins  polyprotein  ORF 6  and hypothetical protein genes  complete cds; and hypothetical protein genes  partial cds","Rubus yellow net virus isolate Canadian 2 complete cds")
        #remove resistance genes from list of results
        raw_data = raw_data[~raw_data["stitle"].str.contains("resistance gene")]
        raw_data = raw_data[~raw_data["stitle"].str.contains("resistance protein")]
        raw_data = raw_data[~raw_data["stitle"].str.contains("pararetrovirus")]
        raw_data = raw_data[~raw_data["stitle"].str.contains("transposon")]
        raw_data = raw_data[~raw_data["stitle"].str.contains("Petunia vein clearing virus like nonautonomous isolate")]
        
        raw_data = pd.merge(raw_data, taxonomy_df, on=["sacc"])
        raw_data["Species"] = raw_data["Species"].str.replace("_", " ")

        raw_data = raw_data.sort_values("stitle")

        print("If present in original nomenclature, add RNA type information to virus standardised species name")

        raw_data["RNA_type"] = np.where(raw_data.stitle.str.contains("RNA1|RNA 1|segment 1|polyprotein P1"), "RNA1",
                            np.where(raw_data.stitle.str.contains("RNA2|RNA 2|segment 2|polyprotein P2"), "RNA2",
                            np.where(raw_data.stitle.str.contains("RNA3|RNA 3|segment 3|polyprotein P3"), "RNA3", "NaN")))
        
        raw_data["Species_updated"] = raw_data[["Species", "RNA_type"]].agg(" ".join, axis=1)
        #final_data = final_data[~((final_data["Species"].duplicated(keep=False))&(final_data["RNA_type"].str.contains("NaN")))]

        raw_data["Species_updated"] = raw_data["Species_updated"].astype(str).str.replace("NaN", "")
        
        #there are some instances where the species generic name incorporates an RNA type, this fix will catch those cases
        raw_data["Species_updated"] = raw_data["Species_updated"].astype(str).str.replace("RNA1 RNA1", "RNA1")
        raw_data["Species_updated"] = raw_data["Species_updated"].astype(str).str.replace("RNA2 RNA2", "RNA2")
        raw_data["Species_updated"] = raw_data["Species_updated"].astype(str).str.replace("RNA3 RNA3", "RNA3")
        
        raw_data["Species_updated"] = raw_data["Species_updated"].astype(str).str.rstrip( )
        raw_data = raw_data.reset_index(drop=True)
        print (len(raw_data.Species.value_counts()))

        if len(raw_data.Species.value_counts()) == 0:
            print ("Dataframe has no targetted viruses or viroids")
            csv_file1 = open(sample + "_" + read_size + "_all_targets_with_scores.txt", "w")
            csv_file1.write("sacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\nSpecies\tRNA_type\tSpecies_updated\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score")
            csv_file1.close()
            csv_file2 = open(sample + "_" + read_size + "_top_scoring_targets.txt", "w")
            csv_file2.write("sacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tSpecies\tRNA_type\tSpecies_updated\tnaccs_score\tlength_score\tavpid_score\tcov_score\tcompleteness_score\ttotal_score")
            csv_file2.close()
            csv_file3 = open(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats.txt", "w")
            if dedup == "true":
                csv_file3.write("Sample\tSpecies\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tmean_read_depth\tread_count\tdedup_read_count\tduplication_rate\tRPM\tFPKM\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            else:
                csv_file3.write("Sample\tSpecies\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tmean_read_depth\tread_count\tRPM\tFPKM\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            csv_file3.close()
            exit ()

        print("Applying scoring to blast results to select best hit")
        raw_data["naccs"] = raw_data["naccs"].astype(int)
        raw_data["naccs_score"] = raw_data.groupby("Species_updated").apply(max_naccs)

        raw_data["length"] = raw_data["length"].astype(int)
        raw_data["length_score"] = raw_data.groupby("Species_updated").apply(max_length)

        raw_data["av-pident"] = raw_data["av-pident"].astype(float)
        raw_data["avpid_score"] = raw_data.groupby("Species_updated").apply(max_avpid)
        
        raw_data["cov"] = raw_data["cov"].astype(float)
        raw_data["cov_score"] = raw_data.groupby("Species_updated").apply(max_cov)

        raw_data["completeness_score"] = raw_data["stitle"].apply(completeness_score)
        
        raw_data["total_score"] = raw_data["length_score"] + raw_data["naccs_score"] + raw_data["avpid_score"] + raw_data["cov_score"] + raw_data["completeness_score"].astype(int)
        
        print("Output all hits that match species of interest")
        raw_data.to_csv(sample + "_" + read_size + "_all_targets_with_scores.txt", index=None, sep="\t" )

        print("Remove seconday hits based on contig name")
        unique_contigs = list(set([i.strip() for i in ",".join(raw_data["qseqids"]).split(",")]))

        filtered_data = pd.DataFrame()
        for contig in unique_contigs:
            selected = pd.DataFrame()
            selected = raw_data[raw_data["qseqids"].str.contains(contig)]
            
            if len(selected) == 1:
                filtered_data = filtered_data.append(selected)
            #If contigs hit to multiple viruses and viroids, choose best hit
            elif len(selected)>1:
                print(selected)
                # Extract list of spp for a given contig
                Species_updated_list = selected["Species_updated"].tolist()
                #This should accomodate several RNAs per virus spp.
                if len(Species_updated_list) == 1:
                    filtered_data = filtered_data.append(selected)
                # If there are several contigs, retain the top hit, remove 2ary hits
                elif len(Species_updated_list) > 1:
                    topmatch = selected["naccs"].max()
                    selected = selected[selected["naccs"] == topmatch]
                    # Check if there is a tie when selecting by max naccs:
                    Species_updated_list = selected["Species_updated"].tolist()
                    if len(Species_updated_list) == 1:
                        filtered_data = filtered_data.append(selected)
                    # If there is a tie, select next based on av-pidentity
                    else:
                        topmatch = selected["av-pident"].max()
                        selected = selected[selected["av-pident"] == topmatch]
                        filtered_data = filtered_data.append(selected)
                
        filtered_data = filtered_data.drop_duplicates()

        print("Only retain the top hits")
        idx = filtered_data.groupby(["Species_updated"])["total_score"].transform(max) == filtered_data["total_score"]
        filtered_data = filtered_data[idx]
        #print(filtered_data.dtypes)

        #select one random hit if tie for top hits:
        print("If there is a tie, select a random sequence out of the top scoring hit")
        filtered_data = filtered_data.drop_duplicates(subset="Species_updated", keep="first")
    
        #By setting keep on False, all duplicates are True
        #If there are duplicates in species name (ie RNA types present), then it will drop NaN
        filtered_data = filtered_data[~((filtered_data["Species"].duplicated(keep=False))&(filtered_data["RNA_type"].str.contains("NaN")))]
        final_data = filtered_data.drop(["Species"], axis=1)
        final_data = final_data.rename(columns={"Species_updated": "Species"})
        final_data.to_csv(sample + "_" + read_size + "_top_scoring_targets.txt", index=None, sep="\t")

        target_dict = {}
        target_dict = pd.Series(filtered_data.Species_updated.values,index=filtered_data.sacc).to_dict()
        print(filtered_data)
        print (target_dict)
        filtered_data = filtered_data[["sacc","Species","Species_updated","naccs","length","slen","cov","av-pident","stitle","qseqids","contig_ind_lengths","cumulative_contig_len","contig_lenth_min","contig_lenth_max","longest_contig_fasta","total_score"]]
        print(filtered_data)
        #cov_stats (blastdbpath, cpus, dedup, fastqfiltbysize, filtered_data, rawfastq, read_size, sample, target_dict, mode, diagno)
        cov_stats (blastdbpath, cpus, dedup, fastqfiltbysize, filtered_data, rawfastq, read_size, sample, target_dict, mode)

    elif mode == "viral_db":
        final_data = pd.read_csv(results_path, header=0, sep="\t",index_col=None)
        final_data = final_data.rename(columns={"Species": "Species_updated"})
        if len(final_data) == 0:
            print("DataFrame is empty!")
            #if diagno == "true":
                #extension = ("_top_scoring_targets_with_cov_stats_viral_db.txt", "_top_scoring_targets_with_cov_stats_viral_db_regulated.txt", "_top_scoring_targets_with_cov_stats_viral_db_endemic.txt")
                #for ext in extension:
            outfile = open(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats_viral_db.txt", 'w')
            if dedup == "true":
                outfile.write("Sample\tSpecies\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tICTV_information\tmean_read_depth\tread_count\tdedup_read_count\tduplication_rate\tRPM\tFPKM\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            else:
                outfile.write("Sample\tSpecies\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tlongest_contig_fasta\tICTV_information\tmean_read_depth\tread_count\tRPM\tFPKM\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_20X\tconsensus_fasta")
            outfile.close()

            exit ()
        final_data["stitle"] = final_data["stitle"].str.replace("\._", "_")
        target_dict = {}
        target_dict = pd.Series(final_data.Species_updated.values,index=final_data.sacc).to_dict()
        print (target_dict)

        cov_stats (blastdbpath, cpus, dedup, fastqfiltbysize, final_data, rawfastq, read_size, sample, target_dict, mode)

def cov_stats(blastdbpath, cpus, dedup, fastqfiltbysize, final_data, rawfastq, read_size, sample, target_dict, mode):
    print("Align reads and derive coverage and depth for best hit")
    rawfastq_read_counts = (len(open(rawfastq).readlines(  ))/4)


    cov_dict = {}
    dedup_read_counts_dict = {}
    dup_pc_dict = {}
    fpkm_dict = {}
    PCT_1X_dict = {}
    PCT_5X_dict = {}
    PCT_10X_dict = {}
    PCT_20X_dict = {}
    read_counts_dict = {}
    rpm_dict = {}
    consensus_dict = {}

    read_counts_dedup_df = pd.DataFrame()
    dup_pc_df = pd.DataFrame()
    cov_df = pd.DataFrame()
    read_counts_df = pd.DataFrame()
    rpm_df = pd.DataFrame()
    fpkm_df = pd.DataFrame()
    PCT_1X_df = pd.DataFrame()
    PCT_5X_df = pd.DataFrame()
    PCT_10X_df = pd.DataFrame()
    PCT_20X_df = pd.DataFrame()
    consensus_df = pd.DataFrame()

    
    for refid, refspname in target_dict.items():
        try:
            print (refid)
            print (refspname)
            combinedid = str(refid + " " + refspname).replace("sp.","sp").replace(" ","_")

            print("Extract sequence from blast database")
            fastafile = (sample + "_" + read_size + "_" + combinedid + ".fa").replace(" ","_")
            
            if mode == "ncbi":
                single_fasta_entry = open(fastafile, "w")
                command_line = ["blastdbcmd","-db", blastdbpath, "-entry", refid, \
                                "-outfmt","'%f'"]
                subprocess.call(command_line, stdout=single_fasta_entry)

                filesize = os.path.getsize(fastafile)

                if filesize == 0:
                    print("Retrieval from blast db failed")
                    single_fasta_entry = open(fastafile, "w")
                    p1 = subprocess.Popen(["esearch", "-db", "nucleotide", "-query", refid], stdout=subprocess.PIPE)
                    p2 = subprocess.run(["efetch", "-format", "fasta"], stdin=p1.stdout, stdout=single_fasta_entry)
                    single_fasta_entry.close()

            elif mode == "viral_db":
                #p1 = subprocess.Popen(["esearch", "-db", "nucleotide", "-query", refid], stdout=subprocess.PIPE)
                #p2 = subprocess.run(["efetch", "-format", "fasta"], stdin=p1.stdout, stdout=single_fasta_entry)
                #p1 = subprocess.Popen(["grep", "-A1", refid, blastdbpath], stdout=single_fasta_entry)
                #single_fasta_entry.close()
                
                bowtie_index = ["grep", "-A1", refid, blastdbpath]
                subprocess.call(bowtie_index, stdout=open(fastafile,"w"))

            print("Building a bowtie index")
            index=(sample + "_" + read_size + "_" + combinedid).replace(" ","_")
            buildindex = ["bowtie-build","-f", fastafile, index]
            subprocess.call(buildindex)

            print("Aligning original reads")
            samoutput = str(index + ".sam")
            bowtie_output = str(index + "_bowtie_log.txt")
            aligning = ["bowtie", "-q", "-v", "2", "-k", "1", "-p", cpus , "-x", index, fastqfiltbysize, "-S", samoutput]
            subprocess.call(aligning, stderr=open(bowtie_output,"w"))

            print("Derive a bam file")
            bamoutput = str(index + ".bam")
            derivebam = ["samtools", "view", "-@", cpus, "-bS", samoutput]
            subprocess.call(derivebam, stdout=open(bamoutput,"w"))

            subprocess.call(["rm","-r", samoutput])

            print("Sorting bam file")
            sortedbamoutput = str(index + ".sorted.bam")
            sorting = ["samtools", "sort", "-@", cpus, bamoutput, "-o", sortedbamoutput]
            subprocess.call(sorting)

            print("Indexing bam file")
            bamindex = str(index + ".sorted.bam.bai")
            indexing = ["samtools", "index", sortedbamoutput]
            subprocess.call(indexing, stdout=open(bamindex,"w"))

            read_counts = ()
            with open(bowtie_output) as bo:
                a = " "
                while(a):
                    a = bo.readline()
                    l = a.find("# reads with at least one alignment:") #Gives a non-negative value when there is a match
                    if ( l >= 0 ):
                        print(a)
                        read_counts = a.split(" ")[7]
            read_counts_dict[refspname] = read_counts

            subprocess.call(["rm","-r", bamoutput])

            #If data needs to be deduplicated
            dedupbamoutput = str(index + ".dedup.bam")
            umi_dedup_log = str(index + "_umi_tools.log")
            dedupbamindex = str(index + ".dedup.bam.bai")
            dedup_read_counts = ()
            final_read_counts = ()
            dup_pc = ()
            finalbamoutput = ()
            finalbamindex = ()
            rpm = ()
            fpkm = ()

            if dedup == "true":
                print("Deduping bam file")
                umitools_dedup = ["umi_tools", "dedup", "-I", sortedbamoutput, "--method", "unique",  "-L", umi_dedup_log]
                subprocess.call(umitools_dedup, stdout=open(dedupbamoutput,"w"))
                
                print("Indexing dedup bam file")
                dedup_indexing = ["samtools", "index", dedupbamoutput]
                subprocess.call(dedup_indexing, stdout=open(dedupbamindex,"w"))
                p = run(["samtools", "view", "-c", "-F", "260", dedupbamoutput], stdout=PIPE, encoding='ascii')
                
                dedup_read_counts = int(p.stdout.replace("\n",""))
                dedup_read_counts_dict[refspname] = dedup_read_counts
                print(dedup_read_counts_dict)
                
                dup_pc = round(100-(int(dedup_read_counts)*100/int(read_counts)))
                dup_pc_dict[refspname] = dup_pc
                
                read_counts_dedup_df = pd.DataFrame(dedup_read_counts_dict.items(),columns=["Species_updated", "dedup_read_count"]) 
                dup_pc_df = pd.DataFrame(dup_pc_dict.items(),columns=["Species_updated", "duplication_rate"])
                finalbamoutput = dedupbamoutput
                finalbamindex = dedupbamindex
                final_read_counts = dedup_read_counts            
                
                subprocess.call(["rm","-r", sortedbamoutput])
                subprocess.call(["rm","-r", bamindex])

            if dedup == "false":
                finalbamoutput = sortedbamoutput
                finalbamindex = bamindex
                final_read_counts = read_counts

            pileup = str(index + ".pileup")
            derivepileup = ["samtools", "mpileup", "-uf", fastafile, finalbamoutput, "-o", pileup]
            #subprocess.call(derivepileup, stdout=pileup)
            subprocess.call(derivepileup)

            #variant calling
            vcfout = str(index + ".vcf.gz")
            vcfcall = ["bcftools", "call", "-c", pileup, "-Oz", "-o", vcfout]
            subprocess.call(vcfcall)
            
            vcfindex = ["bcftools", "index", vcfout]
            subprocess.call(vcfindex)

            # Normalise indels:
            bcfnormout = str(index + "_norm.bcf")
            bcfnorm = ["bcftools", "norm", "-f", fastafile, vcfout, "-Ob", "-o", bcfnormout]
            subprocess.call(bcfnorm)
            bcfnormoutindex = ["bcftools", "index", bcfnormout]
            subprocess.call(bcfnormoutindex)

            # Filter adjacent indels within 5bp
            bcfnormoutfiltout = str(index + "_norm_flt_indels.bcf")
            bcfnormoutfilt = ["bcftools", "filter", "--IndelGap", "5", bcfnormout, "-Ob", "-o", bcfnormoutfiltout]
            subprocess.call(bcfnormoutfilt)
            bcfnormoutfiltindex = ["bcftools", "index", bcfnormoutfiltout]
            subprocess.call(bcfnormoutfiltindex)

            # Convert bcf to vcf
            vcfnormoutfiltout = str(index + "_sequence_variants.vcf.gz")
            vcfnormoutfilt = ["bcftools", "view", "-Oz", "-o", vcfnormoutfiltout, bcfnormoutfiltout]
            subprocess.call(vcfnormoutfilt)
            vcfnormoutfiltindex = ["bcftools", "index", vcfnormoutfiltout]
            subprocess.call(vcfnormoutfiltindex)

            # Get consensus fasta file
            genomecovbed = str(index + "_genome_cov.bed")
            gencovcall = ["bedtools", "genomecov", "-ibam", finalbamoutput, "-bga"]
            subprocess.call(gencovcall, stdout=open(genomecovbed,"w"))

            # Assign N to nucleotide positions that have zero coverage
            zerocovbed = str(index + "_zero_cov.bed")
            zerocovcall = ["awk", "$4==0 {print}", genomecovbed]
            subprocess.call(zerocovcall, stdout=open(zerocovbed,"w"))

            maskedfasta = (sample + "_" + read_size + "_" + combinedid + "_masked.fa").replace(" ","_")
            maskedfastaproc = ["bedtools", "maskfasta", "-fi",  fastafile, "-bed", zerocovbed, "-fo", maskedfasta]
            subprocess.call(maskedfastaproc)

            # Derive a consensus fasta file
            consensus = str(index + ".consensus.fasta")
            consensuscall = ["bcftools", "consensus", "-f",  maskedfasta, vcfout, "-o", consensus]
            subprocess.call(consensuscall)

            consensus_seq = ""
            with open(consensus, 'r') as f:
                for line in f:
                    if line[0] == ">":
                        consensus_seq += line.strip()
                        consensus_seq += ' '
                    else:
                        consensus_seq += line.strip()

            consensus_seq = consensus_seq.replace('"', '')
            consensus_dict[refspname] = consensus_seq
            print(consensus_seq)

            subprocess.call(["rm","-r", pileup])
            subprocess.call(["rm","-r", vcfout])
            subprocess.call(["rm","-r", genomecovbed])
            subprocess.call(["rm","-r", zerocovbed])
            subprocess.call(["rm","-r", maskedfasta])
            subprocess.call(["rm","-r", bcfnormout])
            subprocess.call(["rm","-r", bcfnormout + '.csi'])
            subprocess.call(["rm","-r", bcfnormoutfiltout])
            subprocess.call(["rm","-r", bcfnormoutfiltout + '.csi'])
            

            # Derive Picard statistics 
            print("Running picard")
            picard_output = (index + "_picard_metrics.txt")
            picard = ["picard", "CollectWgsMetrics", "-I", str(finalbamoutput), "-O", str(picard_output), "-R", str(fastafile), "-READ_LENGTH","22", "-COUNT_UNPAIRED", "true"]
            subprocess.call(picard)

            reflen = ()
            cov = ()
            PCT_1X = ()
            PCT_5X = ()
            PCT_10X = ()
            PCT_20X = ()
            
            with open(picard_output) as f:
                a = " "
                while(a):
                    a = f.readline()
                    l = a.find("MEAN_COVERAGE") #Gives a non-negative value when there is a match
                    if ( l >= 0 ):
                        line = f.readline()
                        elements = line.split("\t")
                        reflen, cov, PCT_1X, PCT_5X, PCT_10X, PCT_20X = elements[0], elements[1], elements[13], elements[14], elements[15],elements[17]
            f.close()
            cov_dict[refspname] = cov
            PCT_1X_dict[refspname] = PCT_1X
            PCT_5X_dict[refspname] = PCT_5X
            PCT_10X_dict[refspname] = PCT_10X
            PCT_20X_dict[refspname] = PCT_20X

            fpkm = round(int(final_read_counts)/(int(reflen)/1000*int(rawfastq_read_counts)/1000000))
            rpm = round(int(final_read_counts)*1000000/int(rawfastq_read_counts))

            rpm_dict[refspname] = rpm
            fpkm_dict[refspname] = fpkm

            cov_df = pd.DataFrame(cov_dict.items(),columns=["Species_updated", "mean_read_depth"])
            read_counts_df = pd.DataFrame(read_counts_dict.items(),columns=["Species_updated", "read_count"])
            rpm_df = pd.DataFrame(rpm_dict.items(),columns=["Species_updated", "RPM"])
            fpkm_df = pd.DataFrame(fpkm_dict.items(),columns=["Species_updated", "FPKM"])
            PCT_1X_df = pd.DataFrame(PCT_1X_dict.items(),columns=["Species_updated", "PCT_1X"])
            PCT_5X_df = pd.DataFrame(PCT_5X_dict.items(),columns=["Species_updated", "PCT_5X"])
            PCT_10X_df = pd.DataFrame(PCT_10X_dict.items(),columns=["Species_updated", "PCT_10X"])
            PCT_20X_df = pd.DataFrame(PCT_20X_dict.items(),columns=["Species_updated", "PCT_20X"])
            consensus_df = pd.DataFrame(consensus_dict.items(),columns=["Species_updated", "consensus_fasta"])

            project_files = glob(index + "*ebwt") + glob(index + ".vcf.gz*")
            for fl in project_files:
                subprocess.call(["rm","-r", fl])
        except OSError as err:
            print("OS error: {0}".format(err))

    print("Deriving summary table with coverage statistics")

    if read_counts_dedup_df.empty:
        dfs = [final_data, cov_df, read_counts_df, rpm_df, fpkm_df, PCT_1X_df, PCT_5X_df, PCT_10X_df, PCT_20X_df, consensus_df]
    else:
        dfs = [final_data, cov_df, read_counts_df, read_counts_dedup_df, dup_pc_df, rpm_df, fpkm_df, PCT_1X_df, PCT_5X_df, PCT_10X_df, PCT_20X_df, consensus_df]

    full_table = reduce(lambda left,right: pd.merge(left,right,on=["Species_updated"],how='outer'), dfs)

    full_table["mean_read_depth"] = full_table["mean_read_depth"].astype(float)
    full_table["PCT_1X"] = full_table["PCT_1X"].astype(float)
    full_table["PCT_5X"] = full_table["PCT_5X"].astype(float)
    full_table["PCT_10X"] = full_table["PCT_10X"].astype(float)
    full_table["PCT_20X"] = full_table["PCT_20X"].astype(float)
    if "duplication_rate" in full_table.columns:
        full_table["duplication_rate"] = full_table["duplication_rate"].astype(float)
    full_table.insert(0, "Sample", sample)

    if mode == 'ncbi':
        full_table = full_table.drop(["Species"], axis=1)
        full_table = full_table.rename(columns={"Species_updated": "Species"})
        full_table.to_csv(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats.txt", index=None, sep="\t",float_format="%.2f")
       
    elif mode == 'viral_db':
        full_table = full_table.rename(columns={"Species_updated": "Species"})
        full_table.to_csv(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats_viral_db.txt", index=None, sep="\t",float_format="%.2f")
    

def max_avpid(df):
    max_row = df["av-pident"].max()
    labels = np.where((df["av-pident"] == max_row),
                    "1",
                    "0")
    return pd.DataFrame(labels, index=df.index).astype(int)

def max_length(df):
    max_row = df["length"].max()
    labels = np.where((df["length"] == max_row),
                        "2",
                        "0")
    return pd.DataFrame(labels, index=df.index).astype(int)

def max_naccs(df):
    max_row = df["naccs"].max()
    labels = np.where((df["naccs"] == max_row),
                        "1",
                        "0")
    return pd.DataFrame(labels, index=df.index).astype(int)

def max_cov(df):
    max_row = df["cov"].max()
    labels = np.where((df["cov"] == max_row),
                        "1",
                        "0")
    return pd.DataFrame(labels, index=df.index).astype(int)

def completeness_score(x):
    if "complete sequence" in str(x):
        return 3
    elif "complete genome" in str(x):
        return 3
    elif "polyprotein gene complete cds" in str(x):
        return 3
    elif "polyprotein 1 gene complete cds" in str(x):
        return 3
    elif "polyprotein 2 gene complete cds" in str(x):
        return 3
    elif "nearly complete sequence" in str(x):
        return -3
    if "partial" in str(x):
        return -3
    elif "polymerase protein" in str(x):
        return -3
    elif "RNA-dependent RNA polymerase" in str(x):
        return -3
    else:   
        return 0

if __name__ == "__main__":
    main()
