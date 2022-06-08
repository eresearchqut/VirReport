# VirReport workflow
## Authors
Roberto Barrero, 14/03/2019  
Desmond Schmidt, 2/7/2019  
Converted to Nextflow by Craig Windell, 11/2020  
Modified by Maely Gauthier 12/2021  

## About Pipeline
VirReport pipeline is based on the scientific workflow manager Nextflow.
It is designed to help phytosanitary diagnostics of viruses and viroid pathogens in quarantine facilities. It takes small RNA-Seq samples as input.

# Run the Pipeline

## Test run
Download the pipeline and test it on a minimal dataset with a single command:
```bash
nextflow -c conf/test.config run main.nf -profile test,conda --dedup -resume --contamination_detection
```

This dataset can be run locally. It requires 2 cpus and 8 Gb of memory and should take less than 5 mins to complete.

## Run with your own data

Run the command:
```bash
nextflow run eresearchqut/VirReport -profile {docker or singularity or conda} --indexfile $PBS_O_WORKDIR/index_example.csv
```

Set the profile parameter to one of
```
docker
singularity
conda
```
To suit your environment.

The VirReport workflow will perform the following steps by default:
- Retain reads of a given length (e.g. 21-22 or 24 nt long) from fastq file(s) provided in index.csv file (readprocessing)  
- De novo assembly using kmer 15 and coverage 3 (velvet) and kmer range 19-21 (SPAdes)
- Collapse contigs into scaffolds (min length 30) (cap3)
- Run megablast homology search against NCBI NT database (blastn_nt_cap3)
- Summarise megablast results and restrict to virus and viroid matches (BlastTools_blastn_cap3)
- Derive coverage statistics, consensus sequence and VCF matching to top blast hits (filter_n_cov)

A number of additional optional steps can be run:
```
     --blastlocaldb: Run megablast and blastn homology search on de novo assembly against local  virus and viroid database (blast_nt_localdb_cap3, filter_blast_nt_localdb_cap3). An example of local virus database can be downloaded at wget https://data.researchdatafinder.qut.edu.au/dataset/60eed574-a745-4a0f-ab7c-fb8b3c711695/resource/a17dfa13-a093-407a-a047-27f134f92ac9/download/pvirdbv1.fasta.gz
     
     --blastp: Predict ORF from de novo assembly and run blastP againts NCBI NR (getorf, blastp, blastpdbcmd, BlastToolsp) --blastp

     --contamination_detection: Run cross-sample contamination prediction on blastn top hit for NT homology search (contamination_detection) 
     
     --contamination_detection_localdb: Run cross-sample contamination prediction on blastn top hit for homology search against local virus databsase (contamination_detection) 
     
```
A number of additional options are included:
```
    --qualityfilter: performs a quality filtering on raw fastq files (currently specifically written for samples prepared using QIAGEN QIAseq miRNA library kit). The pipeline will also derive a qc report, a run_read_size_distribution.pdf, a read_RNA_source.pdf and a read_origin_pc_summary.txt.

    --targets: A text file with the taxonomy of the viruses/virioids of interest can be provided and only these will be retained in the megablast summary results derived at the filter_n_cov step for the blastn resulsts against NT.

    --cap3_len specifies the minimal length of contigs to retain after CAP3 scaffolding

    --blastn_evalue and --blastp_evalue specifies the evalue parameter to use during blast analyses

    --orf_minsize correspond to the minimal open reading frame getorf retains

    --virusdetect: runs in parallel ViruDetect v1.7
```
To enable these options, they can either be included in the nextflow run command provided in the PBS script: 
```
nextflow run eresearchqut/VirReport -profile {docker or singularity or conda} --indexfile index_example.csv --blastlocaldb --contamination_detection --virusdetect
```
or update parameter to true in the nextflow.config file. For instance:
```
params {
  blastlocaldb = true
  virusdetect = true
  contamination_detection = true
}
```
## Preparing the data
Preparing an index.csv file

You need to create a TAB delimited text file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the --indexfile [filename] in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,samplepath```

- sampleid will be the sample name that will be given to the files created by the pipeline
- samplepath is the full path to the fastq files that the pipeline requires as starting input

An index_example.csv is included in the base directory:
```
sampleid,samplepath
MT212,/work/diagnostics/2021/MT212_21-22bp.fastq
MT213,/work/diagnostics/2021/MT213_21-22bp.fastq
```
By default the pipeline will only retain 21-22 nt sRNA reads for downstream analysis but you can change this range in the nextflow.config file. For instance:
```
params {
  minlen = '18'
  maxlen = '25'
}
```
You also need to provide the path of your NCBI blast nt and/or nr directory/ies in the nextflow.config file. For instance:
```
params {
  blast_db_dir = '/work/hia_mt18005_db/blastDB/20220408'
}
```
If you want to run a blast analysis against a local database, please ensure you use NCBI BLAST+ makeblastdb to create the database. Then specify the full path to the database files including the prefix in the nextflow.config file. For example:
```
params {
  blast_local_db_path = '/work/hia_mt18005/databases/PVirDB/PVirDB_ver2022_02_09/dev/PVirDB_ver20220209.fasta'
}
```
If you want to run the initial qualityfilter step, you will need to specify in the nextflow.config file the directory which holds the required bowtie indices to: 1) derive the origin of the filtered reads obtained and 2) filter non-informative reads (qc_fastq_filtered and derive_usable_reads processes). For example:
params {
  bowtie_db_dir = '/work/hia_mt18005_db/bowtie_idx'
}
```
If you want to run VirusDetect, then specify the path to the viral database directory in nextflow.config file. These can be downloaded at http://bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/. For example:
```
virusdetect_db_path = '/work/hia_mt18005_db/blastDB/VirusDetect/vrl_Plants_239_U97'
```

# Running the pipeline


## Outputs
The folders are structured as follows:
- 00_quality_filtering/Sample_name: this folder will output FASTQC of raw and filtered fastq files, cutadapt and umi_tools log files, a quality_trimmed.fastq.gz file and by default a fastq.gz file including reads only matching the size specified in the nextflow.config file)
- 00_quality_filter/qc_report: this folder contains summaries for all samples included in the index_example.csv file

- 01_VirReport/Sample_name/Assembly: fasta file which includes the assembled contigs before and after CAP3, for example MT020_velvet_assembly_21-22nt.fasta, MT001_cap3_21-22nt.fasta
- 01_VirReport/Sample_name/blastn/NT: this folder contains all blastn results, filtered results limited to only viruses and viroid top 5 hit matches and the final BlastTools.jar summary output. For example: MT020_cap3_21-22nt_megablast_vs_NT.bls, MT020_cap3_21-22nt_megablast_vs_NT_top5Hits.txt, MT020_cap3_21-22nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt, summary_MT029_cap3_21-22nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt
- 01_VirReport/Sample_name/blastn/localdb: analysis using the local db will be saved in this folder, for example:
MT001_cap3_21-22nt_blastn_vs_localdb.bls, summary_MT001_cap3_21-22nt_blastn_vs_localdb.bls_viruses_viroids_ICTV.txt
- 01_VirReport/Sample_name/blastp: this folder contains getorf and blastp outputs. For example: MT020_cap3_21-22nt_getorf.min50aa.fasta, MT020_cap3_21-22nt_getorf.min50aa_blastp_vs_NR_out_virus_viroid.txt
- 01_VirReport/Sample/alignments: this folder ontains filtered blast summary with various coverage statistics for each detected virus and viroid hit, and associated consensus fasta file and vcf file. For example: MT020_21-22nt_top_scoring_targets_with_cov_stats.txt, MT020_21-22nt_MK929590_Peach_latent_mosaic_viroid.consensus.fasta, MT020_21-22nt_MK929590_Peach_latent_mosaic_viroid_sequence_variants.vcf.gz
- 01_VirReport/summary: this folder contains a summary of results for all samples included in the index.csv file. This summay table includes a cross-contamination prediction flag. For example: run_top_scoring_targets_with_cov_stats_with_cont_flag_21-22nt_0.01.txt

- 02_virusdetect/Sample_name: this folder includes a results folder with blastn and blastx summary. For example: MT016_21-22nt.blastn.summary.txt and MT016_21-22nt.blastx.summary.txt
