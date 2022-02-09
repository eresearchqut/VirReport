# VirReport workflow
Roberto Barrero, 14/03/2019  
Desmond Schmidt, 2/7/2019  
Converted to Nextflow by Craig Windell, 11/2020  
Modified by Maely Gauthier 12/2021  

## About Pipeline
VirReport pipeline based on the scientific workflow manager Nextflow.
It is designed to help phytosanitary diagnostics of viruses and viroid pathogens in quarantine facilities. It takes small RNA-Seq samples as input.

## Run the test

## Run your own analysis

Run the command:
```nextflow run eresearchqut/VirReport -profile {docker or singularity or conda} --indexfile $PBS_O_WORKDIR/index_example.csv```

Set the profile parameter to one of
```
docker
singularity
conda
```
To suite your environment.

The VSD workflow will perform the following steps by default:
- Retain reads of a given length (e.g. 21-22 or 24 nt long) from fastq file(s) provided in index.csv file (readprocessing)  
- De novo assembly using kmer 15 and coverage 3 (velvet) - 
- Collapse contigs into scaffolds (min length 20) (cap3)
- Run megablast homology search against NCBI NT database (blastn_nt_velvet)
- Summarise megablast results and restrict to virus and viroid matches (BlastTools_blastn_velvet)
- Derive coverage statistics, consensus sequence and VCF matching to top blast hits (filter_n_cov)

A number of additional optional steps can be run:
```
     --blastp: Predict ORF from de novo assembly (derived with Velvet) and run blastP againts NCBI NR (getorf, blastp, blastpdbcmd, BlastToolsp) --blastp

     --contamination_detection: Run cross-sample contamination prediction (contamination_detection) 

     --blastlocaldb: Run blastn and megablast homology search on de novo assembly (derived with Velvet) against local  virus and viroid database (blast_nt_localdb_velvet, filter_blast_nt_localdb_velvet)

     --blastn_method: The blastn homology search can be specified as blastn instead of megablast (--blastn_method blastn)

     --spades: Run SPAdes 3.14 de novo assembler and perform blastn homology analysis on the derived de novo contigs (spades, cap3_spades, megablast_nt_spades , BlastToolsn_megablast_spades)
```
A number of additional options are included:
```
    --targets: A text file with the taxonomy of the viruses/virioids of interest can be provided and only these will be retained in the megablast summary results derived at the filter_n_cov step.

    --spadeskmer specifies the range of kmers to use when running spades

    --cap3_len specifies the minimal length of contigs to retain after CAP3 scaffolding

    --blastn_evalue and --blastp_evalue specifies the evalue parameter to use during blast analyses

    --orf_minsize correspond to the minimal open reading frame getorf retains
```
To enable these options, they can either be included in the nextflow run command provided in the PBS script: 
```
nextflow run eresearchqut/VirReport -profile {docker or singularity or conda} --indexfile index_example.csv --blastlocaldb --spades --contamination_detection
```
or update parameter to true in the nextflow.config file. For instance:
```
params {
blastlocaldb = true
spades = true
contamination_detection = true
}
```
## Preparing the data
Preparing an index.csv file

You need to create a TAB delimited text file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the --indexfile [filename] in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,samplepath,minlen,maxlen:```

- sampleid will be the sample name that will be given to the files created by the pipeline
- samplepath is the full path to the quality filtered fastq files that the pipeline requires as starting input
- minlen and maxlen correspond to the read size that will be retained for downstream analyses. 

An index_example.csv is included in the base directory:
```
sampleid,samplepath,minlen,maxlen
MT212,/work/hia_mt18005/diagnostics/2021/14_RAMACIOTTI_LEL9742-LEL9751/results/06_usable_reads/MT212_21-22bp.fastq,21,22
MT213,/work/hia_mt18005/diagnostics/2021/14_RAMACIOTTI_LEL9742-LEL9751/results/06_usable_reads/MT213_21-22bp.fastq,21,22
```
You also need to provide the path of your NCBI blast directory in the nextflow.config file. For instance:
```
params {
  blast_db = '/lustre/work-lustre/hia_mt18005/blastDB/30112021'
  }
```
If you are interested to run a blast analysis against a local database, you also need to specify its path in the nextflow.config file. For example:
blastn_local_db = '/work/hia_mt18005/databases/PVirDB/PVriDB_ver2021_11_09/PVirDB_ver20211109.fasta'
Running the pipeline


## Outputs
The folders are structures as follows (examples of outputs are provided in italics):
- 01_read_size_selection (cutadapt log file and fastq file including reads only matching the size specified in the index.csv file) MT020_21-22nt_cutadapt.log & MT020_21-22nt.fastq
- 02_velvet (velvet results and the fasta file which includes the velvet assembled contigs MT020_velvet_assembly_21-22nt.fasta
- 02a_spades (if spades is additionally run)
- 03_cap3 (fasta file of the scaffolds produced by CAP3 as well as the singletons) MT020_velvet_cap3_21-22nt_rename.fasta
- 04_blastn (all blastn results, filtered results limited to only viruses and viroid top 5 hit matches and their taxonomy) MT020_velvet_21-22nt_megablast_vs_NT.bls, MT020_velvet_21-22nt_megablast_vs_NT_top5Hits.txt, MT020_velvet_21-22nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt MT020_velvet_21-22nt_megablast_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
- 05_blastoutputs (BlastTools.jar summary output which clusters all the contigs matching to a specific hit. summary_MT029_velvet_21-22nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt
- 06_blastp (blastp outputs) MT020_velvet_21-22nt_getorf.min50aa.fasta, MT020_velvet_21-22nt_getorf.min50aa_blastp_vs_NR_out_virus_viroid.txt
- 07_filternstats (filtered blast summary with various coverage statistics for each virus and viroid hit, and associated consensus fasta file and vcf file) MT020_21-22nt_top_scoring_targets_with_cov_stats.txt, MT020_21-22nt_MK929590_Peach_latent_mosaic_viroid.consensus.fasta, MT020_21-22nt_MK929590_Peach_latent_mosaic_viroid_sequence_variants.vcf.gz
- 08_summary (summary of results for all samples included in the index.csv file. This includes a cross-contamination prediction) run_top_scoring_targets_with_cov_stats_with_cont_flag_21-22nt_0.01.txt

Future potential additional features:
- Include a deduplication step for fastq files that have UMIs incorporated
- Incorporate the fastq file initial filtering steps from sRNAqc as option
