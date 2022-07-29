# VirReport


## Introduction
eresearchqut/VirReport is a bioinformatics pipeline based upon the scientific workflow manager Nextflow. It was designed to help phytosanitary diagnostics of viruses and viroid pathogens in quarantine facilities. It takes small RNA-Seq fastq files as input. These can either be in raw format (currently only samples specifically prepared with the QIAGEN QIAseq miRNA library kit can be processed this way) or quality-filtered.

The pipeline can either perform blast homology searches against a virus database or/and a local copy of NCBI nr and nt databases.

**Pipeline prerequisites:** java 11 or later, Nextflow, and Docker/Singularity/Conda to suit your environment.

## Pipeline summary
![diagram pipeline](docs/images/diagram_pipeline.jpeg)

The VirReport workflow will perform the following steps by default:

- Retain reads of a given length (21-22 nt long by default) from fastq file(s) provided in index.csv file (**READPROCESSING**)  

- De novo assembly using both Velvet and SPAdes. The contigs obtained are collapsed into scaffolds using cap3. By default, only contigs > 30 bp will be retained (**DENOVO_ASSEMBLY(**)

- Run megablast homology search against either a local virus database or NCBI NT/NR databases:

- Searches against a local virus database:
  * Run megablast homology searches on de novo assembly against local virus and viroid database. Homology searches against blastn are also run in parallel for comparison with the megablast algorithm (**BLAST_NT_VIRAL_DB_CAP3**)
  * Retain top megablast hit and restrict results to virus and viroid matches. Summarise results by group all the de novo contigs matching to the same viral hit and deriving the cumulative blast coverage and percent ID for each viral hit (**FILTER_BLAST_NT_VIRAL_DB_CAP3**)
  * Align reads to top hit, derive coverage statistics, consensus sequence and VCF matching to top blast hit (**FILTER_BLAST_NT_VIRAL_DB_CAP3, COVSTATS_VIRAL_DB**)
  * Run tblastn homolgy search on predicted ORF >= 90 bp derived using getORF (**TBLASTN_VIRAL_DB**)

The pipeline can perform additional optional steps, which include:
- Searches against local NCBI NT and NR databases:
  * Retain top 5 megablast hits and restrict results to virus and viroid matches. Summarise results by group all the de novo contigs matching to the same viral hit and deriving the cumulative blast coverage and percent ID for each viral hit (**BLATN_NT_CAP3**)
  * Align reads to top hit, derive coverage statistics, consensus sequence and VCF matching to top blast hits (**COVSTATS_NT**)
  * Run blastx homolgy search on contigs >= 90 bp long for which no match was obtained in the megablast search. Summarise the blastx results and restrict to virus and viroid matches (**BLASTX**)
  
- A quality filtering step on raw fastq files (currently the workflow only processes samples prepared using QIAGEN QIAseq miRNA library kit). After performing quality filtering (**FASTQC_RAW, ADAPTER_AND_QUAL_TRIMMING, QC_POST_QUAL_TRIMMING, DERIVE_USABLE_READS**). the pipeline will also derive a qc report (QCREPORT). An RNA souce profile can be included as part of this step (**RNA_SOURCE_PROFILE, RNA_SOURCE_PROFILE_REPORT**)

- VirusDetect version 1.8 can also be run in parallel. A summary of the top virus/viroid blastn hits will be separately output (**VIRUS_DETECT, VIRUS_IDENTIFY, VIRUS_DETECT_BLASTN_SUMMARY, VIRUS_DETECT_BLASTN_SUMMARY_FILTERED**)


## Run the Pipeline
1. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Docker`](https://docs.docker.com/get-docker/), [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) or [`Conda`](https://conda.io/miniconda.html) to suit your environment.

3. Download the pipeline and test it on minimal datatests:

  This command will test your installation on a single quality filtered fastq file derived from a sample infected with citrus exocortis viroid and will run VirReport using a mock ncbi database:

  ```bash
  nextflow -c conf/test.config run eresearchqut/VirReport -profile test,{docker, singularity or conda}
  ```

  This command will test your installation on a pair of raw fastq files derived from a sample infected with citrus tristeza virus and will run VirReport using a mock viral database:

  ```bash
  nextflow -c conf/test2.config run eresearchqut/VirReport -profile test2,{docker, singularity or conda}
  ```

  Running these test datasets requires 2 cpus and 8 Gb mem and should take 5 mins to complete.

4. Run with your own data
- Provide an index.csv file.  
  Create a TAB delimited text file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the --indexfile [filename] in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,samplepath``` 

  **sampleid** will be the sample name that will be given to the files created by the pipeline  
  **samplepath** is the full path to the fastq files that the pipeline requires as starting input  

  An index_example.csv is included in the base directory:
  ```
  sampleid,samplepath
  MT212,/work/diagnostics/2021/MT212_21-22bp.fastq
  MT213,/work/diagnostics/2021/MT213_21-22bp.fastq
  ```
  
- Run the command:
  ```bash
  nextflow run eresearchqut/VirReport -profile {docker, singularity or conda} --indexfile index_example.csv
  ```  
  setting the profile parameter to one of
  ```
  docker
  singularity
  conda
  ```  
  to suit your environment.

- By default the pipeline will only retain 21-22 nt long sRNA reads for downstream analysis but you can change this range in the nextflow.config file. For instance this set of parameters will target 21-25nt long reads:
  ```
  params {
    minlen = '18'
    maxlen = '25'
  }
  ```

- Provide a database
  * By default, the pipeline is set to run homology blast searches against a local plant virus/viroid database (this is set in the nextflow.config file with parameter `--virreport_viral_db = true`. You will need to provide this database to run the pipeline. You can either provide your own or use a curated database provided at https://github.com/maelyg/PVirDB.git. Ensure you use NCBI BLAST+ makeblastdb to create the database. For instance, to set up this database, you would take the following steps:

    ```
    git clone https://github.com/maelyg/PVirDB.git
    cd PVirDB
    gunzip PVirDB_v1.fasta.gz
    makeblastdb -in PVirDB_v1.fasta -parse_seqids -dbtype nucl
    ```

    Then specify the full path to the database files including the prefix in the nextflow.config file. For example:
    ```
    params {
      blast_local_db_path = '/path_to_viral_DB/viral_DB_name'
    }
    ```  

  * If you also want to run homology searches against public NCBI databases, you need to set the parameter `virreport_ncbi` in the nextflow.config file to `true`:
    ```
    params {
      virreport_ncbi = true
    }  
    ```
    or add it in your nextflow command:  
    ```
    nextflow run eresearchqut/VirReport -profile {docker, singularity or conda} --virreport_ncbi
    ```  
    Download these locally, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. 
    Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:  
    ```
    mkdir blastDB/30112021  
    ```
    You will need to use the update_blastdb.pl script from the blast+ version used with the pipeline.  
    For example:  
    ```
    perl update_blastdb.pl --decompress nt [*]
    perl update_blastdb.pl --decompress nr [*]
    perl update_blastdb.pl taxdb
    tar -xzf taxdb.tar.gz
    ```  
    Make sure the taxdb.btd and the taxdb.bti files are present in the same directory as your blast databases.  
    Specify the path of your local NCBI blast nt and nr directories in the nextflow.config file.  
    For instance:
    ```
    params {
      blast_db_dir = '/work/hia_mt18005_db/blastDB/20220408'
    }
    ```  

  * If you want to run VirusDetect in parallel, then either set the parameter `--virusdetect` to `true` in your config file or specify it in your nextflow run command.
    You will need to download the VirusDetect viral database files at http://bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v248. And specify the path to the VirusDetect viral database in the nextflow.config file (using the `--virusdetect_db_path` parameter).  For example:
  ```
  virusdetect_db_path = '/home/gauthiem/code/VirusDetect_v1.8/databases/vrl_plant'
  ```

- If you want to run the initial qualityfilter step on raw fastq files, you will need to specify in the nextflow.config file the directory which holds the required bowtie indices (using the **--bowtie_db_dir** parameter) to: 1) filter non-informative reads (using the blacklist bowtie indices for the DERIVE_USABLE_READS process) and 2) derive the origin of the filtered reads obtained (optional RNA_SOURCE_PROFILE process). Examples of fasta files are available at https://github.com/maelyg/bowtie_indices.git and bowtie indices can be built from these using Bowtie.
For instance to derive the bowtie indices for the blacklist, run the following command:

  ```
  git clone https://github.com/maelyg/bowtie_indices.git
  gunzip blacklist_v2.fasta.gz
  bowtie-build -f blacklist_v2.fasta blasklist
  ```

  The location of the directory in which the bowtie indices are located will need to be specified in the nextflow.config file:
  ```
  params {
    bowtie_db_dir = '/path_to_bowtie_idx_directory'
  }
  ```

  If you are interested to derive the RNA profile of your fastq files you will need to specify:
  ```
  params {
    rna_source_profile = true
  }
  ```
  And build the other indices from the fasta files included in https://github.com/maelyg/bowtie_indices.git (i.e. rRNA, plant_tRNA, plant_noncoding, plant_pt_mt_other_genes, artefacts, plant_miRNA, virus).

- Additional optional parameters available include:
  ```     
    --merge-lane: if several fastq files are provided per sample, these will be collapsed together before performing downstream analyses
   
    --dedup: first umi_tools will be used to extract UMI informations from fast file headers and after alignment, umi_tools will be used to identify and remove duplicate reads.
    --spadesmem specifies the memory usage available for SPAdes (by default 32)
      
    --cap3_len: specifies the minimal length of contigs to retain after CAP3 scaffolding (by default 30)

    --blastn_method: The blastn homology search can be specified as blastn instead of megablast using --blastn_method blastn
      
    --blastn_evalue and --blastp_evalue: specifies the evalue parameter to use during blast analyses (by deafult 0.0001)
     
    --targets: A text file with the taxonomy of the viruses/virioids of interest can be provided and only these will be retained in the megablast summary results derived at the COVSTATS_NT step using NCBI NT database.
    --blastx: if set to false, it will skip the blastx step performed when searching against local NCBI NT and NR databases 
      
    --orf_minsize and --orf_circ_minsize: correspond to the minimal open reading frames getorf retains that will be used in the tblastn homology search against the virus database (by default 90 bp) 
    --contamination_detection and contamination_detection_viral_db: Run cross-sample contamination prediction (CONTAMINATION_DETECTION and/or CONTAMINATION_DETECTION_VIRAL_DB) 
  ```

  To enable these options, they can either be included in the nextflow run command: 
  ```
  nextflow run eresearchqut/VirReport -profile {singularity, docker or conda} --indexfile index_example.csv --contamination_detection --virusdetect
  ```
  or the parameter in the nextflow.config file can be udpated. For instance:
  ```
  params {
    virusdetect = true
    contamination_detection = true
  }
  ```

## Outputs
The folders are structured as follows:
```
results/
├── 00_quality_filtering
│   └── sample_name
│   │   ├── sample_name_18-25nt_cutadapt.log
│   │   ├── sample_name_fastqc.html
│   │   ├── sample_name_fastqc.zip
│   │   ├── sample_name_21-22nt_cutadapt.log
│   │   ├── sample_name_21-22nt.fastq.gz
│   │   ├── sample_name_24nt_cutadapt.log
│   │   ├── sample_name_blacklist_filter.log
│   │   ├── sample_name_fastp.html
│   │   ├── sample_name_fastp.json
│   │   ├── sample_name_qual_filtering_cutadapt.log
│   │   ├── sample_name_quality_trimmed_fastqc.html
│   │   ├── sample_name_quality_trimmed_fastqc.zip
│   │   ├── sample_name_quality_trimmed.fastq.gz
│   │   ├── sample_name_read_length_dist.pdf
│   │   ├── sample_name_read_length_dist.txt
│   │   ├── sample_name_truseq_adapter_cutadapt.log
│   │   └── sample_name_umi_tools.log
│   └── qc_report
│       ├── read_origin_counts.txt
│       ├── read_origin_detailed_pc.txt
│       ├── read_origin_pc_summary.txt
│       ├── read_origin_pc_summary.txt
│       ├── run_qc_report.txt
│       └── run_read_size_distribution.pdf
├── 01_VirReport
│   └── sample_name
│   │   └── alignments
│   │   │   └──NT
│   │   │   │   ├── sample_name_21-22nt_all_targets_with_scores.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_bowtie_log.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.consensus.fasta
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.dedup.bam
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.dedup.bam.bai
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.fa
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.fa.fai
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm.bcf
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm.bcf.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm_flt_indels.bcf
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm_flt_indels.bcf.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_picard_metrics.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_sequence_variants.vcf.gz
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_sequence_variants.vcf.gz.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_umi_tools.log
│   │   │   │   ├── sample_name_21-22nt_top_scoring_targets.txt
│   │   │   │   └── sample_name_21-22nt_top_scoring_targets_with_cov_stats.txt
│   │   │   └──viral_db
│   │   │   │   ├── sample_name_21-22nt_all_targets_with_scores.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_bowtie_log.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.consensus.fasta
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.dedup.bam
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.dedup.bam.bai
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.fa
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name.fa.fai
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm.bcf
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm.bcf.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm_flt_indels.bcf
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_norm_flt_indels.bcf.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_picard_metrics.txt
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_sequence_variants.vcf.gz
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_sequence_variants.vcf.gz.csi
│   │   │   │   ├── sample_name_21-22nt_GenBankID_virus_name_umi_tools.log
│   │   │   │   └── sample_name_21-22nt_top_scoring_targets_with_cov_stats_viraldb.txt
│   │   └── assembly
│   │   │   ├── sample_name_cap3_21-22nt.fasta
│   │   │   ├── sample_name_spades_assembly_21-22nt.fasta
│   │   │   ├── sample_name_spades_log
│   │   │   ├── sample_name_velvet_assembly_21-22nt.fasta
│   │   │   └── sample_name_velvet_log
│   │   └──blastn
│   │   │   └── NT
│   │   │   │   ├── sample_name_cap3_21-22nt_blastn_vs_NT.bls
│   │   │   │   ├── sample_name_cap3_21-22nt_blastn_vs_NT_top5Hits.txt
│   │   │   │   ├── sample_name_cap3_21-22nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt
│   │   │   │   ├── sample_name_cap3_21-22nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
│   │   │   │   └── summary_sample_name_cap3_21-22nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt
│   │   │   └── viral_db
│   │   │       ├── sample_name_cap3_21-22nt_blastn_vs_viral_db.bls
│   │   │       ├── sample_name_cap3_21-22nt_megablast_vs_viral_db.bls
│   │   │       ├── summary_sample_name_cap3_21-22nt_blastn_vs_viral_db.bls_filtered.txt
│   │   │       ├── summary_sample_name_cap3_21-22nt_blastn_vs_viral_db.bls_viruses_viroids_ICTV.txt
│   │   │       ├── summary_sample_name_cap3_21-22nt_megablast_vs_viral_db.bls_filtered.txt
│   │   │       └── summary_sample_name_cap3_21-22nt_megablast_vs_viral_db.bls_viruses_viroids_ICTV.txt
│   │   ├── blastx
│   │   │   └── NT
│   │   │       ├── sample_name_cap3_21-22nt_blastx_vs_NT.bls
│   │   │       ├── sample_name_cap3_21-22nt_blastx_vs_NT_top5Hits.txt
│   │   │       ├── sample_name_cap3_21-22nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt
│   │   │       └── summary_sample_name_cap3_21-22nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt
│   │   └── tblastn
│   │       └── viral_db
│   │           ├── sample_name_cap3_21-22nt_getorf.all.fasta
│   │           ├── sample_name_cap3_21-22nt_getorf.all_tblastn_vs_viral_db_out.bls
│   │           └── sample_name_cap3_21-22nt_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids_final.txt
│   └── Summary
│       ├── run_top_scoring_targets_with_cov_stats_with_cont_flag_FPKM_0.01_21-22nt.txt
│       └── run_top_scoring_targets_with_cov_stats_with_cont_flag_FPKM_0.01_21-22nt_viral_db.txt
└── 02_VirusDetect
    └── sample_name
    │   ├── blastn.reference.fa
    │   ├── blastn_references
    │   ├── blastx.reference.fa
    │   ├── blastx_references
    │   ├── contig_sequences.blastn.fa
    │   ├── contig_sequences.blastx.fa
    │   ├── contig_sequences.fa
    │   ├── contig_sequences.undetermined.fa
    │   ├── sample_name_21-22nt.blastn.html
    │   ├── sample_name_21-22nt.blastn.sam
    │   ├── sample_name_21-22nt.blastn_spp.txt
    │   ├── sample_name_21-22nt.blastn.summary.filtered.txt
    │   ├── sample_name_21-22nt.blastn.summary.txt
    │   ├── sample_name_21-22nt.blastn.txt
    │   ├── sample_name_21-22nt.blastx.html
    │   ├── sample_name_21-22nt.blastx.sam
    │   ├── sample_name_21-22nt.blastx.summary.txt
    │   └── sample_name_21-22nt.blastx.txt
    └── Summary
        ├── run_summary_top_scoring_targets_virusdetect_21-22nt_filtered.txt
        └── run_summary_virusdetect_21-22nt.txt
```     
- 00_quality_filtering/sample_name: this folder will output FASTQC of raw or filtered fastq files, cutadapt and umi_tools log files, a quality_trimmed.fastq.gz file and by default a fastq.gz file including reads only matching the size specified in the nextflow.config file

- 00_quality_filter/qc_report: this folder contains QC summaries for all samples included in the index_example.csv file

- 01_VirReport/sample_name/assembly: fasta file which includes the assembled contigs before and after CAP3. The sample_name_cap3_21-22nt.fasta will be used  for homology searches in the next steps.

- 01_VirReport/sample_name/blastn/NT: this folder contains all megablast results, filtered results limited to only viruses and viroid top 5 hit matches and the final BlastTools.jar summary output. 
For example: sample_name_cap3_21-22nt_blastn_vs_NT.bls, summary_sample_name_cap3_21-22nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt

- 01_VirReport/sample_name/blastn/viral_db: Analysis using the viral db will be saved in this folder
sample_name_cap3_21-22nt_megablast_vs_viral_db.bls, summary_sample_name_cap3_21-22nt_megablast_vs_viral_db.bls_viruses_viroids_ICTV.txt

- 01_VirReport/sample_name/blastx/NT: this folder contains the blastx outputs. 

- 01_VirReport/sample_name/alignments: this folder contains a filtered blast summary with various coverage statistics for each virus and viroid hit, and associated consensus fasta file and vcf file. The results for each sample are saved in a separate folder. 
For example: sample_name_21-22nt_top_scoring_targets_with_cov_stats.txt, sample_name_21-22nt_GenBankID_virus_name.consensus.fasta 

- 01_VirReport/Summary: this folder contains a summary of results for all samples included in the index.csv file. The summay table includes a cross-contamination prediction flag. 
For example: run_top_scoring_targets_with_cov_stats_with_cont_flag_21-22nt_0.01.txt

- 02_VirusDetect/sample_name: this folder includes a results folder with blastn and blastx summary. 
For example: sample_name_21-22nt.blastn.summary.txt and sample_name_21-22nt.blastx.summary.txt

 
## Credits
Roberto Barrero, 14/03/2019  
Desmond Schmidt, 2/7/2019  
Converted to Nextflow by Craig Windell, 11/2020  
Modified by Maely Gauthier 12/2021

# Citations
If you use nf-core/rnaseq for your analysis, please cite it using the following doi: 10.3390/biology11020263
