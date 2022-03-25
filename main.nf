#!/usr/bin/env nextflow
/*
VirReport workflow
Roberto Barrero, 14/03/2019
Desmond Schmidt, 2/7/2019
Converted to Nextflow by Craig Windell, 11/2020
Modified by Maely Gauthier, 12/2021
*/

import java.util.*;
import java.util.stream.IntStream;
import java.util.stream.Collectors;

def helpMessage () {
    log.info """

    VirReport workflow
    Roberto Barrero, 14/03/2019
    Desmond Schmidt, 2/7/2019
    Converted to Nextflow by Craig Windell, 11/2020
    Modified by Maely Gauthier, 12/2021

    Usage:

    Run the command
    nextflow run eresearchqut/virreport -profile ...

    Mandatory arguments:
       -profile '[docker, singularity, conda]'      Profile to use. Choose docker, or singularity, or conda

    Optional arguments:
      --indexfile '[path/to/file]'                  Path to the csv file that contains the list of
                                                    samples to be analysed by this pipeline.
                                                    'index.csv'
      Contents of indexfile csv:
        sampleid,samplepath,params.minlen,params.maxlen
        MT019,MT019_sRNA.fastq,21,22

      --blast_local_nt_db '[path/to/dir]'           Path to the local blast NT database files. Required if --blastlocaldb option is specified
                                                    '/work/hia_mt18005/databases/sequences/PVirDB_20210330'

      --blast_nt_db '[path/to/dir]'                 Path to the blast NT database files
                                                    '/work/eresearch_bio/reference/blastDB/nt'

      --blast_nr_db '[path/to/dir]'                 Path to the blast NR database files
                                                    '/work/eresearch_bio/reference/blastDB/nr'

      --blastlocaldb                                Run blastn and megablast homology search on velvet de novo assembly against local virus and viroid database
                                                    [False]

      --blastn_evalue '[value]'                     Blastn evalue.
                                                    '0.0001'

      --blastn_method ['blastn/megablast']          Run blastn homology search on velvet de novo assembly againts NCBI NT
                                                    [default megablast]

      --blastp [True/False]                         Predict ORF from de novo assembled contigs and run blastP againts NCBI NR
                                                    [False]

      --blastp_evalue '[value]'                     Blastp's evalue. Required if --blastp option is specified
                                                    '0.0001'
                    
      --cap3_len '[value]'                          Trim value used in the CAP3 step.
                                                    '20'

      --contamination_detection [True/False]        Run false positive prediction due to cross-sample contamination
                                                    [False]

      --contamination_detection_method '[value]'    Either use FPKM or RPM for cross-contamination detection 

      
      --contamination_flag '[value]'                Threshold value to predict false positives due to cross-sample contamination. 
                                                    Required if --contamination_detection option is specified
                                                    '0.01'
      --ictvinfo '[path/to/dir]'                    Path to ICTV info file. Required if --blastlocaldb option is specified
                                                    ['ICTV_taxonomy_MinIdentity_Species.tsv']
      
      --orf_minsize '[value]'                       The value of minsize for getorf
                                                    '150'

      --orf_circ_minsize '[value]'                  The value of minsize for getorf -circular
                                                    '150'

      --spades [True/False]                         Run SPAdes 3.14 de novo assembler and perform blastn homology analysis on the derived de novo contigs
                                                    [False]

      --spadeskmer '[value]'                        K-mer range to use when running SPAdes. Required if --spades option is specified
                                                    ['K9-21']

      --targets [True/False]                        Filter the blastn results to viruses/viroids of interest
                                                    [False]

      --targets_file '[path/to/dir]'                File specifying the name of the viruses/viroids of interest to filter from the blast results output
                                                    ['Targetted_Viruses_Viroids.txt']
    
      --qualityfilter [True/False]                  Perform initial quality filtering of fastq file
                                                    [False]
      --minlen '[value]'                            Minimum read length to extract
      ['21']
      --maxlen '[value]'                            Maximum read length to extract
      ['22']
    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

blastn_db_name = "${params.blast_db_dir}/nt"
blastp_db_name = "${params.blast_db_dir}/nr"
negative_seqid_list = "${params.blast_db_dir}/negative_list.txt"
blast_local_db_name = file(params.blast_local_db_path).name
blast_local_db_dir = file(params.blast_local_db_path).parent
blacklist_db_name = file(params.blacklist_db_path).name
blacklist_db_dir = file(params.blacklist_db_path).parent

switch (workflow.containerEngine) {
    case "docker":
        bindOptions = "-v ${params.blast_db_dir}:${params.blast_db_dir} -v ${blast_local_db_dir}:${blast_local_db_dir} -v ${blacklist_db_dir}:${blacklist_db_dir}"
        break;
    case "singularity":
        bindOptions = "-B ${blast_local_db_dir} -B ${params.blast_db_dir} -B ${blacklist_db_dir}"
        break;
    default:
        bindOptions = ""
}

if (params.indexfile) {
    if (params.qualityfilter) {
        Channel
          .fromPath(params.indexfile, checkIfExists: true)
          .splitCsv(header:true)
          .map{ row-> tuple(row.samplepath.substring(row.samplepath.lastIndexOf('/') +1), file(row.samplepath)) }
          .into{ fastqcraw_ch; filter_n_cov_ch; contamination_detection_ch}
        Channel
          .fromPath(params.indexfile, checkIfExists: true)
          .splitCsv(header:true)
          .map{ row-> tuple(row.sampleid, file(row.samplepath)) }
          .groupTuple()
          .set{ samples_ch }
    }
    else {
        Channel
          .fromPath(params.indexfile, checkIfExists: true)
          .splitCsv(header:true)
          .map{ row-> tuple(row.sampleid), file(row.samplepath) }
          .into{ read_size_selection_ch; filter_n_cov_ch; contamination_detection_ch}
        }
    }
    else { exit 1, "Input samplesheet file not specified!" }

if (params.qualityfilter) {
    process fastqc_raw {
        memory '60 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from fastqcraw_ch
        
        output:
        file "*_fastqc.{zip,html}"

        script:
        """
        fastqc --quiet --threads ${task.cpus} ${fastqfile}
        """
    }

    process merge_lanes {
        memory '15 GB'
        cpus 2
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(samplepath) from samples_ch
        
        output:
        //file "${sampleid}_R1.merged.fastq.gz"
        tuple val(sampleid), file("${sampleid}_R1.merged.fastq.gz") into umitools_ch

        script:
        if (params.merge_lane) {
          samplepathList = samplepath.collect{it.toString()}
          if (samplepathList.size > 1 ) {
          """
          cat ${samplepath} > ${sampleid}_R1.merged.fastq.gz
          """
          }
        } else {
          """
          ln ${samplepath} ${sampleid}_R1.merged.fastq.gz
          """
        }
    }

    process umiclean_cutadapt {
        memory '80 GB'
        cpus 8
        time '2h'
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from umitools_ch

        output:
        file "${sampleid}_umi_tools.log"
        file "${sampleid}_truseq_adapter_cutadapt.log"
        file "${sampleid}_qual_filtering_cutadapt.log"
        file "${sampleid}_quality_trimmed.fastq"
        file "${sampleid}_umi_tools.log" into umi_tools_results
        file "${sampleid}_qual_filtering_cutadapt.log" into cutadapt_qual_filt_results

        tuple val(sampleid), file("${sampleid}_quality_trimmed.fastq") into fastqc_filtered_ch, fastp_report_ch, read_size_distribution_ch, rna_source_distribution_ch, derive_usable_reads_ch, read_size_selection_ch

        script:
        """
        cutadapt -j ${task.cpus} \
                --no-indels \
                -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;min_overlap=12" \
                -g "ACACTCTTTCCCTACACGACGCTCTTCCGATCT;min_overlap=9" \
                --times 2 \
                -o ${sampleid}_trimmed.fastq \
                ${fastqfile} > ${sampleid}_truseq_adapter_cutadapt.log

        umi_tools extract --extract-method=regex \
                            --bc-pattern=".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})\$" \
                            -I ${sampleid}_trimmed.fastq \
                            -S ${sampleid}_umi_cleaned.fastq > ${sampleid}_umi_tools.log
        
        cutadapt -j ${task.cpus} \
                --trim-n --max-n 0 -m 5 -q 30 \
                -o ${sampleid}_quality_trimmed.fastq \
                ${sampleid}_umi_cleaned.fastq > ${sampleid}_qual_filtering_cutadapt.log
        """
    }

    process fastqc_filtered {
        memory '60 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from fastqc_filtered_ch
        
        output:
        file "*_fastqc.{zip,html}" into fastqc_qual_results

        script:
        """
        fastqc --quiet --threads ${task.cpus} ${fastqfile}
        """
    }

    process fastp_report {
        memory '60 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from fastp_report_ch
        
        output:
        file "${sampleid}_fastp.json"
        file "${sampleid}_fastp.html"
        file "${sampleid}_fastp.json" into fastp_results

        script:
        """
        fastp --in1=${fastqfile} --out1=${sampleid}_fastp_trimmed.fastq \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --json=${sampleid}_fastp.json \
            --html=${sampleid}_fastp.html \
            --thread=${task.cpus}
        """
    }

    process read_size_distribution {
        memory '60 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from read_size_distribution_ch
  
        output:
        file "${sampleid}_read_length_dist.pdf"
        file "${sampleid}_read_length_dist.txt"
        file "${sampleid}_read_length_dist.txt" into read_length_dist_results

  
        script:
        """
        #derive distribution for quality filtered reads > 5 bp long
        perl ${projectDir}/bin/fastq2fasta.pl ${fastqfile} > ${sampleid}_trimmed.fasta
        python ${projectDir}/bin/read_length_dist.py --input ${sampleid}_trimmed.fasta
        mv ${sampleid}_trimmed.fasta_read_length_dist.txt ${sampleid}_read_length_dist.txt
        mv ${sampleid}_trimmed.fasta_read_length_dist.pdf ${sampleid}_read_length_dist.pdf
        """
    }
    process rna_source_distribution {
        memory '60 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from rna_source_distribution_ch

        output:
        file "${sampleid}_bowtie.log"
        file "${sampleid}_UniVec_cleaned_sRNA.fq"
        file "${sampleid}_final_unaligned_sRNA.fq"
        file "${sampleid}_bowtie.log" into rna_source_bowtie_results

        
        script:
        """
        #derive distribution for quality filtered reads > 5 bp long
        echo ${sampleid} > ${sampleid}_bowtie.log;

        echo 5S rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_5S_rRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/5S_rRNA \
            ${fastqfile} \
            ${sampleid}_5S_rRNA_match 2>>${sampleid}_bowtie.log;


        echo nc SSU and LSU rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_nc_rRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/rRNA \
            ${sampleid}_5S_rRNA_cleaned_sRNA.fq \
            ${sampleid}_rRNA_match 2>>${sampleid}_bowtie.log;

        echo mt rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_mt_rRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_mt_rRNA_genes \
            ${sampleid}_nc_rRNA_cleaned_sRNA.fq \
            ${sampleid}_mt_rRNA_match 2>>${sampleid}_bowtie.log;

        echo pt rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_pt_rRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_pt_rRNA_genes \
            ${sampleid}_mt_rRNA_cleaned_sRNA.fq \
            ${sampleid}_pt_rRNA_match 2>>${sampleid}_bowtie.log;

        echo mt other genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_mt_other_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_mt_other_genes \
            ${sampleid}_pt_rRNA_cleaned_sRNA.fq \
            ${sampleid}_mt_other_match 2>>${sampleid}_bowtie.log;

        echo pt other genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_pt_other_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_pt_other_genes \
            ${sampleid}_mt_other_cleaned_sRNA.fq \
            ${sampleid}_pt_other_match 2>>${sampleid}_bowtie.log;

        echo plant miRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_miRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_miRNA \
            ${sampleid}_pt_other_cleaned_sRNA.fq \
            ${sampleid}_plant_miRNA_match 2>>${sampleid}_bowtie.log;

        echo other miRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_other_miRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/miRBase \
            ${sampleid}_plant_miRNA_cleaned_sRNA.fq \
            ${sampleid}_other_miRNA_match 2>>${sampleid}_bowtie.log;

        echo tRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_tRNA_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_tRNA \
            ${sampleid}_other_miRNA_cleaned_sRNA.fq \
            ${sampleid}_tRNA_match 2>>${sampleid}_bowtie.log;

        echo plant nc genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_nc_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_nc_genes \
            ${sampleid}_plant_tRNA_cleaned_sRNA.fq \
            ${sampleid}_plant_nc_match 2>>${sampleid}_bowtie.log;

        echo plant transposons alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_transposon_cleaned_sRNA.fq \
            -x ${projectDir}/databases/plant_transposons \
            ${sampleid}_plant_nc_cleaned_sRNA.fq \
            ${sampleid}_transposon_match 2>>${sampleid}_bowtie.log;

        echo PhiX alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_PhiX_cleaned_sRNA.fq \
            -x ${projectDir}/databases/phiX174 \
            ${sampleid}_transposon_cleaned_sRNA.fq \
            ${sampleid}_PhiX_sRNA_match 2>>${sampleid}_bowtie.log;

        echo Vector alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_UniVec_cleaned_sRNA.fq \
            -x ${projectDir}/databases/UniVec \
            ${sampleid}_PhiX_cleaned_sRNA.fq \
            ${sampleid}_UniVec_match 2>>${sampleid}_bowtie.log;

        echo Virus and viroid alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_final_unaligned_sRNA.fq \
            -x ${projectDir}/databases/virus \
            ${sampleid}_UniVec_cleaned_sRNA.fq \
            ${sampleid}_viral_match 2>>${sampleid}_bowtie.log;
        """
    }

    process derive_usable_reads {
        memory '80 GB'
        cpus 4
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from derive_usable_reads_ch
        
        output:
        file "${sampleid}*_cutadapt.log"
        file "${sampleid}_blacklist_filter.log"
        file "${sampleid}_${params.minlen}-${params.maxlen}nt.fastq"
        tuple val(sampleid), file(fastqfile), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into velvet_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into spades_ch
        file ("*_18-25nt_cutadapt.log") into cutadapt_18_25nt_results
        file ("*_21-22nt_cutadapt.log") into cutadapt_21_22nt_results
        file ("*_24nt_cutadapt.log") into cutadapt_24nt_results
        file ("*_blacklist_filter.log") into bowtie_usable_read_cat_results

        script:
        """
        bowtie -q -v 1 \
            -k 1 --un ${sampleid}_cleaned.fastq -p ${task.cpus} \
            -x ${projectDir}/databases/blacklist \
            ${fastqfile} \
            ${sampleid}_blacklist_match 2>${sampleid}_blacklist_filter.log


        cutadapt -j ${task.cpus} -m 18 -M 25 -o ${sampleid}_18-25nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_18-25nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 21 -M 22 -o ${sampleid}_21-22nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_21-22nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 24 -M 24 -o ${sampleid}_24nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_24nt_cutadapt.log
        if [[ ${params.minlen} != 21 ]] || [[ ${params.maxlen} != 22 ]]; then
            cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log
        fi
        """
    }

    process qcreport {
        publishDir "${params.outdir}/01_QC_report", mode: 'link'

        input:
        file ('*') from cutadapt_qual_filt_results.collect().ifEmpty([])
        file ('*') from fastp_results.collect().ifEmpty([])
        file ("*") from cutadapt_18_25nt_results.collect().ifEmpty([])
        file ("*") from cutadapt_21_22nt_results.collect().ifEmpty([])
        file ("*") from cutadapt_24nt_results.collect().ifEmpty([])
        file ('*') from bowtie_usable_read_cat_results.collect().ifEmpty([])
        file ('*') from read_length_dist_results.collect().ifEmpty([])
        file ('*') from umi_tools_results.collect().ifEmpty([]) 
        file ('*') from rna_source_bowtie_results.collect().ifEmpty([]) 

        output:
        file "run_report.txt"
        file "run_read_size_distribution.pdf"
        file "read_origin_pc_summary.txt"
        file "read_RNA_source.pdf"
        file "read_origin_detailed_pc.txt"


        script:
        """
        qc_report.py
        grouped_bar_chart.py
        rna_source_summary.py
        """
    }
} 
// If user does not specify qualityfilter parameter, then only read size selection (using the minlen and maxlen params specified in the nextflow.config file) will be performed on the fastq file specified in the index file
else {
    process readprocessing {
        tag "$sampleid"
        publishDir "${params.outdir}/01_read_size_selection", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from read_size_selection_ch

        output:
        file "${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log"
        file "${sampleid}_${params.minlen}-${params.maxlen}nt.fastq"
        tuple val(sampleid), file(fastqfile), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into velvet_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into spades_ch

        script:
        """
        cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq ${fastqfile} > ${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log
        """
    }
}

process velvet {
    publishDir "${params.outdir}/02_velvet/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size) from velvet_ch

    output:
    file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15/*"
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta") into cap3_ch

    script:
    """
    #run velvet de novo assembler
    echo 'Starting velvet de novo assembly';
    velveth ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15 15 -short -fastq ${fastq_filt_by_size}
    velvetg ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15 -exp_cov 2

    #edit contigs name and rename velvet assembly
    sed 's/>/>velvet_/' ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15/contigs.fa > ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta
    """
}

process cap3 {
    label "local"
    publishDir "${params.outdir}/03_cap3/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(scaffolds_fasta) from cap3_ch

    output:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta") into blastn_nt_velvet_ch
    tuple val(sampleid), file("${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta") into blast_nt_localdb_velvet_ch, getorf_ch

    script:
    """
    #collapse velvet contigs
    cap3 ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta -s 300 -j 31 -i 30 -p 90 -o 16
    cat ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta.cap.singlets ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta.cap.contigs > ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt.fasta
    extract_seqs_rename.py ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt.fasta ${params.cap3_len} \
                             | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" \
                             | sed 's/|>/ |/' | awk '{print \$1}'\
                             > ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta
    """
}

process blastn_nt_velvet {
    label "blastn_mem"
    publishDir "${params.outdir}/04_blastn/${sampleid}", mode: 'link'
    tag "$sampleid"
    containerOptions "${bindOptions}"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta") from blastn_nt_velvet_ch

    output:
    file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls"
    file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt"
    file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt") into blastTools_blastn_velvet_ch
    
    script:
    def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
    """
    #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
    cp ${params.blast_db_dir}/taxdb.btd .
    cp ${params.blast_db_dir}/taxdb.bti .
    blastn ${blast_task_param} \
        -query ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta \
        -db ${blastn_db_name} \
        -negative_seqidlist ${negative_seqid_list} \
        -out ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls \
        -evalue 0.0001 \
        -num_threads 4 \
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
        -max_target_seqs 50

    grep ">" ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta | sed 's/>//' > ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.ids
    
    #fetch top blastn hits
    for i in `cat ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.ids`; do
        grep \$i ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls | head -n5 >> ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt;
    done
    
    grep -i "Virus" ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
    grep -i "Viroid" ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt >> ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
    cat ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt
    cut -f3,26 ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt | sort | uniq > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
    """
}

process BlastTools_blastn_velvet {
    label "medium_mem"
    publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(top5Hits_final), file(taxonomy) from blastTools_blastn_velvet_ch

    output:
    file "summary_${top5Hits_final}"
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("summary_${top5Hits_final}"), file(taxonomy) into blastTools_results_ch

    script:
    """
    java -jar ${projectDir}/bin/BlastTools.jar -t blastn ${top5Hits_final}
    """
}

if (params.blastlocaldb) {
    process blast_nt_localdb_velvet {
        label "blastn_mem"
        publishDir "${params.outdir}/04_blastn/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file("${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta") from blast_nt_localdb_velvet_ch
        
        output:
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls"
        tuple val(sampleid), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls") into filter_blast_nt_localdb_velvet_ch

        script:
        """
        #1. blastn search
        blastn -task blastn \
            -query ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50

        #2. megablast search
        blastn -query ${sampleid}_velvet_cap3_${params.minlen}-${params.maxlen}nt_rename.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50
        """
    }

    process filter_blast_nt_localdb_velvet {
        label "local"
        publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls") from filter_blast_nt_localdb_velvet_ch

        output:
        //file "summary_${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls_ENDEMIC_viruses_viroids_ICTV.txt"
        //file "summary_${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls_REGULATED_viruses_viroids_ICTV.txt"
        //file "summary_${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls_ENDEMIC_viruses_viroids_ICTV.txt"
        file "summary_${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls_viruses_viroids_ICTV.txt"
        file "summary_${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls_viruses_viroids_ICTV.txt"
        
        script:
        """
        #retain 1st blast hit
        for var in ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls;
            do 
                cat \${var} | awk '{print \$1}' | sort | uniq > \${var}.top1.ids
                for i in `cat \${var}.top1.ids`; do echo "fetching top hits..." \$i; grep \$i \${var} | head -1 >> \${var}.top1Hits.txt ; done
                cat \${var}.top1Hits.txt | sed 's/ /_/g' > \${var}.txt

                #summarise the blast files
                java -jar ${projectDir}/bin/BlastTools.jar -t blastn \${var}.txt

                #only retain hits to plant viruses (regulated/edemic/LandPlant)
                cat summary_\${var}.txt | grep "regulated\\|endemic\\|higher_plant_viruses" > summary_\${var}_filtered.txt

                if ! [ -s summary_\${var}_filtered.txt ]; then
                    touch summary_\${var}_viruses_viroids.txt;
                    touch summary_\${var}_viruses_viroids_ICTV.txt;
                else
                    #fetch unique virus/viroid species name from Blast summary reports
                    cat summary_\${var}_filtered.txt | awk '{print \$7}' | awk -F "|" '{print \$4}'| sort | uniq | sed 's/Species://' > \${var}_uniq.ids

                    #retrieve the best hit for each unique virus/viroid species name by selecting longest alignment (column 3) and highest genome coverage (column 5)
                    for id in `cat \${var}_uniq.ids`;
                        do
                            grep \${id} summary_\${var}.txt | sort -k3,3nr -k5,5nr | head -1 >> \${var}_filtered.txt
                        done

                    #print the header of the inital summary_blastn file
                    cat summary_\${var}.txt | head -1 > header

                    #report 1
                    cat header \${var}_filtered.txt > summary_\${var}_viruses_viroids.txt
                    
                    #fetch genus names of identified hits
                    awk '{print \$7}' summary_\${var}_viruses_viroids.txt | awk -F "|" '{print \$4}' | sed 's/Species://' | sed 1d > wanted.names
                
                    #add species to report
                    paste wanted.names \${var}_filtered.txt | sort > summary_\${var}_viruses_viroids.MOD

                    #fecth ICTV information
                    grep -w -F -f wanted.names ${projectDir}/bin/${params.ictvinfo} | sort > wanted.ICTV

                    #join reports with ICTV information
                    #join -a 1 -1 1 -2 1 summary_\${var}_viruses_viroids.MOD  wanted.ICTV | tr ' ' '\\t' | awk '\$4>=70' >  summary_\${var}_viruses_viroids_ICTV
                    join -a1 -1 1 -2 1 summary_\${var}_viruses_viroids.MOD  wanted.ICTV | tr ' ' '\\t' >  summary_\${var}_viruses_viroids_ICTV

                    #report 2
                    awk '{print "Species" "\\t" \$0 "\\t" "ICTV_information"}' header > header2
                    cat header2 summary_\${var}_viruses_viroids_ICTV | awk -F"\\t" '\$1!=""&&\$2!=""&&\$3!=""' > summary_\${var}_viruses_viroids_ICTV.txt
                fi
            done
        """
    }
}
process filter_n_cov {
    tag "$sampleid"
    publishDir "${params.outdir}/07_filternstats/${sampleid}", mode: 'link'
    containerOptions "${bindOptions}"
    
    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(samplefile), file(taxonomy) from blastTools_results_ch

    output:
    file "${sampleid}_${params.minlen}-${params.maxlen}*"
    file("${sampleid}_${params.minlen}-${params.maxlen}nt_top_scoring_targets_with_cov_stats.txt") into contamination_flag
    
    script:
    """
    if [[ ${params.targets} == true ]]; then
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --cov --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name} --targets --targetspath ${projectDir}/bin/${params.targets_file}
    else
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --cov --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name}
    fi
    """
}

if (params.contamination_detection) {
    process contamination_detection {
        label "local"
        publishDir "${params.outdir}/08_summary", mode: 'link'
        
        input:
        tuple val(sampleid), file(fastqfile) from contamination_detection_ch
        file ('*') from contamination_flag.collect().ifEmpty([])

        output:
        file "run_top_scoring_targets_with_cov_stats_with_cont_flag*.txt"

        script:
        """
        flag_contamination.py --read_size ${params.minlen}-${params.maxlen}nt --threshold ${params.contamination_flag} --method ${params.contamination_detection_method}
        """
    }
}

if (params.blastp) {
    process getorf {
        label 'local'
        publishDir "${params.outdir}/06_blastp/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(velvet_cap3_rename_fasta) from getorf_ch
        
        output:
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta"
        tuple val(sampleid), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids") into blastp_ch
        
        script:
        """
        getorf -sequence ${velvet_cap3_rename_fasta} -outseq ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta  -minsize ${params.orf_minsize}
        cat ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids
        getorf -sequence ${velvet_cap3_rename_fasta} -circular True -outseq ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta.ids
        """
    }

    process blastp {
        label "xlarge"
        publishDir "${params.outdir}/06_blastp/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fasta), file(fasta_ids) from blastp_ch
        
        output:
        file "${fasta.baseName}_blastp_vs_NR_out.bls"
        tuple val(sampleid), file("${fasta.baseName}_blastp_vs_NR_out.wanted.ids"), file("${fasta.baseName}_blastp_vs_NR_out.bls") into blastpdbcmd_ch
        
        script:
        """         
        blastp -query ${fasta} \
            -db ${blastp_db_name} \
            -evalue ${params.blastp_evalue} \
            -out ${fasta.baseName}_blastp_vs_NR_out.bls \
            -num_threads ${task.cpus} \
            -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid'

        cat ${fasta.baseName}_blastp_vs_NR_out.bls | awk '{print \$2}' | cut -f2 -d '|' | sort | uniq > ${fasta.baseName}_blastp_vs_NR_out.wanted.ids
        """
    }

    process blastpdbcmd {
        label "medium_mem"
        publishDir "${params.outdir}/06_blastp/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(blastp_nr_bls_ids), file(blastp_nr_bls) from blastpdbcmd_ch

        output:
        file "${blastp_nr_bls.baseName}_virus_viroid.txt"
        tuple val(sampleid), file("${blastp_nr_bls.baseName}_virus_viroid.txt") into BlastToolsp_ch

        script:
        """
        blastdbcmd  -db ${blastp_db_name} \
                    -dbtype prot \
                    -entry_batch ${blastp_nr_bls_ids} > ${blastp_nr_bls_ids.baseName}.out

        grep ">" ${blastp_nr_bls_ids.baseName}.out | sed 's/>//' > ${blastp_nr_bls_ids.baseName}.out.header || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header | sed 's/ /__/g'| sed 's/__/ /' | tr ' ' '\\t' | grep "irus" > ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header | sed 's/ /__/g'| sed 's/__/ /' | tr ' ' '\\t' | grep "iroid" >> ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid | awk '{print \$1}' | sort | uniq > ${blastp_nr_bls_ids.baseName}_virus_viroid.ids
        cat ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid | sed 's/__/_/g' | sort -k1 > ${blastp_nr_bls_ids.baseName}_virus_viroid.mod

        grep -F -f ${blastp_nr_bls_ids.baseName}_virus_viroid.ids ${blastp_nr_bls} > ${blastp_nr_bls.baseName}_virus_viroid.bls.txt || [[ \$? == 1 ]]
        cat ${blastp_nr_bls.baseName}_virus_viroid.bls.txt  | sort -k2 >  ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt 
        cut -f2 ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt  | cut -f2 -d '|' > ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction
        paste ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction  | sort -k20 > ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt
        join -1 20 -2 1 ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt ${blastp_nr_bls_ids.baseName}_virus_viroid.mod | sort -u | tr ' ' '\\t' | awk -v OFS='\\t' '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21}' > ${blastp_nr_bls.baseName}_virus_viroid.txt || [[ \$? == 1 ]]
        """
    }   

    process BlastToolsp {
        label "local"
        publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(topHits_blastp_final) from BlastToolsp_ch

        output:
        file "summary_${topHits_blastp_final}"

        script:
        """
        java -jar ${projectDir}/bin/BlastTools.jar -t blastp ${topHits_blastp_final}
        """
    }
}

/*
If the parameter --spades is specified when running nextflow, these processes will derive de novo contigs using SPAdes
*/

if (params.spades) {
    process spades { 
        label "spades_mem_cpu"
        publishDir "${params.outdir}/02a_spades/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(samplefile) from spades_ch

        output:
        file "${sampleid}_${params.minlen}-${params.maxlen}nt_spades/*"
        file "${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta"
        tuple val(sampleid), file("${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta") into cap3_spades_ch
        
        script:
        """
        spades.py -t 2 -k ${params.spadeskmer} --only-assembler -m 180 -s $samplefile -o ${sampleid}_${params.minlen}-${params.maxlen}nt_spades
        
        #edit contigs name and rename velvet assembly
        sed 's/>/>spades_/' ${sampleid}_${params.minlen}-${params.maxlen}nt_spades/contigs.fasta > ${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta
        """
    }

    process cap3_spades {
        label "local"
        publishDir "${params.outdir}/03_cap3/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(scaffolds_fasta) from cap3_spades_ch

        output:
        file "${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.fasta"
        tuple val(sampleid), file("${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.fasta"), file("${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.ids") into blastn_nt_spades_ch

        script:
        """
        cap3 ${scaffolds_fasta} -s 300 -j 31 -i 30 -p 90 -o 16
        cat ${scaffolds_fasta}.cap.singlets ${scaffolds_fasta}.cap.contigs > ${sampleid}_${params.minlen}-${params.maxlen}nt_spades_cap3.fasta
        extract_seqs_rename.py ${sampleid}_${params.minlen}-${params.maxlen}nt_spades_cap3.fasta ${params.cap3_len} | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" | sed 's/|>/ |/' | awk '{print \$1}' > ${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.fasta
        #fetch scaffold IDs
        cat ${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.fasta | grep ">" | sed 's/>//' > ${sampleid}_spades_cap3_${params.minlen}-${params.maxlen}nt.rename.ids
        """
    }

    process blastn_nt_spades {
        label "medium_mem"
        publishDir "${params.outdir}/04_blastn/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(spades_cap3_rename_fasta), file(spades_cap3_rename_fasta_ids) from blastn_nt_spades_ch

        output:
        file "${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT.bls"
        file "${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits.txt"
        file "${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt"
        tuple val(sampleid), file("${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt") into BlastTools_blastn_spades_ch
        script:
        def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
        """
        blastn ${blast_task_param} \
                -query ${spades_cap3_rename_fasta} \
                -db ${blastn_db_name} \
                -out ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT.bls \
                -evalue ${params.blastn_evalue} \
                -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
                -num_threads ${task.cpus} \
                -max_target_seqs 100 \
                -word_size 11

        #fetch top blastn hits
        for i in `cat ${spades_cap3_rename_fasta_ids}`; do
            grep \$i ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT.bls | head -n5 >> ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits.txt;
        done
        grep -i "Virus" ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits.txt > ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        grep -i "Viroid" ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits.txt >> ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        cat ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_spades_${params.minlen}-${params.maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt
        """
    }

    process BlastTools_blastn_spades {
        label "local"
        publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(top5Hits_final) from BlastTools_blastn_spades_ch

        output:
        file "summary_${top5Hits_final}"

        script:
        """
        java -jar ${projectDir}/bin/BlastTools.jar -t blastn ${top5Hits_final}
        """
    }
}
