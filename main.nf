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
negative_seqid_list = "${params.blast_db_dir}/negative_list_out.txt"
blast_local_db_name = file(params.blast_local_db_path).name
blast_local_db_dir = file(params.blast_local_db_path).parent
//bowtie_db_dir = file(params.bowtie_db_dir).name
//blacklist_db_name = file(params.blacklist_db_path).name
//blacklist_db_dir = file(params.blacklist_db_path).parent
virusdetect_db_dir = file(params.virusdetect_db_path).parent

switch (workflow.containerEngine) {
    case "docker":
        bindOptions = "-v ${params.blast_db_dir}:${params.blast_db_dir} -v ${blast_local_db_dir}:${blast_local_db_dir} -v ${bowtie_db_dir}:${bowtie_db_dir} -v ${virusdetect_db_dir}:${virusdetect_db_dir}"
        break;
    case "singularity":
        bindOptions = "-B ${blast_local_db_dir} -B ${params.blast_db_dir} -B ${params.bowtie_db_dir}"
        break;
    default:
        bindOptions = ""
}

if (params.indexfile) {
    if (params.qualityfilter) {
        Channel
          .fromPath(params.indexfile, checkIfExists: true)
          .splitCsv(header:true)
          .map{ row-> tuple(row.sampleid, row.samplepath.substring(row.samplepath.lastIndexOf('/') +1), file(row.samplepath)) }
          .set{ fastqcraw_ch}
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
          .set{ read_size_selection_ch}
        }
    }
    else { exit 1, "Input samplesheet file not specified!" }

if (params.qualityfilter) {
    process fastqc_raw {
        label "setting_1"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile), file(fastqfile_path) from fastqcraw_ch
        
        output:
        file "*_fastqc.{zip,html}"

        script:
        """
        fastqc --quiet --threads ${task.cpus} ${fastqfile_path}
        """
    }

    process merge_lanes {
        tag "$sampleid"

        input:
        tuple val(sampleid), file(samplepath) from samples_ch
        
        output:
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

    process adapter_trimming {
        label "setting_4"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link', overwrite: true, pattern: "*{log,json,html,trimmed.fastq.gz,zip,html}"

        input:
        tuple val(sampleid), file(fastqfile) from umitools_ch

        output:
        file "${sampleid}_umi_tools.log"
        file "${sampleid}_truseq_adapter_cutadapt.log"
        file "${sampleid}_qual_filtering_cutadapt.log"
        file "${sampleid}_quality_trimmed.fastq.gz"
        file "*_fastqc.{zip,html}"
        file "${sampleid}_fastp.json"
        file "${sampleid}_fastp.html"
        file "${sampleid}_fastp.json" into fastp_results
        file "${sampleid}_umi_tools.log" into umi_tools_results
        file "${sampleid}_qual_filtering_cutadapt.log" into cutadapt_qual_filt_results

        tuple val(sampleid), file("${sampleid}_quality_trimmed.fastq") into derive_usable_reads_ch
        tuple val(sampleid), file("${sampleid}_umi_cleaned.fastq") into fastqc_filtered_ch

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
                --trim-n --max-n 0 -m 18 -q 30 \
                -o ${sampleid}_quality_trimmed.fastq \
                ${sampleid}_umi_cleaned.fastq > ${sampleid}_qual_filtering_cutadapt.log

        fastqc --quiet --threads ${task.cpus} ${sampleid}_quality_trimmed.fastq

        fastp --in1=${sampleid}_quality_trimmed.fastq --out1=${sampleid}_fastp_trimmed.fastq \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --json=${sampleid}_fastp.json \
            --html=${sampleid}_fastp.html \
            --thread=${task.cpus}

        pigz --best --force -p ${task.cpus} -r ${sampleid}_quality_trimmed.fastq -c > ${sampleid}_quality_trimmed.fastq.gz
        """
    }

    process qc_fastq_filtered {
        label "setting_2"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link',overwrite: true, pattern: "*{log,pdf,txt}"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fastqfile) from fastqc_filtered_ch
        
        output:
        file "${sampleid}_read_length_dist.pdf"
        file "${sampleid}_read_length_dist.txt"
        file "${sampleid}_bowtie.log"
        file "${sampleid}_bowtie.log" into rna_source_bowtie_results
        file "${sampleid}_read_length_dist.txt" into read_length_dist_results

        script:
        """
        cutadapt -j ${task.cpus} \
                --trim-n --max-n 0 -m 5 -q 30 \
                -o ${sampleid}_quality_trimmed_temp.fastq \
                ${sampleid}_umi_cleaned.fastq
        
        #derive distribution for quality filtered reads > 5 bp long
        perl ${projectDir}/bin/fastq2fasta.pl ${sampleid}_quality_trimmed_temp.fastq > ${sampleid}_quality_trimmed.fasta
        python ${projectDir}/bin/read_length_dist.py --input ${sampleid}_quality_trimmed.fasta
        mv ${sampleid}_quality_trimmed.fasta_read_length_dist.txt ${sampleid}_read_length_dist.txt
        mv ${sampleid}_quality_trimmed.fasta_read_length_dist.pdf ${sampleid}_read_length_dist.pdf
    
        cutadapt -j ${task.cpus} \
                --trim-n --max-n 0 -m 15 -q 30 \
                -o ${sampleid}_quality_trimmed_temp2.fastq \
                ${sampleid}_umi_cleaned.fastq

        #derive distribution for quality filtered reads > 15 bp bp long
        echo ${sampleid} > ${sampleid}_bowtie.log;

        echo 5S rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_5S_rRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/5S_rRNA \
            ${sampleid}_quality_trimmed_temp2.fastq \
            ${sampleid}_5S_rRNA_match 2>>${sampleid}_bowtie.log;


        echo nc SSU and LSU rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_nc_rRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/rRNA \
            ${sampleid}_5S_rRNA_cleaned_sRNA.fq \
            ${sampleid}_rRNA_match 2>>${sampleid}_bowtie.log;

        echo mt rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_mt_rRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_mt_rRNA_genes \
            ${sampleid}_nc_rRNA_cleaned_sRNA.fq \
            ${sampleid}_mt_rRNA_match 2>>${sampleid}_bowtie.log;

        echo pt rRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_pt_rRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_pt_rRNA_genes \
            ${sampleid}_mt_rRNA_cleaned_sRNA.fq \
            ${sampleid}_pt_rRNA_match 2>>${sampleid}_bowtie.log;

        echo mt other genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_mt_other_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_mt_other_genes \
            ${sampleid}_pt_rRNA_cleaned_sRNA.fq \
            ${sampleid}_mt_other_match 2>>${sampleid}_bowtie.log;

        echo pt other genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_pt_other_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_pt_other_genes \
            ${sampleid}_mt_other_cleaned_sRNA.fq \
            ${sampleid}_pt_other_match 2>>${sampleid}_bowtie.log;

        echo plant miRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_miRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_miRNA \
            ${sampleid}_pt_other_cleaned_sRNA.fq \
            ${sampleid}_plant_miRNA_match 2>>${sampleid}_bowtie.log;

        echo other miRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_other_miRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/miRBase \
            ${sampleid}_plant_miRNA_cleaned_sRNA.fq \
            ${sampleid}_other_miRNA_match 2>>${sampleid}_bowtie.log;

        echo tRNA genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_tRNA_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_tRNA \
            ${sampleid}_other_miRNA_cleaned_sRNA.fq \
            ${sampleid}_tRNA_match 2>>${sampleid}_bowtie.log;

        echo plant nc genes alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_plant_nc_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_nc_genes \
            ${sampleid}_plant_tRNA_cleaned_sRNA.fq \
            ${sampleid}_plant_nc_match 2>>${sampleid}_bowtie.log;

        echo plant transposons alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_transposon_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/plant_transposons \
            ${sampleid}_plant_nc_cleaned_sRNA.fq \
            ${sampleid}_transposon_match 2>>${sampleid}_bowtie.log;

        echo PhiX alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_PhiX_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/phiX174 \
            ${sampleid}_transposon_cleaned_sRNA.fq \
            ${sampleid}_PhiX_sRNA_match 2>>${sampleid}_bowtie.log;

        echo Vector alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_UniVec_cleaned_sRNA.fq \
            -x ${params.bowtie_db_dir}/UniVec \
            ${sampleid}_PhiX_cleaned_sRNA.fq \
            ${sampleid}_UniVec_match 2>>${sampleid}_bowtie.log;

        echo Virus and viroid alignment: >> ${sampleid}_bowtie.log;
        bowtie -q -v 1 -k 1 -p ${task.cpus} \
            --un ${sampleid}_final_unaligned_sRNA.fq \
            -x ${params.bowtie_db_dir}/virus \
            ${sampleid}_UniVec_cleaned_sRNA.fq \
            ${sampleid}_viral_match 2>>${sampleid}_bowtie.log;
        """
    }

    process derive_usable_reads {
        label "setting_2"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link', overwrite: true, pattern: "*{.log,.fastq.gz}"
        containerOptions "${bindOptions}"
        
        input:
        tuple val(sampleid), file(fastqfile) from derive_usable_reads_ch
        
        output:
        file "${sampleid}*_cutadapt.log"
        file "${sampleid}_blacklist_filter.log"
        file "${sampleid}_${params.minlen}-${params.maxlen}nt.fastq.gz"
        tuple val(sampleid), file(fastqfile), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into velvet_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into virusdetect_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into spades_ch
        file ("*_18-25nt_cutadapt.log") into cutadapt_18_25nt_results
        file ("*_21-22nt_cutadapt.log") into cutadapt_21_22nt_results
        file ("*_24nt_cutadapt.log") into cutadapt_24nt_results
        file ("*_blacklist_filter.log") into bowtie_usable_read_cat_results

        script:
        """
        bowtie -q -v 1 \
            -k 1 --un ${sampleid}_cleaned.fastq -p ${task.cpus} \
            -x ${params.bowtie_db_dir}/blacklist \
            ${fastqfile} \
            ${sampleid}_blacklist_match 2>${sampleid}_blacklist_filter.log


        cutadapt -j ${task.cpus} -m 18 -M 25 -o ${sampleid}_18-25nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_18-25nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 21 -M 22 -o ${sampleid}_21-22nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_21-22nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 24 -M 24 -o ${sampleid}_24nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_24nt_cutadapt.log
        if [[ ${params.minlen} != 21 ]] || [[ ${params.maxlen} != 22 ]]; then
            cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log
        fi
        pigz --best --force -p ${task.cpus} -r ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq -c > ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq.gz
        """
    }

    process qcreport {
        publishDir "${params.outdir}/00_quality_filtering/qc_report", mode: 'link'

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
        file "run_qc_report.txt"
        file "run_read_size_distribution.pdf"
        file "read_origin_pc_summary.txt"
        file "read_origin_counts.txt"
        file "read_RNA_source.pdf"
        file "read_origin_detailed_pc.txt"

        script:
        """
        seq_run_qc_report.py
        grouped_bar_chart.py
        rna_source_summary.py
        """
    }
}
// If user does not specify qualityfilter parameter, then only read size selection (using the minlen and maxlen params specified in the nextflow.config file) will be performed on the fastq file specified in the index file
else {
    process readprocessing {
        tag "$sampleid"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from read_size_selection_ch

        output:
        file "${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log"
        file "${sampleid}_${params.minlen}-${params.maxlen}nt.fastq"
        tuple val(sampleid), file("unzipped.fastqfile"), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into velvet_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into virusdetect_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into spades_ch

        script:
        """
        if [[ ${fastqfile} == *.gz ]];
        then
            gunzip -c ${fastqfile} > unzipped.fastqfile
        else
            ln ${fastqfile} unzipped.fastqfile
        fi

        cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq unzipped.fastqfile > ${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log
        """
    }
}

process velvet {
    publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link', overwrite: true, pattern: "*{fasta,log}"
    tag "$sampleid"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size) from velvet_ch

    output:
    file "${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta"
    file "${sampleid}_velvet_log"
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta") into velvet_cap3_ch

    script:
    """
    #run velvet de novo assembler
    echo 'Starting velvet de novo assembly';
    velveth ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15 15 -short -fastq ${fastq_filt_by_size}
    velvetg ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15 -exp_cov 2

    #edit contigs name and rename velvet assembly
    sed 's/>/>velvet_/' ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15/contigs.fa > ${sampleid}_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta
    cp ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_k15/Log ${sampleid}_velvet_log
    """
}

process spades {
    publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link', overwrite: true, pattern: "*{fasta,log}"
    tag "$sampleid"
    label "setting_7"

    input:
    tuple val(sampleid), file(fastq_filt_by_size) from spades_ch

    output:
    file "${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta"
    file "${sampleid}_spades_log"
    tuple val(sampleid), file("${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta") into spades_cap3_ch

    script:
    """
    #run spades de novo assembler
    spades.py --rna -t ${task.cpus} -k 19,21 -m ${params.spadesmem} -s ${fastq_filt_by_size} -o ${sampleid}_spades_k19_21
    #edit contigs name and rename spades assembly
    if [[ ! -s ${sampleid}_spades_k19_21/transcripts.fasta ]]
    then
        touch ${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta
    else
        sed 's/>/>spades_/' ${sampleid}_spades_k19_21/transcripts.fasta > ${sampleid}_spades_assembly_${params.minlen}-${params.maxlen}nt.fasta
    fi
    cp ${sampleid}_spades_k19_21/spades.log ${sampleid}_spades_log
    """
}

process cap3 {
    publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link', overwrite: true, pattern: "*{fasta,log}"
    tag "$sampleid"
    

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(velvet_scaffolds_fasta) from velvet_cap3_ch
    tuple val(sampleid), file(spades_assembly_fasta) from spades_cap3_ch

    output:
    file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta"
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta") into blastn_nt_cap3_ch, blast_nt_localdb_cap3_ch
    tuple val(sampleid), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta") into getorf_ch

    script:
    """
    #merge velvet and spades assemblies
    cat ${velvet_scaffolds_fasta} ${spades_assembly_fasta} > ${sampleid}_merged_spades_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta
    #collapse derivedcontigs
    cap3 ${sampleid}_merged_spades_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta -s 300 -j 31 -i 30 -p 90 -o 16
    cat ${sampleid}_merged_spades_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta.cap.singlets ${sampleid}_merged_spades_velvet_assembly_${params.minlen}-${params.maxlen}nt.fasta.cap.contigs > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_temp.fasta
    extract_seqs_rename.py ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_temp.fasta ${params.cap3_len} \
                             | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" \
                             | sed 's/|>/ |/' | awk '{print \$1}'\
                             > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta
    """
}

process blastn_nt_cap3 {
    label "setting_2"
    publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/NT", mode: 'link', overwrite: true, pattern: "*{vs_NT.bls,_top5Hits.txt,_final.txt,taxonomy.txt}"
    tag "$sampleid"
    containerOptions "${bindOptions}"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta") from blastn_nt_cap3_ch

    output:
    file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls"
    file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt"
    file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"
    file "summary_${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"
    //tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt") into blastTools_blastn_velvet_ch
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("summary_${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt") into blastTools_results_ch
    tuple val(sampleid), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta"), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt") into blastx_nt_cap3_ch
    
    script:
    def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
    """
    #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
    if [[ ! -f ${params.blast_db_dir}/taxdb.btd || ! -f ${params.blast_db_dir}/taxdb.bti ]]; then
        perl ${projectDir}/bin/update_blastdb.pl taxdb
        tar -xzf taxdb.tar.gz
    else
        cp ${params.blast_db_dir}/taxdb.btd .
        cp ${params.blast_db_dir}/taxdb.bti .
    fi

    blastn ${blast_task_param} \
        -query ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta \
        -db ${blastn_db_name} \
        -negative_seqidlist ${negative_seqid_list} \
        -out ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls \
        -evalue 0.0001 \
        -num_threads ${task.cpus} \
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
        -max_target_seqs 50

    grep ">" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta | sed 's/>//' > ${sampleid}_cap3_assembly_${params.minlen}-${params.maxlen}nt.ids
    
    #fetch top blastn hits
    for i in `cat ${sampleid}_cap3_assembly_${params.minlen}-${params.maxlen}nt.ids`; do
        grep \$i ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT.bls | head -n5 >> ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt;
    done
    
    grep -i "Virus" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
    grep -i "Viroid" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt >> ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
    cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt
    cut -f3,26 ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt | sort | uniq > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
    
    java -jar ${projectDir}/bin/BlastTools.jar -t blastn ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt

    rm taxdb.btd
    rm taxdb.bti
    """
}

if (params.blastlocaldb) {
    process blast_nt_localdb_cap3 {
        label "setting_2"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/localdb", mode: 'link', overwrite: true, pattern: "*{vs_NT.bls,.txt}"
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta") from blast_nt_localdb_cap3_ch
        
        output:
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls"
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls") into filter_blast_nt_localdb_cap3_ch

        script:
        """
        #1. blastn search
        blastn -task blastn \
            -query ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50

        #2. megablast search
        blastn -query ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50
        """
    }

    process filter_blast_nt_localdb_cap3 {
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/localdb", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls") from filter_blast_nt_localdb_cap3_ch

        output:
        file "summary_${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_*_vs_localdb.bls_viruses_viroids_ICTV*.txt"
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("summary_${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls_viruses_viroids_ICTV.txt") into cov_stats_blast_nt_localdb_ch
        
        script:
        """
        c1grep() { grep "\$@" || test \$? = 1; }
        #retain 1st blast hit
        for var in ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_megablast_vs_localdb.bls ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_localdb.bls;
            do 
                cat \${var} | awk '{print \$1}' | sort | uniq > \${var}.top1.ids
                for i in `cat \${var}.top1.ids`; do echo "fetching top hits..." \$i; grep \$i \${var} | head -1 >> \${var}.top1Hits.txt ; done
                cat \${var}.top1Hits.txt | sed 's/ /_/g' > \${var}.txt

                #summarise the blast files
                java -jar ${projectDir}/bin/BlastTools.jar -t blastn \${var}.txt

                #only retain hits to plant viruses (regulated/edemic/LandPlant)
                c1grep  "virus\\|viroid" summary_\${var}.txt > summary_\${var}_filtered.txt

                if [[ ! -s summary_\${var}_filtered.txt ]]
                then
                    for FILE in summary_\${var}_viruses_viroids_ICTV.txt summary_\${var}_viruses_viroids_ICTV_endemic.txt summary_\${var}_viruses_viroids_ICTV_regulated.txt;
                        do
                            echo -e "Species\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tICTV_information" > "\${FILE}"
                        done
                else
                    #fetch unique virus/viroid species name from Blast summary reports
                    cat summary_\${var}_filtered.txt | awk '{print \$7}' | awk -F "|" '{print \$2}'| sort | uniq | sed 's/Species://' > \${var}_uniq.ids

                    #retrieve the best hit for each unique virus/viroid species name by selecting longest alignment (column 3) and highest genome coverage (column 5)
                    touch \${var}_filtered.txt
                    for id in `cat \${var}_uniq.ids`;
                        do
                            grep \${id} summary_\${var}.txt | sort -k3,3nr -k5,5nr | head -1 >> \${var}_filtered.txt
                        done

                    #print the header of the inital summary_blastn file
                    cat summary_\${var}.txt | head -1 > header

                    #report 1
                    cat header \${var}_filtered.txt > summary_\${var}_viruses_viroids.txt
                    
                    #fetch genus names of identified hits
                    awk '{print \$7}' summary_\${var}_viruses_viroids.txt | awk -F "|" '{print \$2}' | sed 's/Species://' | sed 1d > wanted.names
                
                    #add species to report
                    paste wanted.names \${var}_filtered.txt | sort > summary_\${var}_viruses_viroids.MOD

                    #fetch ICTV information
                    grep -w -F -f wanted.names ${projectDir}/bin/${params.ictvinfo} | sort > wanted.ICTV

                    #join reports with ICTV information
                    #join -a 1 -1 1 -2 1 summary_\${var}_viruses_viroids.MOD wanted.ICTV | tr ' ' '\\t' | awk '\$4>=70' >  summary_\${var}_viruses_viroids_ICTV
                    join -a1 -1 1 -2 1 summary_\${var}_viruses_viroids.MOD wanted.ICTV | tr ' ' '\\t' >  summary_\${var}_viruses_viroids_ICTV

                    #report 2
                    awk '{print "Species" "\\t" \$0 "\\t" "ICTV_information"}' header > header2
                    cat header2 summary_\${var}_viruses_viroids_ICTV | awk -F"\\t" '\$1!=""&&\$2!=""&&\$3!=""' > summary_\${var}_viruses_viroids_ICTV.txt
                    c1grep "ICTV_information\\|regulated" summary_\${var}_viruses_viroids_ICTV.txt > summary_\${var}_viruses_viroids_ICTV_regulated.txt
                    c1grep "ICTV_information\\|endemic" summary_\${var}_viruses_viroids_ICTV.txt > summary_\${var}_viruses_viroids_ICTV_endemic.txt
                fi
            done
        """
    }

    process covstats_localdb {
    tag "$sampleid"
    label "setting_2"
    publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/localdb", mode: 'link'
    containerOptions "${bindOptions}"
    
    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(samplefile) from cov_stats_blast_nt_localdb_ch

    output:
    file "${sampleid}_${params.minlen}-${params.maxlen}*"
    file("${sampleid}_${params.minlen}-${params.maxlen}nt_top_scoring_targets_with_cov_stats_PVirDB.txt") into contamination_flag_localdb
    
    script:
    """
    if [[ ${params.dedup} == true ]]; then
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --blastdbpath ${blast_local_db_dir}/${blast_local_db_name} --dedup true --mode local --cpu ${task.cpus}
    else
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --blastdbpath ${blast_local_db_dir}/${blast_local_db_name} --dedup false --mode local --cpu ${task.cpus}
    fi
    """
    }
}

process filter_n_cov {
    tag "$sampleid"
    label "setting_2"
    publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/NT", mode: 'link'
    containerOptions "${bindOptions}"
    
    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(samplefile), file(taxonomy) from blastTools_results_ch

    output:
    file "${sampleid}_${params.minlen}-${params.maxlen}*"
    file("${sampleid}_${params.minlen}-${params.maxlen}nt_top_scoring_targets_with_cov_stats.txt") into contamination_flag
    
    script:
    """

    if [[ ${params.targets} == true ]]; then
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name} --dedup ${params.dedup} --cpu ${task.cpus} --targets --targetspath ${projectDir}/bin/${params.targets_file} --mode NT 
    else
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${params.minlen}-${params.maxlen}nt --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name} --dedup ${params.dedup} --cpu ${task.cpus} --mode NT
    fi
    """
}

if (params.contamination_detection) {
    process contamination_detection {
        label "local"
        publishDir "${params.outdir}/01_VirReport/summary", mode: 'link'
        
        input:
        file ('*') from contamination_flag.collect().ifEmpty([])

        output:
        file "run_top_scoring_targets_with_cov_stats_with_cont_flag*.txt"

        script:
        """
        flag_contamination.py --read_size ${params.minlen}-${params.maxlen}nt --threshold ${params.contamination_flag} --method ${params.contamination_detection_method}
        """
    }
}

if (params.contamination_detection_PVirDB) {
    process contamination_detection_PVirDB {
        label "local"
        publishDir "${params.outdir}/01_VirReport/summary", mode: 'link'
        
        input:
        file ('*') from contamination_flag_localdb.collect().ifEmpty([])

        output:
        file "run_top_scoring_targets_with_cov_stats_with_cont_flag*PVirDB.txt"

        script:
        """
        flag_contamination.py --read_size ${params.minlen}-${params.maxlen}nt --threshold ${params.contamination_flag} --method ${params.contamination_detection_method} --PVirDB true
        """
    }
}

if (params.blastx) {
    process blastx {
        label "setting_3"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastx/NT", mode: 'link', overwrite: true
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        //tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta") from blastx_nt_cap3_ch
        tuple val(sampleid), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta"), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt") from blastx_nt_cap3_ch
        
        output:
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT.bls"
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits.txt"
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt"
        file "summary_${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt"
        
        
        script:
        """
        #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
        cp ${params.blast_db_dir}/taxdb.btd .
        cp ${params.blast_db_dir}/taxdb.bti .

        #extract contigs with blastn results
        cut -f1 ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits.txt | sort | uniq > denovo_contig_name_ids_with_blastn_hits.txt

        #extract all contigs names from de novo assembly
        grep ">" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta | sed 's/>//' | sort | uniq > denovo_contig_name_ids.txt

        #extract contigs with no blastn results
        grep -v -F -f denovo_contig_name_ids_with_blastn_hits.txt denovo_contig_name_ids.txt | sort  > denovo_contig_name_ids_unassigned.txt || [[ \$? == 1 ]]
        
        perl ${projectDir}/bin/faSomeRecords.pl -f ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta -l denovo_contig_name_ids_unassigned.txt -o ${sampleid}_cap3_no_blastn_hits.fasta

        extract_seqs_rename.py ${sampleid}_cap3_no_blastn_hits.fasta 75 \
                                | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" \
                                > ${sampleid}_cap3_no_blastn_hits_75bp.fasta

        blastx -query ${sampleid}_cap3_no_blastn_hits_75bp.fasta \
            -db ${blastp_db_name} \
            -negative_seqidlist ${negative_seqid_list} \
            -out ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT.bls \
            -evalue 0.0001 \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sseqid pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid sscinames' \
            -max_target_seqs 1

        #grep ">" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt.fasta | sed 's/>//' > ${sampleid}_cap3_assembly_${params.minlen}-${params.maxlen}nt.ids
        cut -f1 ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT.bls  | sed 's/ //' | sort | uniq > ${sampleid}_cap3_assembly_${params.minlen}-${params.maxlen}nt.ids
        
        
        #fetch top blastn hits
        touch  ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits.txt
        for i in `cat ${sampleid}_cap3_assembly_${params.minlen}-${params.maxlen}nt.ids`; do
            grep \$i ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT.bls | head -n5 >> ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits.txt;
        done
        grep -i "Virus" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits.txt > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
        grep -i "Viroid" ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits.txt >> ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
        cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt
        #cut -f3,26 ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt | sort | uniq > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
        
        java -jar ${projectDir}/bin/BlastTools.jar -t blastp ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_blastx_vs_NT_top5Hits_virus_viroids_final.txt
        """
    }
}


if (params.blastp || params.tblastn) {
    process getorf {
        label 'local'
        publishDir "${params.outdir}/01_VirReport/${sampleid}/tblastn/localdb", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(cap3_fasta) from getorf_ch
        
        output:
        //file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"
        //file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta"
        //file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta"
        tuple val(sampleid), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min35aa.fasta") into tblastn_ch, blastp_ch
        
        script:
        """
        getorf -sequence ${cap3_fasta} -outseq ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min35aa.fasta  -minsize ${params.orf_minsize}
        getorf -sequence ${cap3_fasta} -circular True -outseq ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min35aa.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min35aa.fasta ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min35aa.fasta >  ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min35aa.fasta
        #cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min35aa.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min35aa.fasta.ids
        """
    }
}


if (params.blastp) {
    
    process blastp {
        label "setting_3"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastp", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fasta) from blastp_ch
        
        output:
        file "${fasta}"
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
        label "setting_5"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastp", mode: 'link'
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
        cat ${blastp_nr_bls.baseName}_virus_viroid.bls.txt | sort -k2 >  ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt 
        cut -f2 ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt  | cut -f2 -d '|' > ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction
        paste ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction  | sort -k20 > ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt
        join -1 20 -2 1 ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt ${blastp_nr_bls_ids.baseName}_virus_viroid.mod | sort -u | tr ' ' '\\t' | awk -v OFS='\\t' '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21}' > ${blastp_nr_bls.baseName}_virus_viroid.txt || [[ \$? == 1 ]]
        """
    }   

    process BlastToolsp {
        label "local"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastp", mode: 'link'
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
/*process getorf {
        label 'local'
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastp", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(velvet_cap3_fasta) from getorf_ch
        
        output:
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"
        file "${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta"
        tuple val(sampleid), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"), file("${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids") into blastp_ch
        
        script:
        """
        getorf -sequence ${velvet_cap3_fasta} -outseq ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta  -minsize ${params.orf_minsize}
        cat ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta.ids
        getorf -sequence ${velvet_cap3_fasta} -circular True -outseq ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta.ids
        """
    }
*/

if (params.tblastn) {

    process tblastn_localdb {
        label "setting_3"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/tblastn/localdb", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fasta) from tblastn_ch
        
        output:
        file "${fasta.baseName}_tblastn_vs_localdb_out.bls"
        //file "${fasta.baseName}_tblastn_vs_localdb_out.wanted.ids"
        file "${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_final.txt"
        //file "${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_mod.txt"
        //file "summary_${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_final.txt"
        //file "seq_ids.txt"
        //file "species_name_extraction.txt"
        //tuple val(sampleid),  file("${fasta.baseName}_tblastn_vs_localdb_out.bls") into BlastTool_tblastn_localdb_ch
        
        script:
        """         
        tblastn -query ${fasta} \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -evalue ${params.tblastn_evalue} \
            -out ${fasta.baseName}_tblastn_vs_localdb_out.bls \
            -num_threads ${task.cpus} \
            -max_target_seqs 10 \
            -outfmt '6 qseqid sseqid pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid stitle'

            
            ###-outfmt '6 qseqid sgi sacc pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid stitle'
        
        grep ">" ${fasta} | sed 's/>//' | cut -f1 -d ' ' | sort | uniq > ${fasta.baseName}_tblastn_vs_localdb_out.wanted.ids
        for i in `cat ${fasta.baseName}_tblastn_vs_localdb_out.wanted.ids`; do
            grep \$i ${fasta.baseName}_tblastn_vs_localdb_out.bls | head -n5 >> ${fasta.baseName}_tblastn_vs_localdb_top5Hits.txt;
        done
    
        grep -i "Virus" ${fasta.baseName}_tblastn_vs_localdb_top5Hits.txt > ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
        grep -i "Viroid" ${fasta.baseName}_tblastn_vs_localdb_top5Hits.txt >> ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        #modify accordingly depending on version of localdb
        cut -f3 ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids.txt | cut -f2 -d '|' > seq_ids.txt
        cut -f20 ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids.txt | cut -f3 -d '|'  | sed 's/Species://' > species_name_extraction.txt
        paste seq_ids.txt ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids.txt  species_name_extraction.txt > ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_mod.txt
        awk -v OFS='\\t' '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$22}'  ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_mod.txt | sed 's/ /_/g' > ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_final.txt
        
        #java -jar ${projectDir}/bin/BlastTools.jar -t blastp ${fasta.baseName}_tblastn_vs_localdb_top5Hits_virus_viroids_final.txt
        """
    }
}


    /*process getorf_localdb {
        label 'local'
        publishDir "${params.outdir}/01_VirReport/${sampleid}/tblastn/localdb", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(cap3_fasta) from getorf_ch
        
        output:
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta"
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta"
        file "${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta"
        tuple val(sampleid), file("${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta") into tblastn_ch
        
        script:
        """
        getorf -sequence ${cap3_fasta} -outseq ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta  -minsize ${params.orf_minsize}
        getorf -sequence ${cap3_fasta} -circular True -outseq ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.min50aa.fasta ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.circular.min50aa.fasta >  ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta
        #cat ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_cap3_${params.minlen}-${params.maxlen}nt_getorf.all.min50aa.fasta.ids
        """
    }
    */
/*
    process tblastndbcmd {
        label "setting_5"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/tblastn/localdb", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(tblastn_nr_bls_ids), file(tblastn_nr_bls) from tblastndbcmd_ch

        output:
        file "${tblastn_nr_bls.baseName}_virus_viroid.txt"
        tuple val(sampleid), file("${tblastn_nr_bls.baseName}_virus_viroid.txt") into BlastTool_tblastn_localdb_ch

        script:
        """
        blastdbcmd -db ${blast_local_db_dir}/${blast_local_db_name} \
                   -dbtype nucl \
                   -entry_batch ${blastp_nr_bls_ids} > ${blastp_nr_bls_ids.baseName}.out

        grep ">" ${blastp_nr_bls_ids.baseName}.out | sed 's/>//' > ${blastp_nr_bls_ids.baseName}.out.header || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header | sed 's/ /__/g'| sed 's/__/ /' | tr ' ' '\\t' | grep "irus" > ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header | sed 's/ /__/g'| sed 's/__/ /' | tr ' ' '\\t' | grep "iroid" >> ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid || [[ \$? == 1 ]]
        cat ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid | awk '{print \$1}' | sort | uniq > ${blastp_nr_bls_ids.baseName}_virus_viroid.ids
        cat ${blastp_nr_bls_ids.baseName}.out.header_virus_viroid | sed 's/__/_/g' | sort -k1 > ${blastp_nr_bls_ids.baseName}_virus_viroid.mod

        grep -F -f ${blastp_nr_bls_ids.baseName}_virus_viroid.ids ${blastp_nr_bls} > ${blastp_nr_bls.baseName}_virus_viroid.bls.txt || [[ \$? == 1 ]]
        cat ${blastp_nr_bls.baseName}_virus_viroid.bls.txt | sort -k2 >  ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt 
        cut -f2 ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt  | cut -f2 -d '|' > ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction
        paste ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.txt ${blastp_nr_bls.baseName}_virus_viroid_sorted.bls.id.extraction  | sort -k20 > ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt
        join -1 20 -2 1 ${blastp_nr_bls.baseName}_virus_viroid_sorted2.bls.txt ${blastp_nr_bls_ids.baseName}_virus_viroid.mod | sort -u | tr ' ' '\\t' | awk -v OFS='\\t' '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21}' > ${blastp_nr_bls.baseName}_virus_viroid.txt || [[ \$? == 1 ]]
        """
    }   
/*
    process BlastToolsp {
        label "local"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastp", mode: 'link'
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
/*
if (params.spades) {
    process spades { 
        label "setting_7"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link'
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
        publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link'
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
        label "setting_5"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blast", mode: 'link'
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
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blast", mode: 'link'
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
*/
if (params.virusdetect) {
    process virus_detect {
        publishDir "${params.outdir}/02_virusdetect/${sampleid}", mode: 'link'
        tag "$sampleid"
        label "setting_6"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(samplefile) from virusdetect_ch

        output:
        file "${sampleid}_${params.minlen}-${params.maxlen}nt_temp/*"
        file "result_${sampleid}_${params.minlen}-${params.maxlen}nt/*"

        script:
        """
        virus_detect.pl --thread_num ${task.cpus}  \
                        --reference ${params.virusdetect_db_path} \
                        ${samplefile} \
                        --depth_cutoff 2 
        
        cp ${sampleid}_${params.minlen}-${params.maxlen}nt_temp/${sampleid}_${params.minlen}-${params.maxlen}nt.combined .
    
        virus_identify.pl --reference ${params.virusdetect_db_path} \
                        --word-size 11 \
                        --exp-value 1e-05 \
                        --exp-valuex 0.01 \
                        --percent-identity 25 \
                        --cpu-num ${task.cpus}  \
                        --mis-penalty -3 \
                        --gap-cost -1 \
                        --gap-extension -1 \
                        --hsp-cover 0.75 \
                        --diff-ratio 0.25 \
                        --diff-contig-cover 0.5 \
                        --diff-contig-length 100 \
                        --coverage-cutoff 0.1 \
                        --depth-cutoff 2 \
                        --siRNA-percent 0.5 \
                        --novel-len-cutoff 100 \
                        --debug \
                        ${samplefile} \
                        ${sampleid}_${params.minlen}-${params.maxlen}nt.combined
        """
    }
}