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
        sampleid,samplepath
        MT019,/user/folder/MT019_sRNA.fastq

      --blast_db_dir '[path/to/files]'                  Path to the blast NT and/or NR database file base name
                                                        '/work/eresearch_bio/reference/blastDB'

      
      --blast_viral_db                                    Run blastn and megablast homology search on cap3 de novo assembly against a virus and viroid database
                                                        [False]

      --blast_viral_nt_db '[path/to/file]'              Path to the viral nucleotide database file base name. Required if --blast_viral_db option is specified
                                                        '/work/hia_mt18005/databases/PVirDB/PVirDB_ver2022_06_03/PVirDB_ver20220603pub.fasta'
      
      --blastn_evalue '[value]'                         Blastn evalue.
                                                        '0.0001'

      --blastn_method ['blastn/megablast']              Specify blastn homology search on cap3 de novo assembly againts NCBI NT
                                                        [default megablast]

      --blastp [True/False]                             Predict ORF from de novo assembled contigs and run blastP againts NCBI NR
                                                        [False]

      --blastp_evalue '[value]'                         Blastp evalue. Required if --blastp option is specified
                                                        '0.0001'
                    
      --blastx [True/False]                             Run blastX againts NCBI NR
                                                        [False]

      --bowtie_db_dir                                   Path to the bowtie indices (for RNA source step and filtering of unusable reads)
      
      
      --cap3_len '[value]'                              Trim value used in the CAP3 step.
                                                        '30'

      --contamination_detection [True/False]            Run false positive prediction due to cross-sample contamination for detections 
                                                        obtained via blastn search against NT
                                                        [False]

      --contamination_detection_viral_db                 Run false positive prediction due to cross-sample contamination for detections 
                                                        obtained via blastn search against a viral database
                                                        [False]

      --contamination_detection_method '[value]'        Either use FPKM or RPM for cross-contamination detection 

      
      --contamination_flag '[value]'                    Threshold value to predict false positives due to cross-sample contamination. 
                                                        Required if --contamination_detection option is specified
                                                        '0.01'

      --dedup                                           Use UMI-tools dedup to remove duplicate reads  
      
      --ictvinfo '[path/to/dir]'                        Path to ICTV info file. Required if --blast_viral_db option is specified
                                                        ['ICTV_taxonomy_MinIdentity_Species.tsv']

      --maxlen '[value]'                                Maximum read length to extract
      ['22']

      --merge_lane                                      Specify this option if sequencing was peformed on several flow cells and 2 or more fastq files were generated for one sample and require to be merged

      --minlen '[value]'                                Minimum read length to extract
      ['21']
      
      --orf_circ_minsize '[value]'                      The value of minsize for getorf -circular
                                                        '75'
      
      --orf_minsize '[value]'                           The value of minsize for getorf
                                                        '75'

      --qualityfilter [True/False]                      Perform adapter and quality filtering of fastq files
                                                        [False]

      --spadesmem  '[value]'                            Memory usage for SPAdes de novo assembler
                                                        [60]               
      
      --targets [True/False]                            Filter the blastn results to viruses/viroids of interest
                                                        [False]

      --targets_file '[path/to/folder]'                 File specifying the name of the viruses/viroids of interest to filter from the blast results output
                                                        ['Targetted_Viruses_Viroids.txt']
      
      --tblastn                                         tblastn homology search on predicted ORFs from getorf against to a viral database
                                                        [False]
      
      --tblastn_evalue                                  tblastn evalue. Required if --tblatsn option is specified
                                                        '0.0001'
      
      --virusdetect [True/False]                        Run VirusDetect
                                                        [False]
      
      --virusdetect_db_path '[path/to/filebasename]'    Path to the virusdetect blast virus database base name
      
      
    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
if (params.blast_db_dir != null) {
    blastn_db_name = "${params.blast_db_dir}/nt"
    blastp_db_name = "${params.blast_db_dir}/nr"
}
if (params.blast_viral_db_path != null) {
    blast_viral_db_name = file(params.blast_viral_db_path).name
    blast_viral_db_dir = file(params.blast_viral_db_path).parent
}
if (params.virusdetect_db_path != null) {
    virusdetect_db_dir = file(params.virusdetect_db_path).parent
}
size_range = "${params.minlen}-${params.maxlen}nt"
if (params.sampleinfo_path != null) {
    sampleinfo_dir = file(params.sampleinfo_path).parent
    sampleinfo_name = file(params.sampleinfo_path).name
}


switch (workflow.containerEngine) {
    case "docker":
        bindbuild = "";
        if (params.blast_viral_db_path != null) {
            bindbuild = "-v ${blast_viral_db_dir}:${blast_viral_db_dir} "
        }
        if (params.blast_db_dir != null) {
            bindbuild = (bindbuild + "-v ${params.blast_db_dir}:${params.blast_db_dir} ")
        }
        if (params.bowtie_db_dir != null) {
            bindbuild = (bindbuild + "-v ${params.bowtie_db_dir}:${params.bowtie_db_dir} ")
        }
        if (params.virusdetect_db_path != null) {
            bindbuild = (bindbuild + "-v ${virusdetect_db_dir}:${virusdetect_db_dir} ")
        }
        if (params.sampleinfo_path != null) {
            bindbuild = (bindbuild + "-v ${sampleinfo_dir}:${sampleinfo_dir} ")
        }
        bindOptions = bindbuild;
        break;
    case "singularity":
        bindbuild = "";
        if (params.blast_viral_db_path != null) {
            bindbuild = "-B ${blast_viral_db_dir} "
        }
        if (params.blast_db_dir != null) {
            bindbuild = (bindbuild + "-B ${params.blast_db_dir} ")
        }
        if (params.bowtie_db_dir != null) {
            bindbuild = (bindbuild + "-B ${params.bowtie_db_dir} ")
        }
        if (params.virusdetect_db_path != null) {
            bindbuild = (bindbuild + "-B ${virusdetect_db_dir} ")
        }
        if (params.sampleinfo_path != null) {
            bindbuild = (bindbuild + "-B ${sampleinfo_dir} ")
        }
        bindOptions = bindbuild;
        break;
    default:
        bindOptions = "";
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
    process FASTQC_RAW {
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

    process MERGE_LANES {
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

    //This step takes > 1h to run for the large flow cells
    process ADAPTER_AND_QUAL_TRIMMING {
        label "setting_6"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link', overwrite: true, pattern: "*{log,json,html,trimmed.fastq.gz,zip,html,pdf,txt}"

        input:
        tuple val(sampleid), file(fastqfile) from umitools_ch

        output:
        file "${sampleid}_umi_tools.log"
        file "${sampleid}_truseq_adapter_cutadapt.log"
        file "${sampleid}_umi_tools.log" into umi_tools_results
        tuple val(sampleid), file("${sampleid}_umi_cleaned.fastq.gz") into rna_profile_ch, qc_post_qual_trimming_ch
        

        script:
        """
        #Checks Illumina seq adapters have been removed
        cutadapt -j ${task.cpus} \
                --no-indels \
                -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;min_overlap=12" \
                -g "ACACTCTTTCCCTACACGACGCTCTTCCGATCT;min_overlap=9" \
                --times 2 \
                -o ${sampleid}_trimmed.fastq.gz \
                ${fastqfile} > ${sampleid}_truseq_adapter_cutadapt.log

        umi_tools extract --extract-method=regex \
                            --bc-pattern=".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})\$" \
                            -I ${sampleid}_trimmed.fastq.gz \
                            -S ${sampleid}_umi_cleaned.fastq.gz > ${sampleid}_umi_tools.log
        
        rm ${sampleid}_trimmed.fastq.gz
        """
    }

    process QC_POST_QUAL_TRIMMING { 
        label "setting_3"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link', overwrite: true, pattern: "*{log,json,html,trimmed.fastq.gz,zip,html,pdf,txt}"
        
        input:
        tuple val(sampleid), file(fastqfile) from qc_post_qual_trimming_ch

        output:
        file "*_fastqc.{zip,html}"
        file "${sampleid}_fastp.json"
        file "${sampleid}_fastp.html"
        file "${sampleid}_read_length_dist.pdf"
        file "${sampleid}_read_length_dist.txt"
        file "${sampleid}_quality_trimmed.fastq.gz"
        file "${sampleid}_qual_filtering_cutadapt.log"

        file "${sampleid}_qual_filtering_cutadapt.log" into cutadapt_qual_filt_results
        tuple val(sampleid), file("${sampleid}_quality_trimmed.fastq") into derive_usable_reads_ch
        file "${sampleid}_fastp.json" into fastp_results
        file "${sampleid}_read_length_dist.txt" into read_length_dist_results

        script:
        """
        cutadapt -j ${task.cpus} \
                --trim-n --max-n 0 -m 18 -q 30 \
                -o ${sampleid}_quality_trimmed.fastq \
                ${sampleid}_umi_cleaned.fastq.gz > ${sampleid}_qual_filtering_cutadapt.log

        pigz --best --force -p ${task.cpus} -r ${sampleid}_quality_trimmed.fastq -c > ${sampleid}_quality_trimmed.fastq.gz

        fastqc --quiet --threads ${task.cpus} ${sampleid}_quality_trimmed.fastq.gz

        fastp --in1=${sampleid}_quality_trimmed.fastq.gz --out1=${sampleid}_fastp_trimmed.fastq.gz \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --json=${sampleid}_fastp.json \
            --html=${sampleid}_fastp.html \
            --thread=${task.cpus}
        
        #derive distribution for quality filtered reads > 5 bp long
        cutadapt -j ${task.cpus} \
                --trim-n --max-n 0 -m 5 -q 30 \
                -o ${sampleid}_quality_trimmed_temp.fastq \
                ${sampleid}_umi_cleaned.fastq.gz
        
        fastq2fasta.pl ${sampleid}_quality_trimmed_temp.fastq > ${sampleid}_quality_trimmed.fasta
        
        read_length_dist.py --input ${sampleid}_quality_trimmed.fasta
        
        mv ${sampleid}_quality_trimmed.fasta_read_length_dist.txt ${sampleid}_read_length_dist.txt
        mv ${sampleid}_quality_trimmed.fasta_read_length_dist.pdf ${sampleid}_read_length_dist.pdf
        rm ${sampleid}_quality_trimmed_temp.fastq
        """
    }

    if (params.rna_source_profile) {
        process RNA_SOURCE_PROFILE {
            label "setting_2"
            tag "$sampleid"
            publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link',overwrite: true, pattern: "*{log}"
            containerOptions "${bindOptions}"

            input:
            tuple val(sampleid), file(fastqfile) from rna_profile_ch
            
            output:
            file "${sampleid}_bowtie.log"
            file "${sampleid}_bowtie.log" into rna_source_bowtie_results

            script:
            """
            cutadapt -j ${task.cpus} \
                    --trim-n --max-n 0 -m 15 -q 30 \
                    -o ${sampleid}_quality_trimmed_temp2.fastq \
                    ${sampleid}_umi_cleaned.fastq.gz

            #derive distribution for quality filtered reads > 15 bp bp long
            echo ${sampleid} > ${sampleid}_bowtie.log;

            count=1
            for rnatype in rRNA plant_pt_mt_other_genes miRNA plant_tRNA plant_noncoding artefacts plant_virus_viroid; do
                if [[ \${count} == 1 ]]; then
                    fastqfile=${sampleid}_quality_trimmed_temp2.fastq
                fi
                echo \${rnatype} alignment: >> ${sampleid}_bowtie.log;
                bowtie -q -v 2 -k 1 -p ${task.cpus} \
                    --un ${sampleid}_\${rnatype}_cleaned_sRNA.fq \
                    -x ${params.bowtie_db_dir}/\${rnatype} \
                    \${fastqfile} \
                    ${sampleid}_\${rnatype}_match 2>>${sampleid}_bowtie.log
                count=\$((count+1))

                if [[ \${count}  > 1 ]]; then
                    fastqfile=${sampleid}_\${rnatype}_cleaned_sRNA.fq
                fi
                rm ${sampleid}_\${rnatype}_match;
            done
            
            rm *cleaned_sRNA.fq
        """
        }

        process RNA_SOURCE_PROFILE_REPORT {
            publishDir "${params.outdir}/00_quality_filtering/qc_report", mode: 'link'

            input:
            file ('*') from rna_source_bowtie_results.collect().ifEmpty([]) 

            output:
            file "read_origin_pc_summary*.txt"
            file "read_origin_counts*.txt"
            file "read_RNA_source*.pdf"
            file "read_origin_detailed_pc*.txt"

            script:
            """
            rna_source_summary.py
            """
        }
    }

    process DERIVE_USABLE_READS {
        label "setting_3"
        tag "$sampleid"
        publishDir "${params.outdir}/00_quality_filtering/${sampleid}", mode: 'link', overwrite: true, pattern: "*{.log,.fastq.gz}"
        containerOptions "${bindOptions}"
        
        input:
        tuple val(sampleid), file(fastqfile) from derive_usable_reads_ch
        
        output:
        file "${sampleid}*_cutadapt.log"
        file "${sampleid}_blacklist_filter.log"
        file "${sampleid}_${params.minlen}-${params.maxlen}nt.fastq.gz"
        
        tuple val(sampleid), file(fastqfile), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into denovo_ch, synthetic_oligos_ch
        tuple val(sampleid), file("${sampleid}_${params.minlen}-${params.maxlen}nt.fastq") into virusdetect_ch
        
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

        cutadapt -j ${task.cpus} -m 18 -M 25 -o ${sampleid}_18-25nt.fastq.gz ${sampleid}_cleaned.fastq > ${sampleid}_18-25nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 21 -M 22 -o ${sampleid}_21-22nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_21-22nt_cutadapt.log
        cutadapt -j ${task.cpus} -m 24 -M 24 -o ${sampleid}_24nt.fastq.gz ${sampleid}_cleaned.fastq > ${sampleid}_24nt_cutadapt.log
        if [[ ${params.minlen} != 21 ]] || [[ ${params.maxlen} != 22 ]]; then
            cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq ${sampleid}_cleaned.fastq > ${sampleid}_${params.minlen}-${params.maxlen}nt_cutadapt.log
        fi

        rm ${sampleid}_24nt.fastq.gz ${sampleid}_18-25nt.fastq.gz

        pigz --best --force -p ${task.cpus} -r ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq -c > ${sampleid}_${params.minlen}-${params.maxlen}nt.fastq.gz
        """
    }

    process QCREPORT {
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

        output:
        file "run_qc_report*.txt"
        file "run_read_size_distribution*.pdf"
        
        script:
        """
        seq_run_qc_report.py
        
        grouped_bar_chart.py
        """
    }
}
// If user does not specify qualityfilter parameter, then only read size selection (using the minlen and maxlen params specified in the nextflow.config file) will be performed on the fastq file specified in the index file
else {
    process READPROCESSING {
        tag "$sampleid"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link'

        input:
        tuple val(sampleid), file(fastqfile) from read_size_selection_ch

        output:
        file "${sampleid}_${size_range}_cutadapt.log"
        file "${sampleid}_${size_range}.fastq"
        tuple val(sampleid), file("unzipped.fastqfile"), file("${sampleid}_${size_range}.fastq") into denovo_ch, synthetic_oligos_ch
        tuple val(sampleid), file("${sampleid}_${size_range}.fastq") into virusdetect_ch

        script:
        """
        if [[ ${fastqfile} == *.gz ]];
        then
            gunzip -c ${fastqfile} > unzipped.fastqfile
        else
            ln ${fastqfile} unzipped.fastqfile
        fi

        cutadapt -j ${task.cpus} -m ${params.minlen} -M ${params.maxlen} -o ${sampleid}_${size_range}.fastq unzipped.fastqfile > ${sampleid}_${size_range}_cutadapt.log
        """
    }
}

// This process performs separate velvet and SPAdes de novo assembly and after merging the assemblies, the contigs are collapsed using cap3
process DENOVO_ASSEMBLY {
    publishDir "${params.outdir}/01_VirReport/${sampleid}/assembly", mode: 'link', overwrite: true, pattern: "*{fasta,log}"
    tag "$sampleid"
    label "setting_1"

    input:
    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size) from denovo_ch

    output:
    file "${sampleid}_velvet_assembly_${size_range}.fasta"
    file "${sampleid}_velvet_log"
    file "${sampleid}_spades_assembly_${size_range}.fasta"
    file "${sampleid}_spades_log"
    file "${sampleid}_cap3_${size_range}.fasta"

    tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${size_range}.fasta") into blastn_nt_cap3_ch, blast_nt_viral_db_cap3_ch
    tuple val(sampleid), file("${sampleid}_cap3_${size_range}.fasta") into getorf_ch
    
    script:
    """
    #run velvet de novo assembler
    echo 'Starting velvet de novo assembly';
    velveth ${sampleid}_velvet_${size_range}_k15 15 -short -fastq ${fastq_filt_by_size}
    velvetg ${sampleid}_velvet_${size_range}_k15 -exp_cov 2

    #edit contigs name and rename velvet assembly
    sed 's/>/>velvet_/' ${sampleid}_velvet_${size_range}_k15/contigs.fa > ${sampleid}_velvet_assembly_${size_range}.fasta
    cp ${sampleid}_velvet_${size_range}_k15/Log ${sampleid}_velvet_log
    
    #run spades de novo assembler
    spades.py --rna -t ${task.cpus} -k 19,21 -m ${params.spadesmem} -s ${fastq_filt_by_size} -o ${sampleid}_spades_k19_21
    #edit contigs name and rename spades assembly

    if [[ ! -s ${sampleid}_spades_k19_21/transcripts.fasta ]]
    then
        touch ${sampleid}_spades_assembly_${size_range}.fasta
    else
        sed 's/>/>spades_/' ${sampleid}_spades_k19_21/transcripts.fasta > ${sampleid}_spades_assembly_${size_range}.fasta
    fi

    cp ${sampleid}_spades_k19_21/spades.log ${sampleid}_spades_log

    #merge velvet and spades assemblies
    cat ${sampleid}_velvet_assembly_${size_range}.fasta ${sampleid}_spades_assembly_${size_range}.fasta > ${sampleid}_merged_spades_velvet_assembly_${size_range}.fasta
    
    #collapse derived contigs
    cap3 ${sampleid}_merged_spades_velvet_assembly_${size_range}.fasta -s 300 -j 31 -i 30 -p 90 -o 16
    cat ${sampleid}_merged_spades_velvet_assembly_${size_range}.fasta.cap.singlets ${sampleid}_merged_spades_velvet_assembly_${size_range}.fasta.cap.contigs > ${sampleid}_cap3_${size_range}_temp.fasta
    
    #retain only contigs > 30 bp long
    extract_seqs_rename.py ${sampleid}_cap3_${size_range}_temp.fasta ${params.cap3_len} \
                             | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" \
                             | sed 's/|>/ |/' | awk '{print \$1}' \
                             > ${sampleid}_cap3_${size_range}.fasta
    """
}

if (params.virreport_viral_db) {
    process BLAST_NT_VIRAL_DB_CAP3 {
        label "setting_4"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/viral_db", mode: 'link', overwrite: true, pattern: "*{vs_viral_db.bls,.txt}"
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${size_range}.fasta") from blast_nt_viral_db_cap3_ch
        
        output:
        file "${sampleid}_cap3_${size_range}_blastn_vs_viral_db.bls"
        file "${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls"
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${size_range}.fasta"), file("${sampleid}_cap3_${size_range}_blastn_vs_viral_db.bls"), file("${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls") into filter_blast_nt_viral_db_cap3_ch

        script:
        """
        #1. blastn search
        blastn -task blastn \
            -query ${sampleid}_cap3_${size_range}.fasta \
            -db ${blast_viral_db_dir}/${blast_viral_db_name} \
            -out ${sampleid}_cap3_${size_range}_blastn_vs_viral_db.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50

        #2. megablast search
        blastn -query ${sampleid}_cap3_${size_range}.fasta \
            -db ${blast_viral_db_dir}/${blast_viral_db_name} \
            -out ${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50
        """
    }

    process FILTER_BLAST_NT_VIRAL_DB_CAP3 {
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/viral_db", mode: 'link', overwrite: true, pattern: "*{.txt}"
        tag "$sampleid"

        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("${sampleid}_cap3_${size_range}.fasta"), file("${sampleid}_cap3_${size_range}_blastn_vs_viral_db.bls"), file("${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls") from filter_blast_nt_viral_db_cap3_ch

        output:
        file "summary_${sampleid}_cap3_${size_range}_*_vs_viral_db.bls_viruses_viroids_ICTV*.txt"
        file "summary_${sampleid}_cap3_${size_range}_*_vs_viral_db.bls_filtered.txt"
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("summary_${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls_viruses_viroids_ICTV.txt") into cov_stats_blast_nt_viral_db_ch
        
        script:
        """
        c1grep() { grep "\$@" || test \$? = 1; }
        #retain 1st blast hit
        for var in ${sampleid}_cap3_${size_range}_megablast_vs_viral_db.bls ${sampleid}_cap3_${size_range}_blastn_vs_viral_db.bls;
            do 
                cat \${var} | awk '{print \$1}' | sort | uniq > \${var}.top1.ids
                for i in `cat \${var}.top1.ids`; do echo "fetching top hits..." \$i; grep \$i \${var} | head -1 >> \${var}.top1Hits.txt ; done
                cat \${var}.top1Hits.txt | sed 's/ /_/g' > \${var}.txt

                #summarise the blast files
                java -jar ${projectDir}/bin/BlastTools.jar -t blastn \${var}.txt

                sequence_length.py --virus_list summary_\${var}.txt --contig_fasta ${sampleid}_cap3_${size_range}.fasta --sample_name ${sampleid} --read_size ${size_range} --out  summary_\${var}_with_contig_lengths.txt

                #only retain hits to plant viruses
                c1grep  "virus\\|viroid\\|Endogenous" summary_\${var}_with_contig_lengths.txt > summary_\${var}_filtered.txt

                sed -i 's/Elephantopus_scaber_closterovirus/Citrus_tristeza_virus/'  summary_\${var}_filtered.txt
                sed -i 's/Hop_stunt_viroid_-_cucumber/Hop_stunt_viroid/' summary_\${var}_filtered.txt
                
                if [[ ! -s summary_\${var}_filtered.txt ]]
                then
                    if [[ ${params.diagno} == true ]]; then
                        #for FILE in summary_\${var}_viruses_viroids_ICTV.txt summary_\${var}_viruses_viroids_ICTV_endemic.txt summary_\${var}_viruses_viroids_ICTV_regulated.txt;
                        for FILE in summary_\${var}_viruses_viroids_ICTV.txt;
                            do
                                echo -e "Species\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tICTV_information" > "\${FILE}"
                            done
                    else
                        echo -e "Species\tsacc\tnaccs\tlength\tslen\tcov\tav-pident\tstitle\tqseqids\tcontig_ind_lengths\tcumulative_contig_len\tcontig_lenth_min\tcontig_lenth_max\tICTV_information" > summary_\${var}_viruses_viroids_ICTV.txt;   
                    fi

                else
                    #fetch unique virus/viroid species name from Blast summary reports
                    cat summary_\${var}_filtered.txt | awk '{print \$7}' | awk -F "|" '{print \$2}'| sort | uniq | sed 's/Species://' > \${var}_uniq.ids

                    #retrieve the best hit for each unique virus/viroid species name by selecting longest alignment (column 3) and highest genome coverage (column 5)
                    touch \${var}_filtered.txt
                    for id in `cat \${var}_uniq.ids`;
                        do
                            grep \${id} summary_\${var}_filtered.txt | sort -k3,3nr -k5,5nr | head -1 >> \${var}_filtered.txt
                        done

                    #print the header of the inital summary_blastn file
                    cat summary_\${var}_with_contig_lengths.txt | head -1 > header

                    #report 1
                    cat header \${var}_filtered.txt > summary_\${var}_viruses_viroids.txt
                    
                    #fetch genus names of identified hits
                    awk '{print \$7}' summary_\${var}_viruses_viroids.txt | awk -F "|" '{print \$2}' | sed 's/Species://' | sed 1d > wanted.names
                
                    #add species to report
                    paste wanted.names \${var}_filtered.txt | sort > summary_\${var}_viruses_viroids.MOD

                    #fetch ICTV information
                    grep -w -F -f wanted.names ${projectDir}/bin/${params.ictvinfo} | sort > wanted.ICTV

                    #join reports with ICTV information
                    join -a1 -1 1 -2 1 -t '\t' summary_\${var}_viruses_viroids.MOD wanted.ICTV |  awk '\$4>=40' > summary_\${var}_viruses_viroids_ICTV

                    #report 2
                    awk '{print "Species" "\\t" \$0 "\\t"  "ICTV_information" }' header > header2
                    cat header2 summary_\${var}_viruses_viroids_ICTV | awk -F"\\t" '\$1!=""&&\$2!=""&&\$3!=""' > summary_\${var}_viruses_viroids_ICTV.txt

                fi
            done
        """
    }

    process COVSTATS_VIRAL_DB {
        tag "$sampleid"
        label "setting_5"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/viral_db", mode: 'link', overwrite: true
        containerOptions "${bindOptions}"
        
        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(samplefile) from cov_stats_blast_nt_viral_db_ch

        output:
        file "${sampleid}_${params.minlen}-${params.maxlen}*"
        file("${sampleid}_${size_range}_top_scoring_targets_with_cov_stats_viral_db.txt") into contamination_flag_viral_db
        
        script:
        """
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${size_range} --blastdbpath ${blast_viral_db_dir}/${blast_viral_db_name} --dedup ${params.dedup} --mode viral_db --cpu ${task.cpus} --diagno ${params.diagno}
        """
    }
    if (params.contamination_detection_viral_db) {
        process CONTAMINATION_DETECTION_VIRAL_DB {
            label "local"
            publishDir "${params.outdir}/01_VirReport/Summary", mode: 'link', overwrite: true
            containerOptions "${bindOptions}"

            input:
            file ('*') from contamination_flag_viral_db.collect().ifEmpty([])

            output:
            file "VirReport_detection_summary*viral_db*.txt"

            script:
            """
            if ${params.sampleinfo}; then
                flag_contamination.py --read_size ${size_range} --threshold ${params.contamination_flag} --method ${params.contamination_detection_method} --viral_db true --diagno ${params.diagno} --dedup ${params.dedup} --sampleinfo ${params.sampleinfo_path}
            else
                flag_contamination.py --read_size ${size_range} --threshold ${params.contamination_flag} --method ${params.contamination_detection_method} --viral_db true --diagno ${params.diagno} --dedup ${params.dedup}
            fi
            """
        }
    }

    process TBLASTN_VIRAL_DB {
        label "setting_4"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/tblastn/viral_db", mode: 'link', overwrite: true
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(cap3_fasta) from getorf_ch
        
        output:
        file "${sampleid}_cap3_${size_range}_getorf.all.fasta"
        file "${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_out.bls"
        file "${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids_final.txt"

        script:
        """
        getorf -sequence ${cap3_fasta} -outseq ${sampleid}_cap3_${size_range}_getorf.fasta  -minsize ${params.orf_minsize}
        getorf -sequence ${cap3_fasta} -circular True -outseq ${sampleid}_cap3_${size_range}_getorf.circular.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_cap3_${size_range}_getorf.fasta ${sampleid}_cap3_${size_range}_getorf.circular.fasta >  ${sampleid}_cap3_${size_range}_getorf.all.fasta
        #cat ${sampleid}_cap3_${size_range}_getorf.all.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_cap3_${size_range}_getorf.all.fasta.ids

        tblastn -query ${sampleid}_cap3_${size_range}_getorf.all.fasta \
            -db ${blast_viral_db_dir}/${blast_viral_db_name} \
            -evalue ${params.tblastn_evalue} \
            -out ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_out.bls \
            -num_threads ${task.cpus} \
            -max_target_seqs 10 \
            -outfmt '6 qseqid sseqid pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid stitle'
        
        grep ">" ${sampleid}_cap3_${size_range}_getorf.all.fasta | sed 's/>//' | cut -f1 -d ' ' | sort | uniq > ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_out.wanted.ids
        for i in `cat ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_out.wanted.ids`; do
            grep \$i ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_out.bls | head -n5 >> ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits.txt;
        done
    
        grep -i "Virus" ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits.txt > ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
        grep -i "Viroid" ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits.txt >> ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        
        #modify accordingly depending on version of viral_db
        cut -f2 ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids.txt | cut -f2 -d '|' > seq_ids.txt
        cut -f20 ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids.txt | cut -f2 -d '|'  | sed 's/Species://' > species_name_extraction.txt
        paste seq_ids.txt ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids.txt  species_name_extraction.txt > ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids_mod.txt
        awk -v OFS='\\t' '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$22}'  ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids_mod.txt | sed 's/ /_/g' > ${sampleid}_cap3_${size_range}_getorf.all_tblastn_vs_viral_db_top5Hits_virus_viroids_final.txt
        """
    }
}
if (params.virreport_ncbi) {
    process BLASTN_NT_CAP3 {
        label "setting_2"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/blastn/NT", mode: 'link', overwrite: true, pattern: "*{vs_NT.bls,_top5Hits.txt,_final.txt,taxonomy.txt}"
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(cap3_fasta) from blastn_nt_cap3_ch

        output:
        file "${cap3_fasta.baseName}_blastn_vs_NT.bls"
        file "${cap3_fasta.baseName}_blastn_vs_NT_top5Hits.txt"
        file "${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt"
        file "summary_${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_final.txt"
        file "summary_${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt"
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file("summary_${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_final.txt"), file("${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt") into blastTools_results_ch
        tuple val(sampleid), file(cap3_fasta), file("${cap3_fasta.baseName}_blastn_vs_NT_top5Hits.txt") into blastx_nt_cap3_ch
        
        script:
        def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
        """
        #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
        if [[ ! -f ${params.blast_db_dir}/taxdb.btd || ! -f ${params.blast_db_dir}/taxdb.bti ]]; then
            update_blastdb.pl taxdb
            tar -xzf taxdb.tar.gz
        else
            cp ${params.blast_db_dir}/taxdb.btd .
            cp ${params.blast_db_dir}/taxdb.bti .
        fi

        blastn ${blast_task_param} \
            -query ${cap3_fasta} \
            -db ${blastn_db_name} \
            -negative_seqidlist ${params.negative_seqid_list} \
            -out ${cap3_fasta.baseName}_blastn_vs_NT.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
            -max_target_seqs 50 \
            -word_size 24

        grep ">" ${cap3_fasta.baseName}.fasta | sed 's/>//' > ${cap3_fasta.baseName}.ids
        
        #fetch top blastn hits
        for i in `cat ${cap3_fasta.baseName}.ids`; do
            grep \$i ${cap3_fasta.baseName}_blastn_vs_NT.bls | head -n5 >> ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits.txt;
        done
        
        grep -i "Virus" ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits.txt > ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_tmp.txt  || [[ \$? == 1 ]]
        grep -i "Viroid" ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits.txt >> ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_tmp.txt || [[ \$? == 1 ]]
        cat ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_tmp.txt | sed 's/ /_/g' > ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt
        cut -f3,26 ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt | sort | uniq > ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
        
        java -jar ${projectDir}/bin/BlastTools.jar -t blastn ${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt

        rm taxdb.btd
        rm taxdb.bti
        
        sequence_length.py --virus_list summary_${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids.txt --contig_fasta ${cap3_fasta.baseName}.fasta --sample_name ${sampleid} --read_size ${size_range} --out summary_${cap3_fasta.baseName}_blastn_vs_NT_top5Hits_virus_viroids_final.txt
        """
    }

    process COVSTATS_NT {
        tag "$sampleid"
        label "setting_6"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/NT", mode: 'link', overwrite: true, pattern: "*{.fa*,.fasta,metrics.txt,scores.txt,targets.txt,stats.txt,log.txt,.bcf*,.vcf.gz*,.bam*}"
        containerOptions "${bindOptions}"
        
        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size), file(samplefile), file(taxonomy) from blastTools_results_ch

        output:
        file "${sampleid}_${size_range}*"
        file("${sampleid}_${size_range}_top_scoring_targets_*with_cov_stats.txt") into contamination_flag
        
        script:
        """
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${size_range} --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name} --dedup ${params.dedup} --cpu ${task.cpus} --mode ncbi --diagno ${params.diagno}
        
        """
    }

    if (params.contamination_detection) {
        process CONTAMINATION_DETECTION {
            label "local"
            publishDir "${params.outdir}/01_VirReport/Summary", mode: 'link', overwrite: true
            containerOptions "${bindOptions}"

            input:
            file ('*') from contamination_flag.collect().ifEmpty([])

            output:
            file "VirReport_detection_summary*.txt"

            script:
            """
            if [[ ${params.sampleinfo} == true ]]; then
                flag_contamination.py --read_size ${size_range} --threshold ${params.contamination_flag} --method ${params.contamination_detection_method} --dedup ${params.dedup} --diagno ${params.diagno} --targets ${params.targets_file} --sampleinfopath ${params.sampleinfo_path}
            else
                flag_contamination.py --read_size ${size_range} --threshold ${params.contamination_flag} --method ${params.contamination_detection_method} --dedup ${params.dedup} --diagno ${params.diagno} --targets ${params.targets_file}
            fi
            """
        }
    }
    //blastx jobs runs out of memory if only given 64Gb
    if (params.blastx) {process BLASTX {
            label "setting_2"
            publishDir "${params.outdir}/01_VirReport/${sampleid}/blastx/NT", mode: 'link', overwrite: true
            tag "$sampleid"
            containerOptions "${bindOptions}"

            input:
            tuple val(sampleid), file(cap3_fasta), file(top5Hits) from blastx_nt_cap3_ch
            
            output:
            file "${cap3_fasta.baseName}_blastx_vs_NT.bls"
            file "${cap3_fasta.baseName}_blastx_vs_NT_top5Hits.txt"
            file "${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids_final.txt"
            file "summary_${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids_final.txt"
            
            script:
            """
            #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
            if [[ ! -f ${params.blast_db_dir}/taxdb.btd || ! -f ${params.blast_db_dir}/taxdb.bti ]]; then
                perl ${projectDir}/bin/update_blastdb.pl taxdb
                tar -xzf taxdb.tar.gz
            else
                cp ${params.blast_db_dir}/taxdb.btd .
                cp ${params.blast_db_dir}/taxdb.bti .
            fi
            #extract contigs with blastn results
            cut -f1 ${top5Hits} | sort | uniq > denovo_contig_name_ids_with_blastn_hits.txt

            #extract all contigs names from de novo assembly
            grep ">" ${cap3_fasta.baseName}.fasta | sed 's/>//' | sort | uniq > denovo_contig_name_ids.txt

            #extract contigs with no blastn results
            grep -v -F -f denovo_contig_name_ids_with_blastn_hits.txt denovo_contig_name_ids.txt | sort  > denovo_contig_name_ids_unassigned.txt || [[ \$? == 1 ]]
            
            perl ${projectDir}/bin/faSomeRecords.pl -f ${cap3_fasta.baseName}.fasta -l denovo_contig_name_ids_unassigned.txt -o ${cap3_fasta.baseName}_no_blastn_hits.fasta

            extract_seqs_rename.py ${cap3_fasta.baseName}_no_blastn_hits.fasta ${params.blastx_len} \
                                    | sed "s/CONTIG/${sampleid}_${params.minlen}-${params.maxlen}_/" \
                                    > ${cap3_fasta.baseName}_no_blastn_hits_${params.blastx_len}nt.fasta

            blastx -query ${cap3_fasta.baseName}_no_blastn_hits_${params.blastx_len}nt.fasta \
                -db ${blastp_db_name} \
                -out ${cap3_fasta.baseName}_blastx_vs_NT.bls \
                -evalue ${params.blastx_evalue} \
                -num_threads ${task.cpus} \
                -outfmt '6 qseqid sseqid pident nident length mismatch gapopen gaps qstart qend qlen qframe sstart send slen evalue bitscore qcovhsp sallseqid sscinames' \
                -max_target_seqs 1

            #grep ">" ${cap3_fasta} | sed 's/>//' > ${cap3_fasta.baseName}.ids
            cut -f1 ${cap3_fasta.baseName}_blastx_vs_NT.bls  | sed 's/ //' | sort | uniq > ${cap3_fasta.baseName}.ids
            
            #fetch top blastn hits
            touch  ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits.txt
            for i in `cat ${cap3_fasta.baseName}.ids`; do
                grep \$i ${cap3_fasta.baseName}}_blastx_vs_NT.bls | head -n5 >> ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits.txt;
            done
            grep -i "Virus" ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits.txt > ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
            grep -i "Viroid" ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits.txt >> ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
            sed 's/ /_/g' ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids.txt  |  awk -v OFS='\\t' '{ print \$2,\$1,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20}' > ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids_final.txt
            
            java -jar ${projectDir}/bin/BlastTools.jar -t blastp ${cap3_fasta.baseName}_blastx_vs_NT_top5Hits_virus_viroids_final.txt
            rm taxdb.btd
            rm taxdb.bti
            """
        }
    }
}

if (params.virusdetect) {
    process VIRUS_DETECT {
        tag "$sampleid"
        label "setting_6"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(samplefile) from virusdetect_ch

        output:
        //file "${sampleid}_${size_range}_temp/*"
        tuple val(sampleid), \
            file(samplefile), \
            file("${sampleid}_${size_range}.combined") into virus_identify_ch

        script:
        """
        virus_detect.pl --thread_num ${task.cpus}  \
                        --reference ${params.virusdetect_db_path} \
                        ${samplefile} \
                        --depth_cutoff 2 

        cp ${sampleid}_${size_range}_temp/${sampleid}_${size_range}.combined .
        """
    }   

    process VIRUS_IDENTIFY {
        publishDir "${params.outdir}/02_VirusDetect", mode: 'link', overwrite: true, pattern: "*/*{references,combined,fa,html,sam,txt,identified,identified_with_depth}"
        tag "$sampleid"
        label "setting_4"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), \
            file(samplefile), \
            file("${sampleid}_${size_range}.combined") from virus_identify_ch

        output:
        file "${sampleid}/*"
        file("${sampleid}_${size_range}.blastn.summary.filtered.txt") into virusdetectblastnsummaryfiltered_flag
        file("${sampleid}_${size_range}.blastn.summary.spp.txt") into virusdetectblastnsummary_flag

        script:
        """
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
                          ${sampleid}_${size_range}.combined

        mv result_${sampleid}_${size_range} ${sampleid}
        mv ${sampleid}_${size_range}.combined ${sampleid}

        #if VirusDetect does not detect a virus hit via blastn, no summary files will be created
        #Exit process
        if [[ ! -f ${sampleid}/${sampleid}_${size_range}.blastn.summary.txt ]]; then
            touch ${sampleid}/${sampleid}_${size_range}.blastn.summary.txt
            echo -e "Sample\tReference\tLength\t%Coverage\t#contig\tDepth\tDepth_Norm\t%Identity\t%Identity_max\t%Identity_min\tGenus\tDescription\tSpecies" | tee ${sampleid}_${size_range}.blastn.summary.filtered.txt >  ${sampleid}/${sampleid}_${size_range}.blastn.summary.filtered.txt
            echo -e "Sample\tReference\tLength\t%Coverage\t#contig\tDepth\tDepth_Norm\t%Identity\t%Identity_max\t%Identity_min\tGenus\tDescription\tSpecies" | tee ${sampleid}_${size_range}.blastn.summary.spp.txt > ${sampleid}/${sampleid}_${size_range}.blastn.summary.spp.txt 
            exit 0
        else
            cp ${sampleid}/${sampleid}_${size_range}.blastn.summary.txt .
        fi

        cut -f2 ${sampleid}_${size_range}.blastn.summary.txt | grep -v Reference > ${sampleid}_${size_range}.blastn_ids.txt
        cp ${params.blast_db_dir}/taxdb.btd .
        cp ${params.blast_db_dir}/taxdb.bti .
        
        touch ${sampleid}_${size_range}.blastn_spp.txt

        for id in `cat ${sampleid}_${size_range}.blastn_ids.txt`;
            do 
                blastdbcmd -db ${blastn_db_name} -entry \${id} -outfmt '%L' | uniq | sed 's/ /_/g' >>  ${sampleid}_${size_range}.blastn_spp.txt
            done
        sed -i '1 i\\Species' ${sampleid}_${size_range}.blastn_spp.txt
        paste ${sampleid}_${size_range}.blastn.summary.txt ${sampleid}_${size_range}.blastn_spp.txt  > ${sampleid}_${size_range}.blastn.summary.spp.txt
        
        #fetch unique virus/viroid species name from Blast summary reports
        cat ${sampleid}_${size_range}.blastn_spp.txt | grep -v Species | sort | uniq  > ${sampleid}_${size_range}.blastn_unique_spp.txt

        head -n1 ${sampleid}_${size_range}.blastn.summary.spp.txt > ${sampleid}_${size_range}.blastn.summary.tmp.txt
        
        for id in `cat ${sampleid}_${size_range}.blastn_unique_spp.txt`;
            do
                grep \${id} ${sampleid}_${size_range}.blastn.summary.spp.txt | sort -k4,4nr | head -1 >> ${sampleid}_${size_range}.blastn.summary.tmp.txt
            done

        grep -v retrovirus ${sampleid}_${size_range}.blastn.summary.tmp.txt > ${sampleid}_${size_range}.blastn.summary.filtered.txt
        for i in ${sampleid}_${size_range}.blastn.summary.spp.txt ${sampleid}_${size_range}.blastn.summary.filtered.txt;
            do
                sed -i 's/Coverage (%)/%Coverage/' \${i}
                sed -i 's/Depth (Norm)/Depth_Norm/' \${i}
                sed -i 's/Iden Max/Identity_max/' \${i}
                sed -i 's/Iden Min/Identity_min/' \${i}
            done


        rm taxdb.btd
        rm taxdb.bti
        cp ${sampleid}_${size_range}.blastn.summary.spp.txt ${sampleid}/${sampleid}_${size_range}.blastn.summary.spp.txt
        cp ${sampleid}_${size_range}.blastn.summary.filtered.txt ${sampleid}/${sampleid}_${size_range}.blastn.summary.filtered.txt
        """
    }

    process VIRUS_DETECT_BLASTN_SUMMARY {
        publishDir "${params.outdir}/02_VirusDetect/Summary", mode: 'link', overwrite: true
        label "local"

        input:
        file ('*') from virusdetectblastnsummary_flag.collect().ifEmpty([])
        file ('*') from virusdetectblastnsummaryfiltered_flag.collect().ifEmpty([])

        output:
        file ("run_summary_top_scoring_targets_virusdetect_${size_range}*.txt")
        file ("run_summary_top_scoring_targets_virusdetect_filtered_${size_range}*.txt")

        script:
        """
        summary_virus_detect.py --read_size ${size_range}
        """
    }
/*
    process VIRUS_DETECT_BLASTN_SUMMARY_FILTERED {
        publishDir "${params.outdir}/02_VirusDetect/Summary", mode: 'link', overwrite: true
        label "local"

        input:
        file ('*') from virusdetectblastnsummaryfiltered_flag.collect().ifEmpty([])

        output:
        file ("run_summary_top_scoring_targets_virusdetect_21-22nt_filtered.txt")

        script:
        """
        touch run_summary_top_scoring_targets_virusdetect_21-22nt_filtered.txt
        echo "Sample\tReference\tLength\tCoverage (%)\t#contig\tDepth\tDepth (Norm)\t%Identity\t%Iden Max\t%Iden Min\tGenus\tDescription\tSpecies" >> run_summary_top_scoring_targets_virusdetect_21-22nt_filtered.txt
        [ -n "\$(find -name '*nt.blastn.summary.filtered.txt' | head -1)" ] && cat *nt.blastn.summary.filtered.txt | sort -k1,1 -k13,13 | grep -v "%Identity" >> run_summary_top_scoring_targets_virusdetect_21-22nt_filtered.txt
        """
    }
*/
}

if (params.synthetic_oligos) {
    process SYNTHETIC_OLIGOS {
        tag "$sampleid"
        label "setting_6"
        publishDir "${params.outdir}/01_VirReport/${sampleid}/synthetic_oligos", mode: 'link', overwrite: true
        
        input:
        tuple val(sampleid), file(fastqfile), file(fastq_filt_by_size) from synthetic_oligos_ch

        output:
        file ("${sampleid}_${size_range}*")

        script:
        """
        synthetic_oligos.py --sample ${sampleid} --rawfastq ${fastqfile} --fastqfiltbysize ${fastq_filt_by_size} --read_size ${size_range}
        """
    }
}
