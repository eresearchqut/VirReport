#!/usr/bin/env nextflow
/*
Virus Surveillance and Diagnosis (VSD) workflow
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

    Virus Surveillance and Diagnosis (VSD) workflow
    Roberto Barrero, 14/03/2019
    Desmond Schmidt, 2/7/2019
    Converted to Nextflow by Craig Windell, 11/2020
    Modified by Maely Gauthier 12/2021

    Usage:

    Run the command
    nextflow run eresearch/vsd -profile ...

    Mandatory arguments:
       -profile '[docker, singularity, conda]'      Profile to use. Choose docker, or singularity, or conda

    Optional arguments:
      --indexfile '[path/to/file]'                  Path to the csv file that contains the list of
                                                    samples to be analysed by this pipeline.
                                                    'index.csv'
      Contents of indexfile csv:
        sampleid,samplepath,minlen,maxlen
        MT019,MT019_sRNA.fastq,21,22

      --blast_nt_db '[path/to/dir]'                 Path to the blast NT database files
                                                    '/work/eresearch_bio/reference/blastDB/nt'

      --blast_nr_db '[path/to/dir]'                 Path to the blast NR database files
                                                    '/work/eresearch_bio/reference/blastDB/nr'

      --cap3_len '[value]'                          Trim value used in the CAP3 step.
                                                    '20'

      --orf_minsize '[value]'                       The value of minsize for getorf
                                                    '150'

      --orf_circ_minsize '[value]'                  The value of minsize for getorf -circular
                                                    '150'

      --blastn_evalue '[value]'                     Blastn evalue.
                                                    '0.0001'
    
      --targets [True/False]                        Filter the blastn results to viruses/viroids of interest
                                                    [False]

      --targets_file '[path/to/dir]'                File specifying the name of the viruses/viroids of interest to filter from the blast results output
                                                    ['/home/gauthiem/code/vsd-2.0/Targetted_Viruses_Viroids.txt']

      --blastn_method ['blastn/megablast']      Run blastn homology search on velvet de novo assembly againts NCBI NT
                                                [default blastn]

      --blastlocaldb                                Run blastn and megablast homology search on velvet de novo assembly against local virus and viroid database
                                                    [False]

      --blast_local_nt_db '[path/to/dir]'           Path to the local blast NT database files. Required if --blastlocaldb option is specified
                                                    '/work/hia_mt18005/databases/sequences/PVirDB_20210330'

      --ictvinfo '[path/to/dir]'                    Path to ICTV info file. Required if --blastlocaldb option is specified
                                                    ['ICTV_taxonomy_MinIdentity_Species.tsv']

      --blastp [True/False]                         Predict ORF from de novo assembled contigs and run blastP againts NCBI NR
                                                    [False]

      --blastp_evalue '[value]'                     Blastp's evalue. Required if --blastp option is specified
                                                    '0.0001'
                    
      --spades [True/False]                         Run SPAdes 3.14 de novo assembler and perform blastn homology analysis on the derived de novo contigs
                                                    [False]

      --spadeskmer '[value]'                        K-mer range to use when running SPAdes. Required if --spades option is specified
                                                    ['K9-21']

      --contamination_detection [True/False]        Run false positive prediction due to cross-sample contamination
                                                    [False]

      --contamination_flag '[value]'                Threshold value to predict false positives due to cross-sample contamination. 
                                                    Required if --contamination_detection option is specified
                                                    '0.01'
      --contamination_detection_method '[value]'    Either use RPKM or Reads per million for cross-contamination detection 

    Other options

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

blastn_db_name = "${params.blast_db_dir}/nt"
blastp_db_name = "${params.blast_db_dir}/nr"
blast_local_db_name = file(params.blast_local_db_path).name
blast_local_db_dir = file(params.blast_local_db_path).parent

switch (workflow.containerEngine) {
    case "docker":
        bindOptions = "-v ${params.blast_db_dir}:${params.blast_db_dir} -v ${blast_local_db_dir}:${blast_local_db_dir}"
        break;
    case "singularity":
        bindOptions = "-B ${blast_local_db_dir} -B ${params.blast_db_dir}"
        break;
    default:
        bindOptions = ""
}

if (params.indexfile) {
  Channel
    .fromPath(params.indexfile, checkIfExists: true)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleid, file(row.samplepath), row.minlen, row.maxlen) }
    .into{ read_size_selection_ch; filter_n_cov_ch; contamination_detection_ch}
    //.groupTuple()
    //.view()
    //.set { samples_ch }
    //[MT019, /work/eresearch_bio/nextflow/plant_biosecurity/workflow/MT019_sRNA.fastq, 21, 22]
} else { exit 1, "Input samplesheet file not specified!" }

process readprocessing {
    tag "$sampleid"
    publishDir "${params.outdir}/01_read_size_selection", mode: 'link'

    input:
    tuple val(sampleid), file(fastqfile), val(minlen), val(maxlen) from read_size_selection_ch

    output:
    file "${sampleid}_${minlen}-${maxlen}nt_cutadapt.log"
    file "${sampleid}_${minlen}-${maxlen}nt.fastq"
    tuple val(sampleid), file("${sampleid}_${minlen}-${maxlen}nt.fastq"), val(minlen), val(maxlen) into velvet_ch, spades_ch
    tuple val(sampleid), file("${sampleid}_${minlen}-${maxlen}nt.fastq") into fastq_filt_by_size_ch


    script:
    """
    cutadapt -j ${task.cpus} -m ${minlen} -M ${maxlen} -o ${sampleid}_${minlen}-${maxlen}nt.fastq ${fastqfile} > ${sampleid}_${minlen}-${maxlen}nt_cutadapt.log
    """
}

process velvet {
    publishDir "${params.outdir}/02_velvet/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(samplefile), val(minlen), val(maxlen) from velvet_ch

    output:
    file "${sampleid}_velvet_${minlen}-${maxlen}nt_k15/*"
    tuple val(sampleid), file("${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.fasta"), \
        val(minlen), val(maxlen) into cap3_ch

    script:
    """
    #run velvet de novo assembler
    echo 'Starting velvet de novo assembly';
    velveth ${sampleid}_velvet_${minlen}-${maxlen}nt_k15 15 -short -fastq $samplefile
    velvetg ${sampleid}_velvet_${minlen}-${maxlen}nt_k15 -exp_cov 2

    #edit contigs name and rename velvet assembly
    sed 's/>/>velvet_/' ${sampleid}_velvet_${minlen}-${maxlen}nt_k15/contigs.fa > ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.fasta
    """
}

process cap3 {
    label "local"
    publishDir "${params.outdir}/03_cap3/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(scaffolds_fasta), val(minlen), val(maxlen) from cap3_ch

    output:
    tuple val(sampleid), file("${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta"), val(minlen), val(maxlen) into blast_nt_localdb_velvet_ch, blastn_nt_velvet_ch, getorf_ch
   
    script:
    """
    #collapse velvet contigs
    cap3 ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.fasta -s 300 -j 31 -i 30 -p 90 -o 16
    cat ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.fasta.cap.singlets ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.fasta.cap.contigs > ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt.fasta
    extract_seqs_rename.py ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt.fasta ${params.cap3_len} \
                             | sed "s/CONTIG/${sampleid}_${minlen}-${maxlen}_/" \
                             | sed 's/|>/ |/' | awk '{print \$1}'\
                             > ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta
    """
}

process blastn_nt_velvet {
    label "blastn_mem"
    publishDir "${params.outdir}/04_blastn/${sampleid}", mode: 'link'
    tag "$sampleid"
    containerOptions "${bindOptions}"

    input:
    tuple val(sampleid), file("${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta"), val(minlen), val(maxlen) from blastn_nt_velvet_ch

    output:
    file "${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT.bls"
    file "${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits.txt"
    file "${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"
    tuple val(sampleid), file("${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt"), file("${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt") into blastTools_blastn_velvet_ch
    
    script:
    def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
    """
    #To extract the taxonomy, copy the taxonomy databases associated with your blast NT database
    cp ${params.blast_db_dir}/taxdb.btd .
    cp ${params.blast_db_dir}/taxdb.bti .
    blastn ${blast_task_param} \
        -query ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta \
        -db ${blastn_db_name} \
        -out ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT.bls \
        -evalue 0.0001 \
        -num_threads 4 \
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
        -max_target_seqs 50

    grep ">" ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta | sed 's/>//' > ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.ids
    
    #fetch top blastn hits
    for i in `cat ${sampleid}_velvet_assembly_${minlen}-${maxlen}nt.ids`; do
        grep \$i ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT.bls | head -n5 >> ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits.txt;
    done
    
    grep -i "Virus" ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits.txt > ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt  || [[ \$? == 1 ]]
    grep -i "Viroid" ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits.txt >> ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
    cat ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt
    cut -f3,26 ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_final.txt | sort | uniq > ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_NT_top5Hits_virus_viroids_seq_ids_taxonomy.txt
    """
}

process BlastTools_blastn_velvet {
    label "medium_mem"
    publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
    tag "$sampleid"

    input:
    tuple val(sampleid), file(top5Hits_final), file(taxonomy) from blastTools_blastn_velvet_ch

    output:
    file "summary_${top5Hits_final}"
    tuple val(sampleid), "summary_${top5Hits_final}", "${taxonomy}" into blastTools_results_ch

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
        tuple val(sampleid), file("${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta"), val(minlen), val(maxlen) from blast_nt_localdb_velvet_ch
        
        output:
        file "${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls"
        file "${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls"
        tuple val(sampleid), file("${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls"), val(minlen), val(maxlen) into filter_blast_nt_localdb_velvet_ch

        script:
        """
        #1. blastn search
        blastn -task blastn \
            -query ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50

        #2. megablast search
        blastn -query ${sampleid}_velvet_cap3_${minlen}-${maxlen}nt_rename.fasta \
            -db ${blast_local_db_dir}/${blast_local_db_name} \
            -out ${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls \
            -evalue ${params.blastn_evalue} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
            -max_target_seqs 50
        """
    }

    process filter_blast_nt_localdb_velvet {
        label "blastn_mem"
        publishDir "${params.outdir}/05_blastoutputs/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file("${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls"), file("${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls"), val(minlen), val(maxlen) from filter_blast_nt_localdb_velvet_ch

        output:
        file "summary_${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls_ENDEMIC_viruses_viroids_ICTV.txt"
        file "summary_${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls_REGULATED_viruses_viroids_ICTV.txt"
        file "summary_${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls_ENDEMIC_viruses_viroids_ICTV.txt"
        file "summary_${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls_REGULATED_viruses_viroids_ICTV.txt"

        script:
        """
        bash run_VSD_report.sh ${sampleid}_velvet_${minlen}-${maxlen}nt_megablast_vs_localdb.bls ${params.ictvinfo}
        bash run_VSD_report.sh ${sampleid}_velvet_${minlen}-${maxlen}nt_blastn_vs_localdb.bls ${params.ictvinfo}
        """
    }
}
process filter_n_cov {
    tag "$sampleid"
    publishDir "${params.outdir}/07_filternstats/${sampleid}", mode: 'link'
    containerOptions "${bindOptions}"
    
    input:
    tuple val(sampleid), file(rawfastqfile), val(minlen), val(maxlen) from filter_n_cov_ch
    tuple val(sampleid), file(fastq_filt_by_size) from fastq_filt_by_size_ch
    tuple val(sampleid), file(samplefile), file(taxonomy) from blastTools_results_ch

    output:
    //file "${sampleid}_${minlen}-${maxlen}nt*.vcf.gz*"
    //file "${sampleid}_${minlen}-${maxlen}nt*consensus.fasta"
    //file "${sampleid}_${minlen}-${maxlen}nt*.fa"
    //file "${sampleid}_${minlen}-${maxlen}*.txt"
    file "${sampleid}_${minlen}-${maxlen}*"
    file "${sampleid}_${minlen}-${maxlen}nt_top_scoring_targets_with_cov_stats.txt" into contamination_flag
    
    script:
    """
    if [[ ${params.targets} == true ]]; then
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${rawfastqfile} --fastqfiltbysize  ${fastq_filt_by_size} --results ${samplefile} --read_size ${minlen}-${maxlen}nt --cov --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name} --targets --targetspath ${projectDir}/bin/${params.targets_file}
    else
        filter_and_derive_stats.py --sample ${sampleid} --rawfastq ${rawfastqfile} --fastqfiltbysize ${fastq_filt_by_size} --results ${samplefile} --read_size ${minlen}-${maxlen}nt --cov --taxonomy ${taxonomy} --blastdbpath ${blastn_db_name}
    fi
    """
}

if (params.contamination_detection) {
    process contamination_detection {
        label "local"
        publishDir "${params.outdir}/08_summary", mode: 'link'
        
        input:
        tuple val(sampleid), file(fastqfile), val(minlen), val(maxlen) from contamination_detection_ch
        file ('*') from contamination_flag.collect().ifEmpty([])

        output:
        file "run_top_scoring_targets_with_cov_stats_with_cont_flag*.txt"

        script:
        """
        flag_contamination.py --read_size ${minlen}-${maxlen}nt --threshold ${params.contamination_flag} --method ${params.contamination_detection_method}
        """
    }
}

if (params.blastp) {
    process getorf {
        label 'local'
        publishDir "${params.outdir}/06_blastp/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(velvet_cap3_rename_fasta), val(minlen), val(maxlen) from getorf_ch
        
        output:
        file "${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta"
        file "${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.circular.min50aa.fasta"
        tuple val(sampleid), file("${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta"), file("${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta.ids"), val(minlen), val(maxlen) into blastp_ch
        
        script:
        """
        getorf -sequence ${velvet_cap3_rename_fasta} -outseq ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta  -minsize ${params.orf_minsize}
        cat ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta.ids | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.min50aa.fasta.ids
        getorf -sequence ${velvet_cap3_rename_fasta} -circular True -outseq ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.circular.min50aa.fasta -minsize ${params.orf_circ_minsize}
        cat ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.circular.min50aa.fasta | grep ">" | sed 's/>//' | awk '{print \$1}' > ${sampleid}_velvet_${minlen}-${maxlen}nt_getorf.circular.min50aa.fasta.ids
        """
    }

    process blastp {
        label "xlarge"
        publishDir "${params.outdir}/06_blastp/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(fasta), file(fasta_ids), val(minlen), val(maxlen) from blastp_ch
        
        output:
        file "${fasta.baseName}_blastp_vs_NR_out.bls"
        tuple val(sampleid), file("${fasta.baseName}_blastp_vs_NR_out.wanted.ids"), file("${fasta.baseName}_blastp_vs_NR_out.bls"), val(minlen), val(maxlen) into blastpdbcmd_ch
        
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
        tuple val(sampleid), file(blastp_nr_bls_ids), file(blastp_nr_bls), val(minlen), val(maxlen) from blastpdbcmd_ch

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
        tuple val(sampleid), file(samplefile), val(minlen), val(maxlen) from spades_ch

        output:
        file "${sampleid}_${minlen}-${maxlen}nt_spades/*"
        file "${sampleid}_spades_assembly_${minlen}-${maxlen}nt.fasta"
        tuple val(sampleid), file("${sampleid}_spades_assembly_${minlen}-${maxlen}nt.fasta"), val(minlen), val(maxlen) into cap3_spades_ch
        
        script:
        """
        spades.py -t 2 -k ${params.spadeskmer} --only-assembler -m 180 -s $samplefile -o ${sampleid}_${minlen}-${maxlen}nt_spades
        
        #edit contigs name and rename velvet assembly
        sed 's/>/>spades_/' ${sampleid}_${minlen}-${maxlen}nt_spades/contigs.fasta > ${sampleid}_spades_assembly_${minlen}-${maxlen}nt.fasta
        """
    }

    process cap3_spades {
        label "local"
        publishDir "${params.outdir}/03_cap3/${sampleid}", mode: 'link'
        tag "$sampleid"

        input:
        tuple val(sampleid), file(scaffolds_fasta), val(minlen), val(maxlen) from cap3_spades_ch

        output:
        file "${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.fasta"
        tuple val(sampleid), file("${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.fasta"), file("${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.ids"), val(minlen), val(maxlen) into blastn_nt_spades_ch

        script:
        """
        cap3 ${scaffolds_fasta} -s 300 -j 31 -i 30 -p 90 -o 16
        cat ${scaffolds_fasta}.cap.singlets ${scaffolds_fasta}.cap.contigs > ${sampleid}_${minlen}-${maxlen}nt_spades_cap3.fasta
        extract_seqs_rename.py ${sampleid}_${minlen}-${maxlen}nt_spades_cap3.fasta ${params.cap3_len} | sed "s/CONTIG/${sampleid}_${minlen}-${maxlen}_/" | sed 's/|>/ |/' | awk '{print \$1}' > ${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.fasta
        #fetch scaffold IDs
        cat ${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.fasta | grep ">" | sed 's/>//' > ${sampleid}_spades_cap3_${minlen}-${maxlen}nt.rename.ids
        """
    }

    process blastn_nt_spades {
        label "medium_mem"
        publishDir "${params.outdir}/04_blastn/${sampleid}", mode: 'link'
        tag "$sampleid"
        containerOptions "${bindOptions}"

        input:
        tuple val(sampleid), file(spades_cap3_rename_fasta), file(spades_cap3_rename_fasta_ids), val(minlen), val(maxlen) from blastn_nt_spades_ch

        output:
        file "${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT.bls"
        file "${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits.txt"
        file "${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt"
        tuple val(sampleid), file("${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt") into BlastTools_blastn_spades_ch
        script:
        def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
        """
        blastn ${blast_task_param} \
                -query ${spades_cap3_rename_fasta} \
                -db ${blastn_db_name} \
                -out ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT.bls \
                -evalue ${params.blastn_evalue} \
                -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
                -num_threads ${task.cpus} \
                -max_target_seqs 100 \
                -word_size 11

        #fetch top blastn hits
        for i in `cat ${spades_cap3_rename_fasta_ids}`; do
            grep \$i ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT.bls | head -n5 >> ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits.txt;
        done
        grep -i "Virus" ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits.txt > ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        grep -i "Viroid" ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits.txt >> ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt || [[ \$? == 1 ]]
        cat ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids.txt | sed 's/ /_/g' > ${sampleid}_spades_${minlen}-${maxlen}nt_megablast_vs_NT_top5Hits_virus_viroids_final.txt
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
