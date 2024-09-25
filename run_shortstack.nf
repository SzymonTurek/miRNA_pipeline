#!/usr/bin/env nextflow


params.mature_miRNA_list = "$baseDir/Mature.fa"
params.reference_genome = "$baseDir/pb_b10_ill1.fasta"
params.fastq_files = "$baseDir/sample_info.csv"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"
params.fastqc_outdir = "results/fastqc_output"
params.shortstack_output_dir = "results/shortstack_output"
//params.counts = "$baseDir/results/shortstack_output/Counts.txt"
//params.results = "$baseDir/results/shortstack_output/Results.txt"

log.info """\
Input parameters:
$params.mature_miRNA_list
$params.reference_genome
$params.fastq_files
$params.multiqc 
$params.outdir    
"""

log.info """\
    S H O R T S T A C K  P I P E L I N E
    ===================================
    transcriptome: ${params.mature_miRNA_list}
    reads: ${params.reference_genome}
    fastq_files: ${params.fastq_files}
    multiqc: ${params.multiqc}
    outdir: ${params.outdir}
    """
    .stripIndent(true)




process READ_FASTQS {
    debug true
    input:
    tuple val(sampleId), file(read1)

    script:
    """
    echo head --sample $sampleId --reads $read1 
    """
}

 process ECHO_CHANNEL {

        debug true
        input:
        tuple val(sampleid), file(file)

        script:
        """
        echo "File name: $sampleid"
        echo "File path: $file"
        """
    }

process FASTQC {
    tag "FastQC on $sample_id"
    publishDir params.fastqc_outdir, mode: 'copy'


    input:
     tuple val(sample_id), path(reads)
        

    output:
    path "fastqc_${sample_id}"

    script:
    """
    echo ${reads}
    echo ${sample_id}
    mkdir fastqc_${sample_id}
    fastqc -o fastqc_${sample_id} -f fastq -q ${reads}
    """


}


process MULTIQC {

    publishDir params.outdir, mode: 'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'


    script:
    """
    multiqc .
    """

}

process RUN_SHORTSTACK {
    tag "ShortStack on all samples"
    publishDir params.outdir, mode: 'copy'

    label 'big_mem'

    input:
     path(genome_file)
     path(known_mirs)
     path(reads)
       

    output:
    path "./shortstack_output"
    path "./shortstack_output/Counts.txt", emit: counts_txt
    path "./shortstack_output/Results.txt", emit: results_txt


    script:
    """
    ShortStack  --genomefile ${genome_file}  --known_miRNAs ${known_mirs}  --readfile ${reads} --outdir shortstack_output --threads 8
    """


}




//ShortStack  --genomefile ${genome_file}  --known_miRNAs ${known_mirs}  --readfile ${reads} --outdir shortstack_output --threads 7 --dn_mirna

process GREP_Y {
    tag "GREP output files"
    publishDir params.shortstack_output_dir, mode: 'copy'


    input:
     path(counts_file)
     path(results_file)
     
       

    output:
    path "Results_Y.txt"
    path "Counts_Y.txt"


    script:
    """
    grep "Y" ${counts_file} > Counts_Y.txt
    grep "Y" ${results_file} > Results_Y.txt
   
    """
}

process MIRTRACE {
    tag "RUN MIRTRACE ON ALL FILES"
    publishDir params.outdir, mode: 'copy'


    input:
     path(reads)
    
    output:
    path "miRTrace_output"

    script:
    """
    mirtrace qc -s meta_species_all -o miRTrace_output ${reads}
   
    """
}


workflow {
    //Channel.fromPath(params.fastq_files) \
     //   | splitCsv(header:true) \
      //  | map { row-> tuple(row.sample_id, file(row.read1)) } \
       // | READ_FASTQS
    

   // Channel.fromPath(params.fastq_files)
   //     .splitCsv( header: true )
   //     .view { row -> "${row.sample_id} - ${row.read1}" }

    samples_ch = Channel.fromPath(params.fastq_files)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.read1)) }
    //    .view()
    
    collected_samples_ch = Channel.fromPath(params.fastq_files)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.read1)) }
        .collect()
        .view()


   // samples_ch.view()
   //  ECHO_CHANNEL(samples_ch)

   fastqc_ch = FASTQC(samples_ch)
   MULTIQC(fastqc_ch.collect())
   shortstack_ch = RUN_SHORTSTACK(params.reference_genome, params.mature_miRNA_list, collected_samples_ch  )
   grep_ch = GREP_Y(RUN_SHORTSTACK.out.counts_txt, RUN_SHORTSTACK.out.results_txt)
   MIRTRACE(collected_samples_ch)
}
