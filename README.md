# miRNA_pipeline

The prepared pipeline is used for processing files after small-RNA-seq sequencing in order to perform differential expression analysis of miRNA. 
The pipeline was written in Nextflow and consists of the following programs:

1. Fastq
2. MultiQC
3. MiRTrace
4. ShortStack

The pipeline was prepared as part of the National Science Center project UMO-2020/37/B/NZ9/00586.

Instruction:

To run the pipeline, execute the following command:

nextflow run run_shortstack.nf --mature_miRNA_list path_to_miRNA_list --reference_genome path_to_reference_genome --fastq_files path_to_fastq_files_directory -with-report

The required programs are available as a conda environment, which can be created from a .yml file, or as a Docker container.