#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params.read = "/gpfs/data/cbc/workflow_workshop/sample1.fq.gz"

process trimmomatic {
  
  input: 
    path(read)

  publishDir "${params.out_dir}/trimmed_reads", mode: 'copy'
  
  output:
    path("*_tr.fq.gz")
  
  script:
    """
      module load trimmomatic/0.39
      TrimmomaticSE ${read} ${read.getBaseName()}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
    """
}

process fastqc {
  input:
    path(read)

  publishDir "${params.out_dir}/fastqc_original_reads", mode: 'copy'

  container 'biocontainers/fastqc:v0.11.9_cv8'

  output:
    path "*"

  script:
    """
      fastqc ${read}
    """
   
}

process fastqc_trimmed {

  input:
    path(trimmed_read)

  publishDir "${params.out_dir}/fastqc_trimmed_reads", mode: 'copy'

  container 'biocontainers/fastqc:v0.11.9_cv8'

  output:
    path "*"

  script:
    """
      fastqc ${trimmed_read}
    """
}

process alignment {
  
  input:
    path(trimmed_read)

  publishDir "${params.out_dir}/alignments", mode: 'copy'

  output:
   path "*.sam.gz"
   path "*_report.txt"
   path "*_unmapped_reads.fq.gz"
   

  script:
   """
    module load bismark/0.20.0
    module load bowtie2/2.4.2
    bismark -o `pwd` --bowtie2 --genome /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index --un --pbat ${trimmed_read}
   """  

} 

workflow
{
  reads = Channel.fromPath(params.read_dir + "/*.fq.gz", checkIfExists: true)
  fastqc(reads)
  trimmed_result = trimmomatic(reads)
  fastqc_trimmed(trimmed_result)
  alignment(trimmed_result)
}
