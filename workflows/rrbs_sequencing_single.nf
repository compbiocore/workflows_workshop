#!/usr/bin/ nextflow
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
      module load trimmomatic/0.39-w5jnhai
      trimmomatic SE ${read} ${read.getBaseName()}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
    """
}

process fastqc {
  input:
    path(read)

  publishDir "${params.out_dir}/fastqc_original_reads", mode: 'copy'

  output:
    path "*"

  script:
    """
      module load libnsl/2.0.1-ed2i5hn
      module load fastqc/0.11.9-mvd2uhw
      fastqc ${read}
    """
   
}

process fastqc_trimmed {

  input:
    path(trimmed_read)

  publishDir "${params.out_dir}/fastqc_trimmed_reads", mode: 'copy'

  output:
    path "*"

  script:
    """
      module load fastqc/0.11.9-mvd2uhw
      module load libnsl/2.0.1-ed2i5hn
      fastqc ${trimmed_read}
    """
}

process alignment {
  
  input:
    path(trimmed_read)

  publishDir "${params.out_dir}/alignments", mode: 'copy'

  output:
   path "*.bam"
   path "*_report.txt"
   path "*_unmapped_reads.fq.gz"
  
  memory '10.GB' 

  script:
   """
    module load bismark/0.23.0-eoksupu
    module load perl/5.36.0-bt34quz 
    bismark -o `pwd` --bowtie2 --genome /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index --un --pbat ${trimmed_read}
   """  

} 

workflow
{
  fastqc(file(params.read))
  trimmed_result = trimmomatic(file(params.read))
  fastqc_trimmed(trimmed_result)
  alignment(trimmed_result)
}
