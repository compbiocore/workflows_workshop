#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params.name = "World"

process sayHello {
  input: 
    val name
  output:
    path "hello.txt"
  script:
    """
    echo 'Hello ${name}!' > hello.txt
    """
}

process countWords {
  publishDir "${params.out_dir}/", mode: 'copy'

  input: 
    path(file_in)
  output:
    path("count_words.txt")

  script:
   """
   wc -w ${file_in} | awk '{print \$1}' > count_words.txt
   """ 
}

workflow {
  countWords(sayHello(params.name))
}
