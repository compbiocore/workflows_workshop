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
  input: 
    path(file_in)
  output:
    path("count_words.txt")

  publishDir "${params.out_dir}/", mode: 'copy'
  
  script:
   """
   wc -w ${file_in} | awk '{print \$1}' > count_words.txt
   """ 
}

workflow {
  countWords(sayHello(params.name))
}
