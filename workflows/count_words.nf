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
    stdout
  
  script:
   """
   cat ${file_in}
   wc -w ${file_in} | awk '{print \$1}'
   """ 
}

workflow {
  countWords(sayHello(params.name)) | view
}
