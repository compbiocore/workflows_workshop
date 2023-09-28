#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.name = "World"

process sayHello {
  input:
    val name
  output:
    tuple path("hello.txt"), stdout
  script:
    """
    echo 'Hello ${name}!' > hello.txt
    echo $( stat -c '%s' hello.txt ) * 3 ### get the size of the output file and multiply it by 3
    """
}

process countWords {
  input:
    tuple path(file_in), val(file_size)
  output:
    path("count_words.txt")
  memory "${file_size}"

  publishDir "${params.out_dir}/", mode: 'copy'

  script:
   """
   wc -w ${file_in} | awk '{print \$1}' > count_words.txt
   """
}

workflow {
  countWords(sayHello(params.name))
}