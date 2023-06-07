#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params.name = "World"

process sayHello {
  input: 
    val name
  output:
    stdout
  script:
    """
    echo 'Hello ${name}!'
    """
}

workflow {
  sayHello(params.name) | view
}
