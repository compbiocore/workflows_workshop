#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params.name = "World"
params.count_letters = false
params.reverse = false

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

process count {
  input: 
    path(file_in)
  output:
    stdout
  
  script:
   if (params.count_letters) {
    """
    cat ${file_in}
    wc -m ${file_in} | awk '{print \$1}'
    """
   } else {
    """
    cat ${file_in}
    wc -w ${file_in} | awk '{print \$1}'
    """
   }
}

process reverse {
  input: 
    path(file_in)
  output:
    path "reverse_hello.txt"
  
  script:
    """
    cat ${file_in} | rev > "reverse_hello.txt"
    """
}

workflow {
  hello_result = null

  if (params.reverse) {
    hello_result = reverse(sayHello(params.name))
  } else {
    hello_result = sayHello(params.name)
  }

  count(hello_result) | view
}
