#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process count_document {
  input: 
    tuple file(file_in), val(count_letters)
  output:
    stdout
  
  script:
   if (count_letters) {
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

// convert each row in the samplesheet into a tuple of File object and a Boolean of count_letter
def get_document_info(LinkedHashMap sample) {
    document  = sample.document ? file(sample.document, checkIfExists: true) : null
    count_letters = (sample.count_letters == "T") || (sample.count_letters == "true") ? true : false   

    return [ document, count_letters ]
}

workflow COUNT_DOCUMENT {
    take:
        input_ch

    main:
        count_document(input_ch) | view
}

workflow {
     // create a Channel of documents to count either by word or by letter, from the samplesheet
     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .map { get_document_info(it) }.set { documents_ch }
    
    // launch a sub-workflow of COUNT_DOCUMENT on each document
     COUNT_DOCUMENT(documents_ch)
}