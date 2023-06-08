#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process count_document {
  input: 
    tuple file(file_in), val(count_letters)
  output:
    file "${file_in.baseName}_count.txt"
  
  script:
   if (count_letters) {
    """
    wc -m ${file_in} | awk '{print \$1}' > ${file_in.baseName}_count.txt
    """
   } else {
    """
    wc -w ${file_in} | awk '{print \$1}' > ${file_in.baseName}_count.txt
    """
   }
}

process sum_all_results {
  input: 
    path (files)
  output:
    stdout

  script:
    """
    cat ${files} > combined.txt
    awk '{ sum += \$1 } END { print sum }' combined.txt
    """
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
        word_or_letter_counts = count_document(input_ch)
    
    emit:
        word_or_letter_counts
}

workflow {
     // create a Channel of documents to count either by word or by letter, from the samplesheet
     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .map { get_document_info(it) }.set { documents_ch }
    
    // launch a sub-workflow of COUNT_DOCUMENT on each document
    count_results = COUNT_DOCUMENT(documents_ch).collect()
    sum_all_results(count_results) | view
}