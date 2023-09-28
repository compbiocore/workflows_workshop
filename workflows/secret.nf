#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process sayHello {
  input:
    secret 'dscov_secret_location'
  output:
    stdout
  script:
    """
    curl https://api.weather.gov/alerts/active?area=\$dscov_secret_location
    """
}

workflow {
  sayHello() | view
}