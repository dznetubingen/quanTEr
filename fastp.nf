#!/usr/bin/env nextflow

/*
================================================================================
                           quanTEfication/main
================================================================================
  Nextflow pipeline for RNA-seq analysis and quantefication of the transposons.
  fastp.nf is used for preprocessing of raw .fq files before their alignment and
  quantification. Produces also multiQC report and report for read mates synch
  of the files that have been processed by fastp to ensure the order of reads is
  correct and maximize the mapping rate in the following steps.
  PLEASE CHECK BOTH REPORTS BEFORE THE MAPPING AND QUANTIFICATION STAGE!
--------------------------------------------------------------------------------
 @Repository
 https://github.com/savytskanatalia/quanTEfication
 @Author
 Savytska Natalia, 2020. Genome Biology of Neurodegenerative Diseases, DZNE TU
--------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run fastp.nf -c fastp.config --reads '/absolute/path/to/files/*_R{1,2}.fastq.gz' -with-docker savytskanatalia/quantefication2
    Mandatory arguments:
      --reads [str]                  Path to input data (must be surrounded with quotes). Default: 'input/*_R{1,2}*.fastq.gz'
      --outdir [str]                 Output directory path. Default: 'output/trimmed'
      --fastpf [str]                 Additional fastp flags for processing. Default additional options: ''. Default hardcoded options: '--detect_adapter_for_pe -w 1 -p -j sampleID.json -h  sampleID.html'.
                                     Example usage. To work with FOUNDIN RNAB data, add flag '-f 3 -F 0' to cut offtemplate switching oligo from the first sequecing read.
      --thr [int]                    Number of threads to allocate for fastp processing. Default: 1. Warning! Be cautious with multithreading! It results in not so deterministic outcome (e.g. read order in output).
    """.stripIndent()
}


threads=Channel.value(params.thr)
fastpopt=Channel.value(params.fastpf)
fasta_files=Channel.fromFilePairs(params.reads)
outdir=params.outdir









process fastp_trim {
  container='savytskanatalia/quantefication2'
  echo true
  publishDir "${params.outdir}/trimmed", mode: 'move', pattern: "*_R{1,2}_PROC.fastq.gz"
  publishDir "${params.outdir}/failed", mode: 'move', pattern: "*_failed.fastq.gz"
  publishDir "${params.outdir}/reports", mode: 'copy', pattern: "*fastp*"
  publishDir "${params.outdir}/reports", mode: 'copy', pattern: "*.json"
  input:
  set sampleId, file(reads) from fasta_files
  val i from threads
  val x from fastpopt
  output:
  file "${sampleId}_R{1,2}_PROC.fastq.gz"
  file "${sampleId}_failed.fastq.gz"
  file "*fastp.json" into multiqc_ch
  file "*fastp.html"
  file "${sampleId}.json" into synch_ch
  script:
  """
  echo HI THERE ${reads[0]} $sampleId
  fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleId}_R1_PROC.fastq.gz -O ${sampleId}_R2_PROC.fastq.gz --detect_adapter_for_pe -w $i -p -j ${sampleId}_fastp.json -h  ${sampleId}_fastp.html  --failed_out ${sampleId}_failed.fastq.gz $x
  pecheck -i ${sampleId}_R1_PROC.fastq.gz -I ${sampleId}_R2_PROC.fastq.gz -j ${sampleId}.json

  """
}

process multiqc_aggregation {
  container='savytskanatalia/quantefication2'
  echo true
  publishDir "${params.outdir}/reports", mode: 'copy'
  input:
  file ("fastp/*") from multiqc_ch.collect().ifEmpty([])
  output:
  file "multiqc*"
  script:
  """
  multiqc . -m fastp 

  """
}

process synch_check {
  container='savytskanatalia/quantefication2'
  echo true
  publishDir "${params.outdir}/reports", mode: 'move'
  input:
  file ("pecheck/*") from synch_ch.collect().ifEmpty([])
  output:
  file "synch.stat"
  shell:
  $/
  echo "Pair mates synchronicity check" >> synch.stat
  pwd >> synch.stat
  for f in pecheck/*.json; do echo "" >> synch.stat  && echo "$$f" | tr "\n" "\t" >> synch.stat && grep -Po '"result":.*?[^\\]"|"message":.*?[^\\]"' $$f | tr "\n" "\t" >> synch.stat ; done
  /$

}






