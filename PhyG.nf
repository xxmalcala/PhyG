#!/usr/bin/env nextflow

/*
* pipeline input parameters
*/
params.reads = ""
params.out_dir = ""
params.num_cpus = 4
params.rrna_outdir = "$projectDir/$params.out_dir/rRNA_Reads/"

log.info """


        P h y G - P i p e l i n e
   ===================================
   projdir      : ${projectDir}
   out_dir      : ${params.out_dir}
   reads        : ${params.reads}
   threads      : ${params.num_cpus}
   ===================================


   """
   .stripIndent(true)

/*
* Use RiboDetector to extract putative rRNA reads from a set of RNA-seq datasets


process DETECT_RRNA_READS {
  tag "RiboDector"
  publishDir params.rrna_outdir, mode: 'copy'
  cpus = params.num_cpus

  input:
  path reads

  script:
  """
  ribodetector -t $task.cpus -i ${reads[0]} ${reads[1]} -e both -o ${outdir}.FPE.Non_rRNA.fastq.gz ${outdir}.RPE.Non_rRNA.fastq.gz -r ${outdir}.FPE.rRNA.fastq.gz ${outdir}.RPE.rRNA.fastq.gz
  """

}
*/
