#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
         MET analysis pipeline
========================================================================================
 #### Documentation
 #### Authors
 Valentine Murigneux <v.murigneux@uq.edu.au>
 Dean Basic <d.basic@uq.edu.au>
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	MET analysis pipeline v${workflow.manifest.version}
	=========================================
	Usage:
	nextflow main.nf --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/

	Required arguments:
		--fqdir					Path to the directory containing the PacBio subreads bam files
		--outdir				Path to the output directory to be created
    
	Optional parameters:
		--threads				Number of threads (default=16)
		--porechop_args				Porechop optional parameters (default=""), see https://github.com/rrwick/Porechop#full-usage
		--porechop_threads			Number of threads for Porechop (default=4) (default=4)
		--skip_porechop				Skip the Porechop trimming step (default=false)
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


process porechop {
	cpus "${params.porechop_threads}"
	tag "${sample}"
	label "cpu"
	label "big_mem"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*_version.txt"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: '*fastq.gz', saveAs: { filename -> "${sample}.$filename" }
	input:
		tuple val(sample), file(reads)
	output:
		tuple val(sample), file("trimmed.fastq.gz"), emit: trimmed_fastq
		path("porechop.log")
		path("porechop_version.txt")
		path("*fastq.gz")
	when:
	!params.skip_porechop
	script:
	"""
	set +eu
	porechop -i ${reads} -t ${params.porechop_threads} -o trimmed.fastq.gz ${params.porechop_args}
	cp .command.log porechop.log
	porechop --version > porechop_version.txt
	"""
}

workflow {
	ch_fastq=Channel.fromPath( "${params.fastqdir}/*.fastq.gz" ). map { file -> tuple(file.simpleName, file) } 
	ch_fastq.view()	
	porechop(ch_fastq)
}

process minimap2 {
    cpus "${params.threads}"
    tag "${minimap2}"
    label "cpu"
    label "big_mem"
    publishDir "$params.outdir/$minimap2/
    
    input:
        tuple val(), file()
    output:
        tuple val()

    script:
    """
    set +eu
    minimap2 -i $   -t ${params.threads} -o *.sam 
    cp .command.log minimap2.log 
    minimap2 --version > minimap2_version.txt 
    









