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
		--skip_extract_adaptive			Skip the adaptive/non-adaptive read extraction step (default=false)
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
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*_version.txt"
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: '*fastq.gz', saveAs: { filename -> "${sample}.$filename" }
	input:
		tuple val(sample), file(reads), file(csv)
	output:
		tuple val(sample), file("trimmed.fastq.gz"), file(csv),  emit: trimmed_fastq
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

process extract_adaptive_readID {
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), file(reads), file(csv)
        output:
                tuple val(sample), file(reads), file("adaptive_reads.txt"), file("non_adaptive_reads.txt"), emit: extracted_readID
                path("extract_adaptive_readID.log")
                path("*txt")
        when:
        !params.skip_extract_adaptive
        shell:
        '''
        set +eu
	awk -F, '$7 = "no_decision" {print $0}' !{csv} | cut -d" " -f5 | tail -n +2 | sort | uniq > adaptive_reads.txt	
	seqkit fx2tab !{reads} | awk '{print $1, $5}' - | sed 's/=/ /' | cut -d" " -f1,3 | awk '$2 > 256 {print $1}' - | sort | uniq > non_adaptive_reads.txt   
	cp .command.log extract_adaptive_readID.log
        '''
}

process extract_adaptive_fastq {
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*_version.txt"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: '*fastq', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), file(reads), file(readID_adaptive), file(readID_nonadaptive)
        output:
                tuple val(sample), file("adaptive.fastq"), file("non_adaptive.fastq"), emit: extracted_fastq
                path("extract_adaptive_fastq.log")
                path("*fastq")
        when:
        !params.skip_extract_adaptive
        shell:
        '''
        set +eu
        seqtk subseq !{reads} !{readID_adaptive} > adaptive.fastq
        seqtk subseq !{reads} !{readID_nonadaptive} > non_adaptive.fastq
        cp .command.log extract_adaptive_fastq.log
        '''
}

workflow {
	//ch_fastq=Channel.fromPath( "${params.fastqdir}/*.fastq.gz" ). map { file -> tuple(file.simpleName, file) } 
	//ch_fastq.view()	
	Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
	.splitCsv(header:true, sep:',')
	.map { row -> tuple(row.sample_id, file(row.fastq, checkIfExists: true), file(row.csv, checkIfExists: true)) }
	.set { ch_samplesheet }
	ch_samplesheet.view()
	porechop(ch_samplesheet)
	extract_adaptive_readID(porechop.out.trimmed_fastq)
	extract_adaptive_fastq(extract_adaptive_readID.out.extracted_readID)
}
