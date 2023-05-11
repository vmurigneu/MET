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

/* 
- Will I be mapping to the adaptive and non-adaptive files separately?
- Will I directly pass the minimap output into samtools etc. for sorting? 
- Will I need two separate processes for each adaptive non_adaptive subset?  
*/ 

Channel.fromPath
    adaptive_fq = Channel.fromPath()
    non_adaptive_fq = Channel.fromPath()
    genome = Channel.fromPath( '/data/refgenomes/*.{fa,fasta,fna}', checkIfExists: true ) 

params.adaptive_reads = "path to adaptive" 
params.nonadaptive_reads = "path to nonadaptive"


process minimap {
    cpus "${params.minimap_threads}"
    tag "${samples}"
    label "cpu"
    publishDir "$params.outdir/$sample/"3_mapping", mode: 'copy' pattern:
    "*.log", saveAs: { filename -> "${sample}_$filename" }  
    publishDir "$params.outdir/$sample/"3_mapping", mode: 'copy' pattern:
    "*_version.txt" 
    publishDir "$params.outdir/$sample/"3_mapping", mode: 'copy', pattern:
    '*fastq.gz', saveAs: { filename -> "${sample}.${filename" } 
    
    input:
        tuple val(sample), file(fastq_adaptive), file(fastq_non_adaptive) 
    output:
        tuple val(sample), file("adaptive_bac.fastq"),
        file("non_adaptive_bac.fastq"), emit: bacterial_fastq  
        path("minimap.log")   
        path("*txt")   
        path("*fastq")   
    when:
    !params.skip_remove_human_reads 
    shell:
    script:
    '''
    set +eu
    minimap2 -i !{   -t ${params.threads} -o *.sam 
    cp .command.log minimap2.log 
    minimap2 --version > minimap2_version.txt
    minimap2 -i 
    ''' 

// Parameters for flye 
Assembly:
		--flye_args				Flye optional parameters (default="--plasmids")
		--flye_threads 				Number of threads for Flye (default=4)

process flye {
	cpus "${params.flye_threads}"
	tag "${sample}"
	label "cpu"
    label "big_mem" 
	publishDir "$params.outdir/$sample/4_flye",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_flye",  mode: 'copy', pattern:
    "assembly*", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_flye",  mode: 'copy', pattern: "*txt", saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), file(fastq_adaptive_bac), file(fastq_non_adaptive_bac)
	output:
		tuple val(sample), file(""), path("adaptive_assembly_bac.fasta"),
        path("adaptive_assembly_info_bac.txt"),
        path("adaptive_assembly_graph_bac.gfa"),
        path("adaptive_assembly_graph_bac.gv"), path("non_adaptive_assembly_bac.fasta"),
        path("non_adaptive_assembly_info_bac.txt"), path("non_adaptive_assembly_graph_bac.gfa"),
        path("non_adaptive_assembly_graph_bac.gv")
        file("non_adaptive_bac.fasta"), emit: bacterial_assembly
		path("flye.log")
		path("flye_version.txt")
	
    when:
	!params.skip_metagenome_assembly
	shell:
	'''
	set +eu
	singularity exec /scratch/project/gihcomp/sw/flye_2.9.1--py38hf4f3596_0.sif
flye --nano-raw ${type}.fastq --threads !{params.metaflye_threads} --out-dir  
	'''
    
    for type in adaptive non_adaptive; do
		samtools sort -o ${type}.bam -@ !{params.minimap_threads} ${type}.sam
		samtools index ${type}.bam 
		samtools flagstat ${type}.bam > ${type}.flagstat.txt
		samtools view -S -f 4 -b ${type}.bam -o ${type}_non_human.unsorted.bam
		samtools sort -o ${type}_non_human.bam -@ !{params.minimap_threads} ${type}_non_human.unsorted.bam
  		samtools index ${type}_non_human.bam
 		samtools flagstat ${type}_non_human.bam > ${type}_non_human.flagstat.txt
        	samtools view ${type}_non_human.bam | cut -f1 | sort | uniq > ${type}_non_human_readID.lst
	done
	seqtk subseq !{fastq_adaptive} adaptive_non_human_readID.lst > adaptive_bac.fastq
	seqtk subseq !{fastq_non_adaptive} non_adaptive_non_human_readID.lst > non_adaptive_bac.fastq
	cp .command.log minimap.log
	'''
}



