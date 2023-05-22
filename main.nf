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
	
	Porechop: 
        	--porechop_args				Porechop optional parameters (default=""), see https://github.com/rrwick/Porechop#full-usage
		--porechop_threads			Number of threads for Porechop (default=4) (default=4)
		--skip_porechop				Skip the Porechop trimming step (default=false)
    
	Adaptive Read Sequencing: 
	
	Mapping: 
        	--minimap_threads			Number of threads for Porechop (default=12)
    
	Flye Assembly: 
        	--flye_threads          		Number of threads for Flye (default=?)
        	--memory                		Memory usage for Flye (default=0)
		--skip_assembly				Skipt the metagenome assembly (default=false)

	Centrifuge taxonomy classification:
		--skip_download_centrifuge_db		Skip the centrifuge database downloading step (default=false)
		--skip_centrifuge			Skip the centrifuge taxonomy classification step (default=false)
		--centrifuge_threads			Number of threads for Centrifuge (default=12)

	Polishing:
		--skip_polishing			Skip the Racon and Medaka polishing step (default=false)
		--racon_nb				Number of Racon long-read polishing iterations (default=4)
		--racon_args 				Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500")
		--racon_threads 			Number of threads for Racon (default=4)
		--medaka_threads 			Number of threads for Medaka (default=4)
		--medaka_model				Medaka model (default=r1041_e82_400bps_sup_g615)
	
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
	label "high_memory"
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*_version.txt"
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: '*fastq.gz', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(reads), path(csv)
	output:
		tuple val(sample), path("trimmed.fastq.gz"), path(csv),  emit: trimmed_fastq
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
        label "high_memory"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads), path(csv)
        output:
                tuple val(sample), path(reads), path("adaptive_reads.txt"), path("non_adaptive_reads.txt"), emit: extracted_readID
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
	label "high_memory"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*_version.txt"
        publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: '*fastq', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads), path(readID_adaptive), path(readID_nonadaptive)
        output:
		tuple val(sample),path("adaptive.fastq"),path("non_adaptive.fastq"), emit: extracted_fastq
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

process minimap {
	cpus "${params.minimap_threads}"
	tag "${sample}"
	label "high_memory"
	label "cpu"
	publishDir "$params.outdir/$sample/3_minimap",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_minimap",  mode: 'copy', pattern: "*fastq", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_minimap",  mode: 'copy', pattern: "*txt", saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(fastq_adaptive), path(fastq_non_adaptive)
	output:
		tuple val(sample), path("adaptive_bac.fastq"), path("non_adaptive_bac.fastq"), emit: bacterial_fastq
		path("minimap.log")
		path("*fastq")
		path("*txt")
	when:
	!params.skip_remove_human_reads
	shell:
	'''
	set +eu
	module load samtools/1.13-gcc-10.3.0 seqtk/1.3-gcc-10.3.0
	/scratch/project/gihcomp/sw/minimap2/minimap2 -t !{params.minimap_threads} -ax map-ont !{params.ref_human} !{fastq_non_adaptive} > non_adaptive.sam
	/scratch/project/gihcomp/sw/minimap2/minimap2 -t !{params.minimap_threads} -ax map-ont !{params.ref_human} !{fastq_adaptive} > adaptive.sam
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

prefix="flye"
prefix_lr="flye_polished"
raconv="racon"
medakav="medaka"

process flye {
	cpus "${params.flye_threads}"
	tag "${sample}"
	label "high_memory" 
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "adaptive_assembly*", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "non_adaptive_assembly*", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "*txt", saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(fastq_non_adaptive_bac)
	output:
		tuple val(sample), path(fastq_adaptive_bac), path("adaptive_assembly_bac.fasta"), path(fastq_non_adaptive_bac), path("non_adaptive_assembly_bac.fasta"), emit: bacterial_assembly_fasta
		tuple val(sample), path("adaptive_assembly_info_bac.txt"), path("adaptive_assembly_graph_bac.gfa"),path("adaptive_assembly_graph_bac.gv"), path("non_adaptive_assembly_info_bac.txt"), path("non_adaptive_assembly_graph_bac.gfa"),path("non_adaptive_assembly_graph_bac.gv"), emit: bacterial_assembly_graph
		path("flye.log")
		path("flye_version.txt")
	when:
	!params.skip_assembly
	shell:
	'''
	set +eu
	flye --nano-hq !{fastq_adaptive_bac} --threads !{params.flye_threads} --out-dir \$PWD !{params.flye_args}
	mv assembly.fasta adaptive_assembly_bac.fasta
	mv assembly_info.txt adaptive_assembly_info_bac.txt
	mv assembly_graph.gfa adaptive_assembly_graph_bac.gfa
	mv assembly_graph.gv adaptive_assembly_graph_bac.gv
	flye --nano-hq !{fastq_non_adaptive_bac} --threads !{params.flye_threads} --out-dir \$PWD !{params.flye_args} 
	mv assembly.fasta non_adaptive_assembly_bac.fasta
	mv assembly_info.txt non_adaptive_assembly_info_bac.txt
	mv assembly_graph.gfa non_adaptive_assembly_graph_bac.gfa
	mv assembly_graph.gv non_adaptive_assembly_graph_bac.gv
	flye -v 2> flye_version.txt
	cp .command.log flye.log
	'''  
}

process racon {
	cpus "${params.racon_threads}"
	tag "${sample}"
	label "racon"
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_$filename"}
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(adaptive_assembly), path(fastq_non_adaptive_bac), path(non_adaptive_assembly) 
	output:
		tuple val(sample), path(fastq_adaptive_bac), path("adaptive_${prefix}_${raconv}_${params.racon_nb}.fasta"), path(fastq_non_adaptive_bac), path("non_adaptive_${prefix}_${raconv}_${params.racon_nb}.fasta"), emit: polished_racon
	path("racon.log")
	path("racon_version.txt")
	when:
	!params.skip_polishing
	script:
	"""
	set +eu
	ln -s ${adaptive_assembly} adaptive_${prefix}_${raconv}_0.fasta
	for i in `seq 1 ${params.racon_nb}`; do
 		ii=\$((\$i-1))
		minimap2 -t ${params.racon_threads} -ax map-ont adaptive_${prefix}_${raconv}_\$ii.fasta ${fastq_adaptive_bac} > adaptive_${prefix}.gfa\$i.sam
		racon ${params.racon_args} -t ${params.racon_threads} ${fastq_adaptive_bac} adaptive_${prefix}.gfa\$i.sam adaptive_${prefix}_${raconv}_\$ii.fasta --include-unpolished > adaptive_${prefix}_${raconv}_\$i.fasta
		rm adaptive_${prefix}.gfa\$i.sam
	done
	ln -s ${non_adaptive_assembly} non_adaptive_${prefix}_${raconv}_0.fasta
	for i in `seq 1 ${params.racon_nb}`; do
		ii=\$((\$i-1))
		minimap2 -t ${params.racon_threads} -ax map-ont non_adaptive_${prefix}_${raconv}_\$ii.fasta ${fastq_non_adaptive_bac} > non_adaptive_${prefix}.gfa\$i.sam
		racon ${params.racon_args} -t ${params.racon_threads} ${fastq_non_adaptive_bac} non_adaptive_${prefix}.gfa\$i.sam non_adaptive_${prefix}_${raconv}_\$ii.fasta --include-unpolished > non_adaptive_${prefix}_${raconv}_\$i.fasta
		rm non_adaptive_${prefix}.gfa\$i.sam	
	done
	cp .command.log racon.log
	racon --version > racon_version.txt 
	"""
}

process medaka {
	cpus "${params.medaka_threads}"
	tag "${sample}"
	label "medaka"
	
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_$filename"}
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_assembly",  mode: 'copy', pattern: "*_version.txt" 
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(adaptive_draft), path(fastq_non_adaptive_bac), path(non_adaptive_draft)
	output:
		tuple val(sample), path(adaptive_draft), path ("adaptive_flye_polished.fasta"), path(non_adaptive_draft), path ("non_adaptive_flye_polished.fasta"), emit: polished_medaka
	path("medaka.log")
	path("medaka_version.txt")
	when:
	!params.skip_polishing	
	script:
	"""
	set +eu
	medaka_consensus -i ${fastq_adaptive_bac} -d ${adaptive_draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
	rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
	cp consensus.fasta adaptive_flye_polished.fasta
	medaka_consensus -i ${fastq_adaptive_bac} -d ${non_adaptive_draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
	rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
	cp consensus.fasta non_adaptive_flye_polished.fasta 
	cp .command.log medaka.log
	medaka --version > medaka_version.txt
 	"""
}

process centrifuge_download_db {
	cpus 1
	label "high_memory"
	publishDir "$params.outdir/centrifuge_database",  mode: 'copy', pattern: "*.cf"
	input:
		val(db)		
	output:
		tuple path("*.1.cf"), path("*.2.cf"), path("*.3.cf"), path("*.4.cf"), emit: centrifuge_db
	when:
	!params.skip_download_centrifuge_db
	script:
	"""
	echo ${db}
	wget ${db}
	tar -xvf nt_2018_3_3.tar.gz
	"""
}

process centrifuge {
	cpus "${params.centrifuge_threads}"
	tag "${sample}"
	label "very_high_memory"
	publishDir "$params.outdir/$sample/4_centrifuge",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_centrifuge",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(fastq_non_adaptive_bac), path(db1), path(db2), path(db3), path(db4)
	output:
		tuple val(sample), path(fastq_adaptive_bac), path(fastq_non_adaptive_bac), emit: bacterial_fastq
		tuple path("adaptive_centrifuge_report.tsv"), path("adaptive_centrifuge_species_report.tsv"), path("non_adaptive_centrifuge_report.tsv"), path("non_adaptive_centrifuge_species_report.tsv"), emit: centrifuge_reports
		path("centrifuge.log")
	when:
	!params.skip_centrifuge
	script:
	"""
	centrifuge -x nt -U ${fastq_adaptive_bac} -S adaptive_centrifuge_species_report.tsv --report-file adaptive_centrifuge_report.tsv --threads ${params.centrifuge_threads}
	centrifuge -x nt -U ${fastq_non_adaptive_bac} -S non_adaptive_centrifuge_species_report.tsv --report-file non_adaptive_centrifuge_report.tsv --threads ${params.centrifuge_threads}
	cp .command.log centrifuge.log
	"""
}

workflow {
	ch_centrifuge_db=Channel.value( "${params.centrifuge_db}")
	ch_centrifuge_db.view()
	Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
	.splitCsv(header:true, sep:',')
	.map { row -> tuple(row.sample_id, file(row.fastq, checkIfExists: true), file(row.csv, checkIfExists: true)) }
	.set { ch_samplesheet }
	ch_samplesheet.view()
	if (!params.skip_porechop) {
		porechop(ch_samplesheet)
		extract_adaptive_readID(porechop.out.trimmed_fastq)
		extract_adaptive_fastq(extract_adaptive_readID.out.extracted_readID)
		minimap(extract_adaptive_fastq.out.extracted_fastq)
		if (!params.skip_assembly) {
			flye(minimap.out.bacterial_fastq)
			if (!params.skip_polishing) {
				racon(flye.out.bacterial_assembly_fasta)
				medaka(racon.out.polished_racon)
			}
		}
		if (!params.skip_centrifuge) {
			if (!params.skip_download_centrifuge_db) {
				centrifuge_download_db(ch_centrifuge_db)
				centrifuge(minimap.out.bacterial_fastq.combine(centrifuge_download_db.out.centrifuge_db))
			} else if (params.skip_download_centrifuge_db) {
				ch_centrifuge_db=Channel.fromPath( "${params.outdir}/centrifuge_database/*.cf" ).collect()
				centrifuge(minimap.out.bacterial_fastq.combine(ch_centrifuge_db))
			}
		}
	} else if (params.skip_porechop) {	
		extract_adaptive_readID(ch_samplesheet)
		extract_adaptive_fastq(extract_adaptive_readID.out.extracted_readID)
		minimap(extract_adaptive_fastq.out.extracted_fastq)
	}
}
