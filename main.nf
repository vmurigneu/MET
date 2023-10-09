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
		--porechop_threads			Number of threads for Porechop (default=4) 
		--skip_porechop				Skip the Porechop trimming step (default=false)
    
	Adaptive Read Sequencing: 
		--skip_adaptive_sampling_metrics	Skip the Adaptive sampling metrics step (default=false)
		--nanocomp_threads			Number of threads for NanoComp (default=4)
	Mapping: 
        	--minimap_threads			Number of threads for Minimap2 (default=12)
    
	Flye Assembly: 
		--flye_args				Flye optional parameters (default=`--flye_args "-meta"), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
		--flye_threads          		Number of threads for Flye (default=?)
        	--memory                		Memory usage for Flye (default=0)
		--skip_assembly				Skip the metagenome assembly step (default=false)

	Centrifuge taxonomy classification:
		--skip_download_centrifuge_db		Skip the centrifuge database downloading step (default=false)
		--skip_centrifuge			Skip the centrifuge taxonomy classification step (default=false)
		--centrifuge_db				Path to download the centrifuge database (default='https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz')
		--centrifuge_reference_tax_ID		Taxonomy ID for reference ID (default="")
		--skip_centrifuge_remove_contaminated	Skip the removal of centrifuge reads from contaminated reference
		--centrifuge_threads			Number of threads for Centrifuge (default=12)
		--skip_krona				Skip the generation of Krona plots (default=false)
	
	Polishing:
		--skip_polishing			Skip the Racon and Medaka polishing step (default=false)
		--racon_nb				Number of Racon long-read polishing iterations (default=4)
		--racon_args 				Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500")
		--racon_threads 			Number of threads for Racon (default=4)
		--medaka_threads 			Number of threads for Medaka (default=4)
		--medaka_model				Medaka model (default=r1041_e82_400bps_sup_g615)

	Virus and Plasmid classification:
		--skip_download_genomad_db		Skip the genomad database download if it is already present locally (default=false)
		--skip_genomad				Skip genomad classification (default=false)

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
	awk -F, '$7 == "stop_receiving" {print $0}' !{csv} | cut -d"," -f5 | sort | uniq > adaptive_reads.txt	
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

process extract_adaptive_sampling_reads {
	tag "${sample}"
	label "cpu"
	label "high_memory"
	publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: "*_version.txt"
	publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: '*fastq', saveAs: { filename -> "${sample}_$filename" }
	input:	
		tuple val(sample), path(reads), path(csv)
	output:
		tuple val(sample),path("stopreceiving.fastq"),path("unblock.fastq"),path("nodecision.fastq"), emit: extracted_fastq
		path("extract_adaptive_fastq.log")
		path("*nodecision.fastq")
		path("unblock.fastq")
	when:
	!params.skip_adaptive_sampling_metrics
	shell:
	'''
	set +eu
	awk -F, '$7 == "stop_receiving" {print $0}' !{csv} | cut -d"," -f5 | sort | uniq | sort > readID_stop_receiving.txt
	awk -F, '$7 == "no_decision" {print $0}' !{csv} | cut -d"," -f5 | sort | uniq | sort > readID_no_decision.txt
	awk -F, '$7 == "unblock" {print $0}' !{csv} | cut -d"," -f5 | sort | uniq | sort > readID_unblock.txt
	comm -12 readID_no_decision.txt readID_unblock.txt > comm_nodecision_unblock.txt 
	awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' readID_no_decision.txt readID_unblock.txt > readID_unblock_not_no_decision.txt
	cat readID_unblock_not_no_decision.txt comm_nodecision_unblock.txt | sort | uniq > readID_unblock_unique.txt 
	awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' readID_unblock.txt readID_no_decision.txt > readID_unique_no_decision_unblock.txt
	awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' readID_stop_receiving.txt readID_unique_no_decision_unblock.txt > readID_no_decision_unique.txt
	seqtk subseq !{reads} readID_stop_receiving.txt > stopreceiving.fastq
	seqtk subseq !{reads} readID_no_decision_unique.txt > nodecision.fastq
	seqtk subseq !{reads} readID_unblock_unique.txt > unblock.fastq
	nb_reads=`cut -d"," -f5  !{csv} | tail -n +2 | sort | uniq | sort | wc -l`
	echo -e !{sample}\\t$nb_reads >> nbReads_AS_csv.txt
 	cp .command.log extract_adaptive_fastq.log
	'''
}

process compute_adaptive_sampling_metrics {
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/2_adaptive",  mode: 'copy', pattern: '*txt', saveAs: { filename -> "${sample}_$filename" }
	input:  
		tuple val(sample), path(fq_stopreceiving), path(fq_unblock), path(fq_nodecision)
	output:
		tuple val(sample),path("NanoStats.txt")
	when:
	!params.skip_adaptive_sampling_metrics
	shell:
	'''
	set +eu
	NanoComp -o \$PWD --fastq !{fq_stopreceiving} !{fq_unblock} !{fq_nodecision} -t !{params.nanocomp_threads} -n stop_receiving unblock no_decision
	grep "N50" NanoStats.txt | tr -d ' ' | sed s/ReadlengthN50:// | sed s/\\.0/\\t/g >> ReadN50.txt
 	grep "Number" NanoStats.txt| grep reads | grep -v percentage | tr -d ' ' | sed s/Numberofreads://| sed s/\\.0/\\t/g  >> NbReads.txt
	grep "Median" NanoStats.txt | grep length |tr -d ' ' | sed s/Medianreadlength:// | sed s/\\.0/\\t/g >> MedianReadLength.txt
	grep "Median" NanoStats.txt | grep quality |tr -s ' ' | sed s/^Median\\sread\\squality:// | sed s/\\s/\\t/g | sed s/^\\s//g >> MedianReadQuality.txt
	echo !{sample} >> samples.txt
	cp .command.log compute_adaptive_sampling_metrics.log
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
	!params.skip_remove_host_reads
	shell:
	'''
	set +eu
	module load samtools/1.13-gcc-10.3.0 seqtk/1.3-gcc-10.3.0
	/scratch/project/gihcomp/sw/minimap2/minimap2 -t !{params.minimap_threads} -ax map-ont !{params.ref_genome} !{fastq_non_adaptive} > non_adaptive.sam
	/scratch/project/gihcomp/sw/minimap2/minimap2 -t !{params.minimap_threads} -ax map-ont !{params.ref_genome} !{fastq_adaptive} > adaptive.sam
	for type in adaptive non_adaptive; do
		samtools sort -o ${type}.bam -@ !{params.minimap_threads} ${type}.sam
		samtools index ${type}.bam 
		samtools flagstat ${type}.bam > ${type}.flagstat.txt
		samtools view -S -f 4 -b ${type}.bam -o ${type}_non_host.unsorted.bam
		samtools sort -o ${type}_non_host.bam -@ !{params.minimap_threads} ${type}_non_host.unsorted.bam
  		samtools index ${type}_non_host.bam
 		samtools flagstat ${type}_non_host.bam > ${type}_non_host.flagstat.txt
        	samtools view ${type}_non_host.bam | cut -f1 | sort | uniq > ${type}_non_host_readID.lst
	done
	seqtk subseq !{fastq_adaptive} adaptive_non_host_readID.lst > adaptive_bac.fastq
	seqtk subseq !{fastq_non_adaptive} non_adaptive_non_host_readID.lst > non_adaptive_bac.fastq
	cp .command.log minimap.log
	'''
}

prefix="assembly"
prefix_lr="assembly_polished"
raconv="racon"
medakav="medaka"

process flye {
	cpus "${params.flye_threads}"
	tag "${sample}"
	label "high_memory" 
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "adaptive_assembly*", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "non_adaptive_assembly*", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "*txt", saveAs: { filename -> "${sample}_$filename" }
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
	if [ -f "assembly.fasta" ]; then
		mv assembly.fasta adaptive_assembly_bac.fasta
		mv assembly_info.txt adaptive_assembly_info_bac.txt
		mv assembly_graph.gfa adaptive_assembly_graph_bac.gfa
		mv assembly_graph.gv adaptive_assembly_graph_bac.gv
	else
		touch adaptive_assembly_bac.fasta adaptive_assembly_info_bac.txt adaptive_assembly_graph_bac.gfa adaptive_assembly_graph_bac.gv
	fi
	flye --nano-hq !{fastq_non_adaptive_bac} --threads !{params.flye_threads} --out-dir \$PWD !{params.flye_args} 
	if [ -f "assembly.fasta" ]; then
		mv assembly.fasta non_adaptive_assembly_bac.fasta
		mv assembly_info.txt non_adaptive_assembly_info_bac.txt
		mv assembly_graph.gfa non_adaptive_assembly_graph_bac.gfa
		mv assembly_graph.gv non_adaptive_assembly_graph_bac.gv
	else
		touch non_adaptive_assembly_bac.fasta non_adaptive_assembly_info_bac.txt non_adaptive_assembly_graph_bac.gfa non_adaptive_assembly_graph_bac.gv
	fi
	flye -v 2> flye_version.txt
	cp .command.log flye.log
	'''  
}

process racon {
	cpus "${params.racon_threads}"
	tag "${sample}"
	label "racon"
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_$filename"}
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "*_version.txt"
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
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_$filename"}
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_assembly",  mode: 'copy', pattern: "*_version.txt" 
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(adaptive_draft), path(fastq_non_adaptive_bac), path(non_adaptive_draft)
	output:
		tuple val(sample), path(fastq_adaptive_bac), path ("adaptive_flye_polished.fasta"), path(fastq_non_adaptive_bac), path ("non_adaptive_flye_polished.fasta"), emit: polished_medaka
	path("medaka.log")
	path("medaka_version.txt")
	when:
	!params.skip_polishing	
	script:
	"""
	set +eu
	medaka_consensus -i ${fastq_adaptive_bac} -d ${adaptive_draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
	rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
	if [ -f "consensus.fasta" ]; then
		mv consensus.fasta adaptive_flye_polished.fasta
	else
		touch adaptive_flye_polished.fasta
	fi
	medaka_consensus -i ${fastq_adaptive_bac} -d ${non_adaptive_draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
	rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
	if [ -f "consensus.fasta" ]; then
		mv consensus.fasta non_adaptive_flye_polished.fasta 
	else
		touch non_adaptive_flye_polished.fasta
	fi
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
		tuple val(sample), path(fastq_adaptive_bac), path("adaptive_centrifuge_species_report.tsv"), path(fastq_non_adaptive_bac), path("non_adaptive_centrifuge_species_report.tsv"), emit: bacterial_fastq
		tuple val(sample), path("adaptive_centrifuge_species_report.tsv"), path("non_adaptive_centrifuge_species_report.tsv"), emit: centrifuge_species_report
		tuple val(sample), path("adaptive_centrifuge_report.tsv"), path("non_adaptive_centrifuge_report.tsv"), emit: centrifuge_report
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

process remove_centrifuge_contaminated {
    tag "${sample}"
    //label "very_high_memory"
    publishDir "$params.outdir/$sample/5_centrifuge_bac_reads",  mode: 'copy', pattern: "*.lst", saveAs: { filename -> "${sample}_$filename" }
    publishDir "$params.outdir/$sample/5_centrifuge_bac_reads",  mode: 'copy', pattern: "*.fastq", saveAs: { filename -> "${sample}_$filename" }
    publishDir "$params.outdir/$sample/5_centrifuge_bac_reads",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
    input:
        tuple val(sample), path(fastq_adaptive_bac), path("adaptive_centrifuge_species_report.tsv"), path(fastq_non_adaptive_bac), path("non_adaptive_centrifuge_species_report.tsv")
    output:
        tuple val(sample), path("adaptive_centrifuge_bac_readID.lst"), path("non_adaptive_centrifuge_bac_readID.lst"), path("adaptive_bacterial.fastq"), path("non_adaptive_bacterial.fastq"), emit: bac_fastq_readID
	tuple val(sample), path("adaptive_centrifuge_species_report_filtered.tsv"), path("non_adaptive_centrifuge_species_report_filtered.tsv"), emit: input_krona
	tuple val(sample), path("adaptive_bacterial.fastq"), path("non_adaptive_bacterial.fastq"), emit: bac_fastq
    when:
    !skip_centrifuge_remove_contaminated
    shell:
    '''
    set +eu
    awk '$3 ~ !{params.centrifuge_reference_tax_ID}' !{"adaptive_centrifuge_species_report.tsv"} | cut -f1 | sort | uniq > adaptive_centrifuge_host_readID.lst
    awk '$3 !~ !{params.centrifuge_reference_tax_ID}' !{"adaptive_centrifuge_species_report.tsv"} | cut -f1 | sort | uniq > adaptive_centrifuge_bac_host_readID.lst
    awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' adaptive_centrifuge_host_readID.lst adaptive_centrifuge_bac_host_readID.lst > adaptive_centrifuge_bac_readID.lst
    awk '$3 ~ !{params.centrifuge_reference_tax_ID}' !{"non_adaptive_centrifuge_species_report.tsv"} | cut -f1 | sort | uniq > non_adaptive_centrifuge_host_readID.lst
    awk '$3 !~ !{params.centrifuge_reference_tax_ID}' !{"non_adaptive_centrifuge_species_report.tsv"} | cut -f1 | sort | uniq > non_adaptive_centrifuge_bac_host_readID.lst
    awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' non_adaptive_centrifuge_host_readID.lst non_adaptive_centrifuge_bac_host_readID.lst > non_adaptive_centrifuge_bac_readID.lst
    seqtk subseq !{fastq_adaptive_bac} adaptive_centrifuge_bac_readID.lst > adaptive_bacterial.fastq
    seqtk subseq !{fastq_non_adaptive_bac} non_adaptive_centrifuge_bac_readID.lst > non_adaptive_bacterial.fastq
    awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' adaptive_centrifuge_host_readID.lst adaptive_centrifuge_species_report.tsv > adaptive_centrifuge_species_report_filtered.tsv
    awk -v FS="[\t= ]" ' FNR==NR { a[$1]=$1; next } !($1 in a){print $0}' non_adaptive_centrifuge_host_readID.lst non_adaptive_centrifuge_species_report.tsv > non_adaptive_centrifuge_species_report_filtered.tsv
    cp .command.log remove_centrifuge_contaminated.log
    '''
}

process krona {
	cpus 1
	tag "${sample}"
	publishDir "$params.outdir/$sample/5_centrifuge_bac_reads",  mode: 'copy', pattern: "*krona.html", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_centrifuge_bac_reads",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(adaptive_species_report), path(non_adaptive_species_report), path(krona_database)
	output:
		tuple val(sample), path("adaptive_centrifuge_taxonomy.krona.html"), path("non_adaptive_centrifuge_taxonomy.krona.html"), emit: krona_html
		path("krona.log")
	when:
	!params.skip_krona | !params.skip_centrifuge
	script:
	"""
	set +eu
	cat ${adaptive_species_report} | cut -f 1,3 > adaptive_centrifuge_species_report.krona
	cat ${non_adaptive_species_report} | cut -f 1,3 > non_adaptive_centrifuge_species_report.krona
	ktImportTaxonomy adaptive_centrifuge_species_report.krona -o adaptive_centrifuge_taxonomy.krona.html -tax \$PWD
	ktImportTaxonomy non_adaptive_centrifuge_species_report.krona -o non_adaptive_centrifuge_taxonomy.krona.html -tax \$PWD
	cp .command.log krona.log
	"""
}

process whokaryote {
	cpus "${params.whokaryote_threads}"
	tag "${sample}"
	publishDir "$params.outdir/$sample/7_whokaryote",  mode: 'copy', pattern: "{*sv}", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/7_whokaryote",  mode: 'copy', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(fastq_adaptive_bac), path(adaptive_assembly), path(fastq_non_adaptive_bac), path(non_adaptive_assembly)
	output:
		tuple val(sample), path("adaptive_whokaryote_predictions_T.tsv"), path("non_adaptive_whokaryote_predictions_T.tsv"), emit: whokaryote_prediction
		tuple val(sample), path("adaptive_featuretable_predictions_T.tsv"), path("adaptive_featuretable.csv"), path("adaptive_contigs_genes.gff"), path("adaptive_contigs_proteins.faa"), path("adaptive_tiara_pred.txt"), path("non_adaptive_featuretable_predictions_T.tsv"), path("non_adaptive_featuretable.csv"), path("non_adaptive_contigs_genes.gff"), path("non_adaptive_contigs_proteins.faa"), path("non_adaptive_tiara_pred.txt"), emit: whokaryote_results
		path("whokaryote.log")
	when:
	!params.skip_whokaryote | !params.skip_assembly
	script:
	"""
	set +eu
	whokaryote.py --contigs ${adaptive_assembly} --outdir \$PWD --threads ${params.whokaryote_threads} 
	mv whokaryote_predictions_T.tsv adaptive_whokaryote_predictions_T.tsv
	mv featuretable_predictions_T.tsv adaptive_featuretable_predictions_T.tsv
	mv featuretable.csv adaptive_featuretable.csv
	mv contigs_genes.gff adaptive_contigs_genes.gff
	mv contigs_proteins.faa adaptive_contigs_proteins.faa
	mv tiara_pred.txt adaptive_tiara_pred.txt
	mv contigs5000.fasta adaptive_contigs5000.fasta
	whokaryote.py --contigs ${non_adaptive_assembly} --outdir \$PWD --threads ${params.whokaryote_threads}
	mv whokaryote_predictions_T.tsv non_adaptive_whokaryote_predictions_T.tsv
	mv featuretable_predictions_T.tsv non_adaptive_featuretable_predictions_T.tsv
	mv featuretable.csv non_adaptive_featuretable.csv
	mv contigs_genes.gff non_adaptive_contigs_genes.gff
	mv contigs_proteins.faa non_adaptive_contigs_proteins.faa
	mv tiara_pred.txt non_adaptive_tiara_pred.txt
 	cp .command.log whokaryote.log
	"""
}

process download_genomad_db {
    cpus 1
    label "high_memory"
    publishDir "$params.outdir/genomad_database",  mode: 'copy'
    input:
        val(db)
    output:
        path("genomad_db/*"),  emit: genomad_db
    when:
    !params.skip_download_genomad_db | !params.skip_assembly
    script:
    """
    echo ${db}
    wget ${db}
    tar -xvf genomad_db_v1.5.tar.gz
    """
}

process genomad {
    cpus 1
    tag "${sample}"
    label "high_memory"
    //publishDir "$params.outdir/$sample/8_genomad",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
    publishDir "$params.outdir/$sample/8_genomad",  mode: 'copy', saveAs: { filename -> "${sample}_$filename" }
    input:
        tuple val(sample), path(fastq_adaptive_bac), path(adaptive_assembly), path(fastq_non_adaptive_bac), path(non_adaptive_assembly), path(genomad_db)
    output:
    path("*_aggregated_classification/*_aggregated_classification.tsv")    , emit: aggregated_classification
    path("*_annotate/*_taxonomy.tsv")                                      , emit: taxonomy
    path("*_find_proviruses/*_provirus.tsv")                               , emit: provirus
    path("*_score_calibration/*_compositions.tsv")                         , emit: compositions                
    path("*_score_calibration/*_calibrated_aggregated_classification.tsv") , emit: calibrated_classification   
    path("*_summary/*_plasmid.fna")                                        , emit: plasmid_fasta
    path("*_summary/*_plasmid_genes.tsv")                                  , emit: plasmid_genes
    path("*_summary/*_plasmid_proteins.faa")                               , emit: plasmid_proteins
    path("*_summary/*_plasmid_summary.tsv")                                , emit: plasmid_summary
    path("*_summary/*_virus.fna")                                          , emit: virus_fasta
    path("*_summary/*_virus_genes.tsv")                                    , emit: virus_genes
    path("*_summary/*_virus_proteins.faa")                                 , emit: virus_proteins
    path("*_summary/*_virus_summary.tsv")                                  , emit: virus_summary
    //path "versions.yml"                                                                     , emit: versions
    when:
    !params.skip_genomad | !params.skip_assembly
    script:
    """
    genomad end-to-end --cleanup --splits 4 ${adaptive_assembly} \$PWD ${genomad_db}
    genomad end-to-end --cleanup --splits 4 ${non_adaptive_assembly} \$PWD ${genomad_db}
    """
}



workflow {
	ch_centrifuge_db=Channel.value( "${params.centrifuge_db}")
	ch_centrifuge_db.view()
	ch_genomad_db=Channel.value( "${params.genomad_db}" )
    	ch_genomad_db.view()
	Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
	.splitCsv(header:true, sep:',')
	.map { row -> tuple(row.sample_id, file(row.fastq, checkIfExists: true), file(row.csv, checkIfExists: true)) }
	.set { ch_samplesheet }
	ch_samplesheet.view()
	if (!params.skip_porechop) {
		porechop(ch_samplesheet)
		extract_adaptive_readID(porechop.out.trimmed_fastq)
	} else if (params.skip_porechop) {	
		extract_adaptive_readID(ch_samplesheet)
	}
	extract_adaptive_fastq(extract_adaptive_readID.out.extracted_readID)
	if (!params.skip_adaptive_sampling_metrics) {
		extract_adaptive_sampling_reads(ch_samplesheet)
		compute_adaptive_sampling_metrics(extract_adaptive_sampling_reads.out.extracted_fastq)
	}
	minimap(extract_adaptive_fastq.out.extracted_fastq)
	if (!params.skip_centrifuge) {
		if (!params.skip_download_centrifuge_db) {
			centrifuge_download_db(ch_centrifuge_db)
			centrifuge(minimap.out.bacterial_fastq.combine(centrifuge_download_db.out.centrifuge_db))
		} else if (params.skip_download_centrifuge_db) {		        
			ch_centrifuge_db=Channel.fromPath( "${params.outdir}/centrifuge_database/*.cf" ).collect()
			centrifuge(minimap.out.bacterial_fastq.combine(ch_centrifuge_db))			
		}
		remove_centrifuge_contaminated(centrifuge.out.bacterial_fastq)
		if (!params.skip_krona) {
			ch_krona_db=Channel.value( "${params.krona_db}")
			krona(remove_centrifuge_contaminated.out.input_krona.combine(ch_krona_db))
		}
		if (!params.skip_assembly) {
			flye(remove_centrifuge_contaminated.out.bac_fastq)
			if (!params.skip_polishing) {
				racon(flye.out.bacterial_assembly_fasta)
				medaka(racon.out.polished_racon)
				whokaryote(medaka.out.polished_medaka)
                if (!params.skip_download_genomad_db) {
                    download_genomad_db(ch_genomad_db)
                    genomad(medaka.out.polished_medaka.combine(download_genomad_db.out.genomad_db))
                } else if (params.skip_download_genomad_db) {
                    ch_genomad_db=Channel.fromPath( "${params.outdir}/genomad_database/genomad_db/" ).collect()
                    genomad(medaka.out.polished_medaka.combine(ch_genomad_db))
}
            } else if (params.skip_polishing) {
                whokaryote(flye.out.bacterial_assembly_fasta)
                if (!params.skip_download_genomad_db) {
                    download_genomad_db(ch_genomad_db)
                    genomad(flye.out.bacterial_assembly_fasta.combine(download_genomad_db.out.genomad_db))
               } else if (params.skip_download_genomad_db) {
                    ch_genomad_db=Channel.fromPath( "${params.outdir}/genomad_database/genomad_db/" ).collect()
                    genomad(flye.out.bacterial_assembly_fasta.combine(ch_genomad_db)) }
            }
        }
        }  else if (params.skip_centrifuge) {
        if (!params.skip_assembly) {
            flye(minimap.out.bacterial_fastq)
            if (!params.skip_polishing) {
                racon(flye.out.bacterial_assembly_fasta)
                medaka(racon.out.polished_racon)
                whokaryote(medaka.out.polished_medaka)
                if (!params.skip_download_genomad_db) {
                    download_genomad_db(ch_genomad_db)
                    genomad(medaka.out.polished_medaka.combine(download_genomad_db.out.genomad_db))
               } else if (params.skip_download_genomad_db) {
                    ch_genomad_db=Channel.fromPath( "${params.outdir}/genomad_database/genomad_db/" ).collect() }
                    genomad(medaka.out.polished_medaka.combine(download_genomad_db.out.genomad_db))
            } else if (params.skip_polishing) {
                whokaryote(flye.out.bacterial_assembly_fasta)
                if (!params.skip_download_genomad_db) {
                    download_genomad_db(ch_genomad_db)
                    genomad(flye.out.bacterial_assembly_fasta.combine(download_genomad_db.out.genomad_db))
               } else if (params.skip_download_genomad_db) {
                    ch_genomad_db=Channel.fromPath( "${params.outdir}/genomad_database/genomad_db/" ).collect()
		genomad(flye.out.bacterial_assembly_fasta.combine(ch_genomad_db)) }
            }
        }
    }
}
