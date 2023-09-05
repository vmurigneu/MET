# Adaptive sequencing for metagenomic sequencing of low biomass samples

The project aimed to develop a field-friendly adaptive sequencing pipeline to increase sensitivity of metagenomic sequencing of low biomass samples.
A bioinformatic pipeline was developed to automate the analysis of the metagenomic samples sequenced on the ONT MinION platform. 

Read more about the project developed at the University of Queensland Genome Innovation Hub: https://gih.uq.edu.au/project/adaptive-sequencing-metagenomic-sequencing-low-biomass-samples

## Overall pipeline 

### 1. Porechop

The 

### 2. 

The

## Usage

### 0. Required input files

**a) Basecalled read files (fastq)**


**b) Adaptive sampling report file (csv)**

An example of adaptive sampling file is provided [here]()

**c) Samplesheet file (csv)**

Prepare a samplesheet containing one line for each sample with the following information: the sampleID (sample_id), the path to the corresponding basecalled read file (fastq) and the path to the corresponding adaptive sampling report (csv).     
```
sample_id,fastq,csv
cowmilk2,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk2_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV67209_fa368ac1_d487a5df.csv
cowmilk12,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk12_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV81134_516a2fa0_5cce3b1f.csv
cowmilk30,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk30_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV64853_8f20760e_1cbbd159.csv
```

**d) Reference genome sequence file (fasta)**

For instance for the human genome, the file can be downloaded from: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.

**e) Nextflow configuration file (nextflow.config)**

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html). 

The pipeline uses separated Singularity containers for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the nextflow.config file ([see documentation](https://www.nextflow.io/docs/latest/singularity.html)). 

An example configuration file can be found [here](https://github.com/vmurigneu/DIS/blob/main/nextflow.config). 

**f) Nextflow main script (main.nf)**

The main.nf script contains the pipeline code and is generally not user-modifiable. 

## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
The main output folder (`--outdir`) will contain a folder per sample (the folder is named as in the column sample_id in the samplesheet file)

Each sample folder will contain the following folders:
* **1_trimming:** Fastq files containing filtered reads (sample_id_trimmed.fastq.gz) 
* **2_adaptive:** Fastq files containing the adaptive (sample_id_adaptive.fastq) and non-adaptive reads (sample_id_non_adaptive.fastq)
* **3_minimap:** 
* **4_centrifuge:** 
* **5_centrifuge_bac_reads:** QUAST quality assessment report, see [details](http://quast.sourceforge.net/docs/manual.html)
* **6_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output). The final polished asssembly fasta files are sample_id_adaptive_flye_polished.fasta and sample_id_non_adaptive_flye_polished.fasta.  

