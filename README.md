# Adaptive sequencing for metagenomic sequencing of low biomass samples


Read more about the project developed at the University of Queensland Genome Innovation Hub: https://gih.uq.edu.au/project/adaptive-sequencing-metagenomic-sequencing-low-biomass-samples

## Overall pipeline 

### 1. C

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
