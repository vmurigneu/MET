# Adaptive sequencing for metagenomic sequencing of low biomass samples

The project aimed to develop a field-friendly adaptive sequencing pipeline to increase sensitivity of metagenomic sequencing of low biomass samples.
A bioinformatic pipeline was developed to automate the analysis of the metagenomic samples sequenced on the ONT MinION platform. 

Read more about the project developed at the University of Queensland Genome Innovation Hub: https://gih.uq.edu.au/project/adaptive-sequencing-metagenomic-sequencing-low-biomass-samples

## Overall pipeline 

### 1. Porechop

The reads are trimmed for adapters using [Porechop](https://github.com/rrwick/Porechop). 

### 2. Extraction of adaptive and non-adaptive reads

The trimmed reads are splitted into two fastq files: the adaptive reads (Adaptive sampling channel number 1-256) and the non-adaptive reads (Regular sequencing channel : 257-512). The adaptive reads are selected based on the adaptive sampling report (classified as "stop_receiving reads").   

### 3. Adaptive sampling metrics 

[Nanocomp](https://github.com/wdecoster/nanocomp) is used to compute metrics (Median Read Length, Read N50, Median Read Quality) for the different categories of adaptive sampling reads: stop_receiving, unblock and no_decision. 

### 4. Minimap2 mapping

The adaptive and non-adaptive reads are mapped to the reference genome provided by the user using [Minimap2](https://github.com/lh3/minimap2). Host reads identified by Minimap2 are excluded from the fastq files. 

### 5. Centrifuge taxonomy classification

Host removed reads from the previous step are used as input to the classifier [Centrifuge](https://ccb.jhu.edu/software/centrifuge/). Host reads identified by Centrifuge are subsequently excluded from the fastq files. The taxonomy ID of the reference genome is a required parameter of the pipeline: "centrifuge_reference_tax_ID" (e.g. "9606" for Homo sapiens or "9913" for Bos taurus). [Krona](https://github.com/marbl/Krona/wiki) is then used to visualise the taxonomy results as pie charts. Users can search for their host taxonomy ID on the [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) -- please note that since the centrifuge database is not updated regularly, certain taxa listed on NCBI will not have a match in centrifuge and will thus produce an error.    

### 6. Flye assembly and polishing

The host removed reads (adaptive and non-adaptive) are assembled using the software [Flye](https://github.com/fenderglass/Flye) (metagenome mode). The draft assemblies are subsequently polished using [Racon](https://github.com/isovic/racon) and [Medaka](https://github.com/nanoporetech/medaka). The model parameter selected to run Medaka (e.g. r1041_e82_400bps_sup_g615) should correspond to the model used for the basecalling (e.g. dna_r10.4.1_e8.2_400bps_sup.cfg).  

### 7. Eukaryote and prokaryote classification

The assembly contigs (adaptive and non-adaptive) are classified into prokaryote or eukaryote contigs using the software [whokaryote](https://github.com/LottePronk/whokaryote). Note that contigs with less than 2 genes or shorter than the minimum size are not classified and do not appear in the output.

### 8. 	Virus and plasmid classification

Viruses and plasmids are predicted using the software [genomad](https://github.com/apcamargo/genomad). If either of the assemblies are empty, genomad will be skipped on the respective assembly. If both adaptive and non-adaptive assemblies are empty, than genomad will error out and the pipeline should be rerun with --skip_genomad and -resume options. 

### 9. 	Aviary Recover MAGs
   
This step will recover MAGs from provided assembly (adaptive and non-adaptive) using a variety of binning algorithms, as implemented in the module recover from the software [aviary](https://github.com/rhysnewell/aviary/).

### Additional notes on donwloading databases 
As the size of the centrifuge and genomad databases is large and the download takes up considerable time, if repeated runs of the pipeline are required, it is recommended that users move the databases from the results folder of one analysis into the results folder of the subsequent analysis, providing the options to skip the download of the databases on the subsequent run. Please be aware that this will only work if the different analyses are not run simultaneously, in which case a database for each analysis is required to be downloaded.

## Usage

### 0. Required input files

**a) Basecalled read files (fastq)**

The basecalling step is not included in the pipeline. Raw ONT fast5 files needs to be converted into basecalled reads using Guppy. 

**b) Adaptive sampling report file (csv)**

An example of adaptive sampling file is provided [here]()

**c) Samplesheet file (csv)**

Prepare a samplesheet containing one line for each sample with the following information: the sampleID (sample_id), the path to the corresponding basecalled read file (fastq) and the path to the corresponding adaptive sampling report (csv). File paths are given in relation to the workflow base directory, they are not absolute paths.        
```
sample_id,fastq,csv
cowmilk2,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk2_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV67209_fa368ac1_d487a5df.csv
cowmilk12,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk12_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV81134_516a2fa0_5cce3b1f.csv
cowmilk30,fastq/Milk_samples_MET/2023_05_22_MET_cowmilk30_guppy6.3.8_sup.fastq.gz,fastq/Milk_samples_MET/adaptive_sampling_FAV64853_8f20760e_1cbbd159.csv
```

**d) Reference genome sequence file (fasta)**

For instance for the human genome, the file can be downloaded from: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz. The file can be provided as a gzip.   

**e) Nextflow configuration file (nextflow.config)**

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html).  

The pipeline uses separated Singularity containers for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the nextflow.config file ([see documentation](https://www.nextflow.io/docs/latest/singularity.html)). Ensure that you have sufficient space in your assigned singularity directory as images can be large.   

An example configuration file can be found [here](https://github.com/vmurigneu/DIS/blob/main/nextflow.config). 

**f) Nextflow main script (main.nf)**

The main.nf script contains the pipeline code and is generally not user-modifiable. 

## Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

Adapter trimming:
* `--porechop_args`: Porechop optional parameters (default=""), see [details](https://github.com/rrwick/Porechop#full-usage)
* `--porechop_threads`: number of threads for Porechop (default=4)
* `--skip_porechop`: skip the Porechop trimming step (default=false)

Adaptive Sampling Read Sequencing: :
* `--skip_adaptive_sampling_metrics`: skip the Adaptive sampling metrics step (default=false)
* `--nanocomp_threads`: Number of threads for NanoComp (default=4)

Mapping: 
* `--minimap_threads`: Number of threads for Minimap2 (default=12)

Taxonomy classification:
* `--skip_download_centrifuge_db`: skip the centrifuge database downloading step (default=false)
* `--skip_centrifuge`: skip the centrifuge taxonomy classification step (default=false)
* `--centrifuge_db`: path to download the centrifuge database (default='https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz')
* `--centrifuge_reference_tax_ID`: Taxonomy ID for reference genome (default='9606' for homo sapiens)
* `--skip_centrifuge_remove_contaminated`: skip the emoval of centrifuge reads from contaminated reference step (default=false)
* `--centrifuge_threads`: number of threads for Centrifuge (default=12)
* `--skip_krona`: skip the generation of Krona plots (default=false)
 
Assembly:
* `--flye_args`: Flye optional parameters (default=`--flye_args "--meta"`), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* `--flye_threads`: number of threads for Flye (default=4)
* `--memory`: Memory usage for Flye (default=0)
* `--skip_assembly`: skip the metagenome assembly step (default=false)

Polishing:
* `--racon_nb`: number of Racon long-read polishing iterations (default=4)
* `--racon_args`: Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500"), see [details](https://github.com/isovic/racon#usage)
* `--racon_threads`: number of threads for Racon (default=4)
* `--medaka_threads`: number of threads for Medaka (default=8)
* `--medaka_model`: name of the Medaka model (default=r1041_e82_400bps_sup_g615, see [details](https://github.com/nanoporetech/medaka#models)

Eukaryote and prokaryote classification:
* `--skip_whokaryote`:	Skip whokaryote classification (default=false)
* `--whokaryote_threads`:	Number of threads for Whokaryote (default=8)

Virus and Plasmid classification:
* `--skip_download_genomad_db`:	Skip the genomad database download if it is already present locally (default=false)
* `--skip_genomad`:	Skip genomad classification (default=false)

Aviary Recover MAGs:
* `--skip_aviary`: Skip aviary recover (default=false)
* `--aviary_threads`: Number of threads for Aviary (default=8)
* `--pplacer_threads`: Number of threads for Aviary (default=8)
* `--max_memory_aviary`:	Maximum memory for Aviary (default=500)
* `--checkm_db`:	Path to the CheckM2 database
* `--gtdb_path`:	Path to the GTDB database
* `--eggnog_db`:	Path to the eggnog-mapper database

## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
The main output folder (`--outdir`) will contain a folder per sample (the folder is named as in the column sample_id in the samplesheet file)

Each sample folder will contain the following folders:
* **1_trimming:** Fastq files containing trimmed reads (sample_id_trimmed.fastq.gz) 
* **2_adaptive:** Fastq files containing the adaptive reads (sample_id_adaptive.fastq) and non-adaptive reads (sample_id_non_adaptive.fastq) in separate files
* **3_minimap:** Fastq files containing the host removed reads after the Minimap2 mapping step for the adaptive (sample_id_adaptive_bac.fastq) and non-adaptive reads (sample_id_non_adaptive_bac.fastq)  
* **4_centrifuge:** Centrifuge taxonomy classification results for the host removed reads for both adaptive and non-adaptive reads, see [details](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output)  
  * Centrifuge classification output: classification assignments for a read (*_centrifuge_species_report.tsv)  
  * Centrifuge summary output: classification summary for each genome or taxonomic unit (*_centrifuge_report.tsv)   
* **5_centrifuge_bac_reads:**  
  * Fastq files containing the host removed reads after both Minimap2 and Centrifuge steps (sample_id_adaptive_bacterial.fastq and sample_id_non_adaptive_bacterial.fastq)
  * List of host removed read identifiers (sample_id_adaptive_centrifuge_bac_readID.lst and sample_id_non_adaptive_centrifuge_bac_readID.lst)
  * Krona pie chart HTML report for the host removed reads (after both Minimap2 and Centrifuge steps)  
* **6_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output). The final polished asssembly fasta files are sample_id_adaptive_flye_polished.fasta and sample_id_non_adaptive_flye_polished.fasta.  
* **7_whokaryote:** Whokaryote output files (.csv, .txt, .tsv, .fasta, .gff, .faa), see [details](https://github.com/LottePronk/whokaryote) for adaptive and non-adaptive assemblies.
* **8_genomad:** Genomad output files (), see [details]()
* **9_aviary:** Aviary recover output files contained in the folders benchmarks/, bins/, data/, diversity/, taxonomy/ for adaptive and non-adaptive bins. 
