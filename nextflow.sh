#!/bin/bash
#
#SBATCH --time=20:00:00
#SBATCH --job-name=MET_pipeline
#SBATCH --output=./s%j_job.pipeline_test.out
#SBATCH --error=./s%j_job.pipeline_test.error
#SBATCH --account=a_gih
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4

module load nextflow/22.10.1  

#directory containing the nextflow.config file and the main.nf script
dir=/scratch/project/gihcomp/MET/pipeline
cd ${dir}
datadir=${dir}/data
out_dir=${dir}/results_test
fqdir=${dir}/fastq

#copy raw data to scratch
#mkdir ${bamdir}
#no control for that dataset, cannot use the pipeline as it is
#recall_medici /QRISdata/Q1207/E_coli/20210701_Sequel64123_0032/raw_data/3_C01_long/m64123_210703_141742.subreads.bam*
#cp /QRISdata/Q1207/E_coli/20210701_Sequel64123_0032/raw_data/3_C01_long/m64123_210703_141742.subreads.bam* ${bamdir}

nextflow main.nf --outdir ${out_dir} --fastqdir ${fqdir}


