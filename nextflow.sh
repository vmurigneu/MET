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
samplesheet=${dir}/samplesheet/test.csv
out_dir=${dir}/results_test
fqdir=${dir}/fastq

#cp /QRISdata/Q4860/MinION_Full-Scale/Raw-Data/Mock_community/20230118_comm_1Pro_BB/20230118_1539_MN33062_FAS65356_85211c55/other_reports/adaptive_sampling_FAS65356_85211c55_7bbc80e1.csv ${dir}/fastq
#cp /QRISdata/Q4860/MinION_Full-Scale/Raw-Data/Mock_community/20230120_comm_2Un_BB/20230120_1227_MN41645_FAS67931_f441b463/other_reports/adaptive_sampling_FAS67931_f441b463_c4415279.csv ${dir}/fastq
nextflow main.nf --outdir ${out_dir} --fastqdir ${fqdir} --samplesheet ${samplesheet}


