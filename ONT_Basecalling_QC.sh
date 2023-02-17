#!/bin/sh
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=2:mem=4gb
#PBS -N NF_ONT_Basecalling
#PBS -j oe

Project_Dir=/rds/general/ephemeral/user/dthorbur/ephemeral/00_Work/05_Nanopore/03_Nanopore_QC
cd $Project_Dir

module load nextflow/22.04.4
module unload -f java/jdk-12
module load java/jdk-16

echo "Starting:`date`"

nextflow run ONT_Basecalling_QC.nf \
	-c ONT_Basecalling_QC.config \
	--profile imperial \
	--Input "/rds/general/ephemeral/user/dthorbur/ephemeral/00_Work/05_Nanopore/01_Tarballs" \
	--Container "/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/ONT_Guppy_GPU.sif" && echo "Finished:`date`"
