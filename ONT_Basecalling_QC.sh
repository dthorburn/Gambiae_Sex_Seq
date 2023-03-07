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

## Defining the help message. To use run 'sh ONT_Basecalling_QC.sh --help'
case $1 in
 -[h?] | --help)
	nextflow run ONT_Basecalling_QC.nf -c ONT_Basecalling_QC.config --help
	exit 0;;
esac

echo "Starting:`date`"
nextflow run ONT_Basecalling_QC.nf \
	-c ONT_Basecalling_QC.config \
	--profile imperial \
	--Input "/rds/general/ephemeral/user/dthorbur/ephemeral/00_Work/05_Nanopore/01_Tarballs" \
	--Container "/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/ONT_Guppy_GPU.sif" && echo "Finished:`date`"
