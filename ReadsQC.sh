#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=2:mem=4gb
#PBS -N NF_ReadsQC
#PBS -j oe

Project_Dir=/path/to/project/directory
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting:`date`"

## To extract reads from Anopheles gambiae after QC completes, uncomment the last 2 arguments and run again. 
nextflow run ReadsQC.nf \
	-c ReadsQC.config \
	--profile imperial \
	--ReadType "IlluminaPE" \
	--InDir "/path/to/fastqs" \
	--KrakenDB "/path/to/Kraken2DB" ## --KrakenExtract --KNE_TaxID "7165"
