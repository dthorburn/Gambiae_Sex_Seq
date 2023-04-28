#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N NF_RNA_Isoform
#PBS -j oe


module load nextflow/20.10.0
## Added a help message that is easier to get to than before
## To print help message use 'sh RNA_Isoform.sh --help'
case $1 in
 -[h?] | --help)
	nextflow run RNA_Isoform.nf -c RNA_Isoform.config --help
	exit 0;;
esac

Project_Dir=/rds/general/project/tmstorage/live/matteo_sexseq/RNA_Seq/14_Parallel_Testing
cd $Project_Dir

## If the run timed out, resume the pipeline by removing the '#' before -resume on line 31. Note the cache can be a bit funny sometimes and 
## some processes may restart despite being completed.
## If you are using uLTRA be aware there are Tb of temporary files that only delete if uLTRA is successful. Use '-w /path/to/ephemeral' can
## work around this issue. Oddly using '-w ${EPHEMERAL}' throws permission errors.  
echo "Starting:`date`"
nextflow run RNA_Isoform.nf \
	-c RNA_Isoform.config \
	-w /rds/general/ephemeral/user/dthorbur/ephemeral/work2 \
	--Mapper minimap2 \
	--profile imperial \
	--MAPDir /rds/general/user/dthorbur/home/tmstorage/live/matteo_sexseq/RNA_Seq/14_Parallel_Testing/02_Mapped \
	--RefGen /rds/general/project/tmstorage/live/matteo_sexseq/RNA_Seq/12_Testing_Tools/00_References/VectorBase-61_AgambiaePEST_Genome.fasta \
	--Annots /rds/general/project/tmstorage/live/matteo_sexseq/RNA_Seq/12_Testing_Tools/00_References/VectorBase-61_AgambiaePEST.gff \
	--ESPPath /rds/general/user/dthorbur/home/00_Private_Tools/espresso/src --Email d.thorbur@ic.ac.uk # -resume 

#	--InDir /rds/general/ephemeral/user/dthorbur/ephemeral/decompressed_fastqs \