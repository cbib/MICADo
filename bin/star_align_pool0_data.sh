#!/bin/sh

#############################
# les directives PBS vont ici:

# Your job name (displayed by the queue)
#PBS -N STAR_2.4.0g1_GRCh37


# Put stout and stderr to homedir
#PBS -k oe

# Specify the working directory
#PBS -d .

# walltime (hh:mm::ss)
#PBS -l walltime=300:00:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(ppn=) to be used
#PBS -l nodes=1:ppn=2

## Mail directive
#PBS -m abe
#PBS -M justinerudewicz@hotmail.fr

# fin des directives PBS
#############################

# modules cleaning
module purge

# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "PBS_JOBID:" $PBS_JOBID
echo "PBS_O_WORKDIR:" $PBS_O_WORKDIR
echo "PBS_NODEFILE: " `cat $PBS_NODEFILE | uniq`
echo "#############################" 

#############################

module load STAR/2.4.0g1

for file in $(ls /home/jrudewicz/MICADo/data/fastq/pool0/*.fastq)
do

	BNfile=$(basename $file .fastq)

	## Alignement + Sam to Bam + Sort + Index	
	STAR --runThreadN 16 --genomeDir /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37 --sjdbGTFfile /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/Homo_sapiens.GRCh37.75.gtf --readFilesIn /home/jrudewicz/MICADo/data/fastq/pool0/"$BNfile".fastq --outSAMattributes All --outFileNamePrefix /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile" --outSAMattrRGline ID:EORTC10994 SM:"$BNfile" PL:PACBIO --outSAMmapqUnique 40

done

# all done
echo "Job finished" 
echo "Date:" `date`
