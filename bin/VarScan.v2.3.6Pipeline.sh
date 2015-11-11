#!/bin/sh

#############################
# les directives PBS vont ici:

# Your job name (displayed by the queue)
#PBS -N VarScan.v2.3.6

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

module load samtools/1.2

coorN = "17:7,578,208-7,579,527"
coorC = "17:7,573,983-7,578,258"

for file in $(ls /home/jrudewicz/MICADo/data/fastq/pool0/*.fastq)
do

	BNfile=$(basename $file .fastq)

	## Alignement (done with STAR_v.2.4.0g1)	
	samtools view -bS /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"Aligned.out.sam > /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".bam
	samtools sort /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".bam /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile"
	samtools index /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".bam

	if [[ $BNfile =~ ^'N' ]]
	then 
		samtools mpileup -B -f /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -r "$coorN" -Q 3 -d 20000 /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".bam > /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".txt							
	else
		samtools mpileup -B -f /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -r "$coorC" -Q 3 -d 20000 /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".bam > /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".txt							
	fi
	java -jar /module/apps/varscan/2.3.6/VarScan.v2.3.6.jar mpileup2cns /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".txt > /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/VarScan/"$BNfile".vcf --strand-filter 0 --min-coverage 5 --min-reads2 5 --min-avg-qual 60 --min-var-freq 0.05 --output-vcf 1 --variants 1

done
# all done
echo "Job finished" 
echo "Date:" `date`
