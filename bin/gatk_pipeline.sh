#!/bin/sh

#############################
# les directives PBS vont ici:

# Your job name (displayed by the queue)
#PBS -N GATK_3.4-0

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
module load GATK/3.4-0 

for file in $(ls /home/jrudewicz/MICADo/data/fastq/pool0/*.fastq)
do

	BNfile=$(basename $file .fastq)

	## Alignement (done with STAR_v.2.4.0g1)	
	samtools view -bS /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"Aligned.out.sam > /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile".bam
	samtools sort /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile".bam /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"
	samtools index /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile".bam

	## Split'N'Trim
	GATK  -T SplitNCigarReads -R /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -I /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile".bam -o /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"Nsplited.bam -U ALLOW_N_CIGAR_READS --allow_potentially_misencoded_quality_scores

	# ## RealignerTargetCreator 
	GATK -T RealignerTargetCreator -R /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -I /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"Nsplited.bam -o /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"forIndelRealigner.intervals --allow_potentially_misencoded_quality_scores

	# ## IndelRealigner
	GATK -T IndelRealigner -R /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -I /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"Nsplited.bam -targetIntervals /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"forIndelRealigner.intervals -o /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"realignedBam.bam --allow_potentially_misencoded_quality_scores

	## HaplotypeCaller
	GATK -T HaplotypeCaller -R /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta -I /home/jrudewicz/GATK/PhasePaper/AlignmentFiles_GRCh37_v2.4.0g1/"$BNfile"realignedBam.bam -o /home/jrudewicz/GATK/PhasePaper/vcf/"$BNfile".raw.snps.indels.vcf --intervals 17 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --allow_potentially_misencoded_quality_scores -dontUseSoftClippedBases --maxReadsInRegionPerSample 10000 --min_base_quality_score 30 --forceActive

	ToGGVCFs=$ToGGVCFs" --variant /home/jrudewicz/GATK/PhasePaper/vcf/"$BNfile".raw.snps.indels.vcf"

done

## GenotypeGVCFs
GATK -T GenotypeGVCFs -R /home/jrudewicz/GATK/PhasePaper/Genome/GRCh37/GRCh37.fasta $ToGGVCFs --intervals 17 -o /home/jrudewicz/GATK/PhasePaper/vcf/PilotSamplesPB.vcf

# all done
echo "Job finished" 
echo "Date:" `date`
