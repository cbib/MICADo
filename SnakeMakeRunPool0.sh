#!/bin/sh

#############################
# les directives PBS vont ici:

# Your job name (displayed by the queue)
#PBS -N PB-pool0ll

# Specify the working directory
#PBS -d .

# walltime (hh:mm::ss)
#PBS -l walltime=300:00:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(ppn=) to be used
#PBS -l nodes=1:ppn=15

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
#echo "PBS_NODEFILE: " `cat $PBS_NODEFILE | uniq`
echo "#############################" 

#############################

# What you actually want to launch

if hash module 2>/dev/null; then
        module load python/2.7.9
        module load python/3.4.3
       `snakemake --bash-completion`
       echo "Loaded python\n"
fi

snakemake --snakefile ./SnakefilePool0 --cores 15 some_samples

# all done
echo "Job finished"