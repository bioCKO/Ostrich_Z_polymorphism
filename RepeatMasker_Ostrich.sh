#!/bin/bash -l

#SBATCH -A b2016026
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 07:00:00
#SBATCH -J RepeatMasker
#SBATCH --mail-user homa.papoli@ebc.uu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools
module load RepeatMasker

RepeatMasker -species aves -pa 16 -gccalc -gff -dir /proj/b2016026/private/Ostrich_Z_LinkageMap_project/Ostrich_gigadb/Ostrich_repeatMask Struthio_camelus.20130116.OM.fa
