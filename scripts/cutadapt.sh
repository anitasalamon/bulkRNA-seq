#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=300gb
#SBATCH --time=7-00:00:00
#SBATCH -p standard
#SBATCH -A owens_rivanna

module load cutadapt/3.4
module load gparallel/20170822

cat sampleinfo.txt | parallel -j 4 "cutadapt -q 10 -u 3 -U 3 -g ^AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG --max-n 0.1 --pair-filter=both --discard-casava -o /project/Owens_Rivanna/00.Raw.and.Aligned.Files/bulk-RNA-seq/2022.09.AS_VS_M001/usftp21.novogene.com/02_filtr/{}_1.trimmed.fq.gz -p /project/Owens_Rivanna/00.Raw.and.Aligned.Files/bulk-RNA-seq/2022.09.AS_VS_M001/usftp21.novogene.com/02_filtr/{}_2.trimmed.fq.gz {}_1.fq.gz {}_2.fq.gz"
