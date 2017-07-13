#!/bin/bash

BAM=$1
BED=$2

qsub -cwd -b y -N sid$BAM -l h_vmem=8g -e ${BAM}.genotype.log -o ${BAM}.genotype.log "module load samtools; samtools view -F 4 $BAM -L $BED | /u/rdenroche/git/spb-analysis-tools/phoenixPipe/sampleIdentity/robCaller_v2.0.pl -q 0 -d 0 -f 0 | /u/rdenroche/git/spb-analysis-tools/phoenixPipe/sampleIdentity/robGenotyper.pl $BED > ${BAM}.genotype" > ${BAM}.genotype.log


