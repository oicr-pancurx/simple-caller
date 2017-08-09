#!/bin/bash

BAM=$1
BED=$2

qsub -cwd -b y -N sid$BAM -l h_vmem=8g -e ${BAM}.genotype.log -o ${BAM}.genotype.log "module load samtools; samtools view -F 4 $BAM -L $BED | /.mounts/labs/PCSI/production/simple-caller/robCaller_v2.0.pl -q 0 -d 0 -f 0 | /.mounts/labs/PCSI/production/simple-caller/robGenotyper.pl $BED > ${BAM}.genotype" > ${BAM}.genotype.log


