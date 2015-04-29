#!/bin/bash

for J in {6..8}; do ${BASEDIR}/Scripts/Support/bam2raw.py ${BASEDIR}/Data/Timing/SRR44388${J}_sub_1.bam ${BASEDIR}/Data/Timing/SRR44388${J}_sub_2.bam ${BASEDIR}/Data/Timing/SRR44388${J}_sub.raw; done
