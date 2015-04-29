#!/bin/bash

${BASEDIR}/bin/hifive hic-normalize binning -m 500000 -x 0 -o ${BASEDIR}/Data/Timing/hifive_bin.hcp -r 1000 -t 1.0 -y cis -v gc,len,mappability -s 20,20,5 -u even,even,fixed-const ${BASEDIR}/Data/Timing/hifive_nodist.hcp