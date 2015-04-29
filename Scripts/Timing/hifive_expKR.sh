#!/bin/bash

${BASEDIR}/bin/hifive hic-normalize express -m 500000 -x 0 -o ${BASEDIR}/Data/Timing/hifive_expKR.hcp -e 1000 -g 0.000000000001 -w cis -k -z -f 10 ${BASEDIR}/Data/Timing/hifive_nodist.hcp