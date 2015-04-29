#!/bin/bash

${BASEDIR}/bin/hifive hic-normalize express -m 500000 -x 0 -o ${BASEDIR}/Data/Timing/hifive_exp.hcp -e 1000 -g 0.000005 -d -w cis -k -f 10 ${BASEDIR}/Data/Timing/hifive.hcp