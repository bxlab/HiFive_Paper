#!/bin/bash

${BASEDIR}/bin/hifive hic-normalize probability -m 500000 -x 0 -o ${BASEDIR}/Data/Timing/hifive_prob.hcp -b 1000 -g 0.0005 -p -l 0.5 ${BASEDIR}/Data/Timing/hifive.hcp
