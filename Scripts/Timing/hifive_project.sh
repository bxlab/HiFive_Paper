#!/bin/bash

${BASEDIR}/bin/hifive hic-project -f 10 -m 500000 -x 0 -j 1000 -n 50 ${BASEDIR}/Data/Timing/hifive/hifive.hcd ${BASEDIR}/Data/Timing/hifive.hcp
