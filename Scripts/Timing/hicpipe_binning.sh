#!/bin/bash

${BASEDIR}/bin/Rscript ${BASEDIR}/Scripts/HiCPipe/R/model_preprocess.r ${BASEDIR}/Data/Timing/hicpipe fends ${BASEDIR}/Scripts/HiCPipe/models/map_len_gc.mdl 3 ${BASEDIR}/Data/Timing/hicpipe.nm far_cis 500000 0 1 ${BASEDIR}/bin/Rscript