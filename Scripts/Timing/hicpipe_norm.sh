#!/bin/bash

${BASEDIR}/bin/Rscript ${BASEDIR}/Scripts/HiCPipe/R/learn_model.r ${BASEDIR}/Data/Timing/hicpipe hicpipe ${BASEDIR}/Scripts/HiCPipe/models/map_len_gc.mdl
touch ${BASEDIR}/Data/Timing/hicpipe.model
