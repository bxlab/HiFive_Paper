#!/bin/bash

mkdir -p ${BASEDIR}/Data/Timing/tmp
mkdir -p ${BASEDIR}/Data/Timing/log
ln -sf ${BASEDIR}/Scripts/HiCPipe/const_correction/map_bin.f ${BASEDIR}/Data/Timing/hicpipe_map_bin.f
ln -sf ${BASEDIR}/Scripts/HiCPipe/const_correction/map.bin_ranges ${BASEDIR}/Data/Timing/hicpipe_map.bin_ranges
head -n 1 ${BASEDIR}/Data/Timning/mm9_NcoI_sub.fend > ${BASEDIR}/Data/Timing/hicpipe.fends
cat ${BASEDIR}/Data/Timing/mm9_NcoI_sub.fend | awk '$5 == 1' >> ${BASEDIR}/Data/Timing/hicpipe.fends

${BASEDIR}/Scripts/HiCPipe/lscripts/coords2fends.pl ${BASEDIR}/Data/Timing/mm9_NcoI_sub.fend ${BASEDIR}/Data/Timing/hicpipe 0 500 ${BASEDIR}/Data/Timing/hicpipe.random ${BASEDIR}/Data/Timing/SRR443886_sub.raw ${BASEDIR}/Data/Timing/SRR443887_sub.raw ${BASEDIR}/Data/Timing/SRR443888_sub.raw
