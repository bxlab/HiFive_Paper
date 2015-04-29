#!/bin/bash

ln -sf ${BASEDIR}/Data/Timing/hicpipe.mat ${BASEDIR}/Data/Timing/hicpipe_10000.mat
ln -sf ${BASEDIR}/Data/Timing/hicpipe_frag_len_bin.f ${BASEDIR}/Data/Timing/hicpipe_10000_frag_len_bin.f
ln -sf ${BASEDIR}/Data/Timing/hicpipe_frag_gc_bin.f ${BASEDIR}/Data/Timing/hicpipe_10000_frag_gc_bin.f
ln -sf ${BASEDIR}/Data/Timing/hicpipe_map_bin.f ${BASEDIR}/Data/Timing/hicpipe_10000_map_bin.f
ln -sf ${BASEDIR}/Data/Timing/hicpipe.prior ${BASEDIR}/Data/Timing/hicpipe_10000.prior
${BASEDIR}/Scripts/HiCPipe/lscripts/bin_coords.pl ${BASEDIR}/Data/Timing/hicpipe.binned ${BASEDIR}/Data/Timing/hicpipe_10000.cbinned ${BASEDIR}/Data/Timing/hicpipe_10000.cbins 10000

${BASEDIR}/Scripts/HiCPipe/lscripts/observed_cbin_counts.pl ${BASEDIR}/Data/Timing/hicpipe_10000 0 ${BASEDIR}/Data/Timing/hicpipe_10000.o_contact

${BASEDIR}/bin/Rscript ${BASEDIR}/Scripts/HiCPipe/R/model_wrapper.r ${BASEDIR}/Data/Timning/hicpipe_10000 cbinned ${BASEDIR}/Scripts/HiCPipe/models/map_len_gc.mdl both 0 ${BASEDIR}/Data/Timing/hicpipe_10000.e_contact 0 1 1 coord_bin 1 coord_bin
sed -i 's/value/expected_count/g' ${BASEDIR}/Data/Timing/hicpipe_10000.e_contact
