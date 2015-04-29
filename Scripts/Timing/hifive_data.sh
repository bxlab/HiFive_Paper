#!/bin/bash

${BASEDIR}/bin/hifive fends -B ${BASEDIR}/Data/Timing/mm9_NcoI_sub.bed ${BASEDIR}/Data/Timing/hifive.fends

${BASEDIR}/bin/hifive hic-data --bam ${BASEDIR}/Data/Timing/SRR443886_sub_1.bam ${BASEDIR}/Data/Timing/SRR443886_sub_2.bam --bam ${BASEDIR}/Data/Timing/SRR443887_sub_1.bam ${BASEDIR}/Data/Timing/SRR443887_sub_2.bam --bam ${BASEDIR}/Data/Timing/SRR443888_sub_1.bam ${BASEDIR}/Data/Timing/SRR443888_sub_2.bam -i 500 ${BASEDIR}/Data/Timing/hifive.fends ${BASEDIR}/Data/Timing/hifive.hcd
