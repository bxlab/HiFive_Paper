SHELL:=/bin/bash

###########################################################################
# user defined parameters are in the cfg file
###########################################################################
BASEDIR=$(CURDIR)
CFG=$(BASEDIR)/Settings/dependencies.txt
EMPTY=
include $(CFG)
export


###########################################################################
all:check folders genome mapping hiclib hifive hicpipe hicnorm matrixbalancing
###########################################################################

###########################################################################
# check that all required software is present
###########################################################################
check:
ifeq ($(shell which python2.7),$(EMPTY))
	$(error No python found.)
endif
ifeq ($(shell which awk),$(EMPTY))
	$(error No awk found.)
endif
ifeq ($(shell which wget),$(EMPTY))
	$(error No wget found.)
endif
ifndef BOWTIE
	$(error No BOWTIE defined.)
endif
ifeq ($(shell which ${BOWTIE}),$(EMPTY))
	$(error No bowtie found.)
endif
ifndef BOWTIE_BUILD
	$(error No BOWTIE_BUILD defined.)
endif
ifeq ($(shell which ${BOWTIE_BUILD}),$(EMPTY))
	$(error No bowtie-build found.)
endif
ifndef BOWTIE_INDEX_DIR
	$(error No BOWTIE_INDEX_DIR defined.)
endif
ifndef SAMTOOLS
	$(error No SAMTOOLS defined.)
endif
ifeq ($(shell which ${SAMTOOLS}),$(EMPTY))
	$(error No samtools found.)
endif
ifndef FASTQ_DUMP
	$(error No FASTQ_DUMP defined.)
endif
ifeq ($(shell which ${FASTQ_DUMP}),$(EMPTY))
	$(error No fastq-dump found.)
endif

###########################################################################
# create needed folder structure
###########################################################################
folders:
	@mkdir -p $(BASEDIR)/Data/Genome
	@mkdir -p $(BASEDIR)/Data/HiC
	@mkdir -p $(BASEDIR)/Scripts/bin
	@mkdir -p $(BASEDIR)/Scripts/downloads
	@mkdir -p $(BASEDIR)/Data/FiveC/Fastq
	@mkdir -p $(BASEDIR)/Data/FiveC/Mapping
	@mkdir -p $(BASEDIR)/Data/HiC/Fastq
	@mkdir -p $(BASEDIR)/Data/HiC/Mapping
	@mkdir -p $(BASEDIR)/Data/HiC/Data
	@mkdir -p $(BASEDIR)/bin
	@ln -s $(abspath $(BOWTIE)) $(BASEDIR)/bin/bowtie 
	@ln -s $(abspath $(BOWTIE_BUILD)) $(BASEDIR)/bin/bowtie_build 
	@ln -s $(abspath $(SAMTOOLS)) $(BASEDIR)/bin/samtools
	@ln -s $(abspath $(FASTQ_DUMP)) $(BASEDIR)/bin/fastq-dump
	@ln -s $(abspath $(BOWTIE_INDEX_DIR)) $(BASEDIR)/Data/Genome/bowtie_indices
	@ln -s $(abspath $(RSCRIPT)) $(BASEDIR)/bin/Rscript
	@ln -s $(abspath $(TMPDIR)) $(BASEDIR)/tmp
	@if [ ! -d "$(BASEDIR)/Data/Genome/MM9_Fasta" ]; ln -s $(abspath $(MM9_FA_PATH)) $(BASEDIR)/Data/Genome/mm9_fasta; fi
	@if [ ! -d "$(BASEDIR)/Data/Genome/HG19_Fasta" ]; ln -s $(abspath $(HG19_FA_PATH)) $(BASEDIR)/Data/Genome/hg19_fasta; fi
	@export PYTHONPATH=$(BASEDIR)/lib/python2.7/site-packages:$(PYTHONPATH)


###########################################################################
# download and uncompress genome partition data
###########################################################################
genome:
	cd $(BASEDIR)/tmp && wget http://files.figshare.com/2043822/hg19_data.tar.bz2
	cd $(BASEDIR)/Data/Genome && tar xjf $(BASEDIR)/tmp/hg19_data.tar.bz2
	cd $(BASEDIR)/tmp && wget http://files.figshare.com/2043821/mm9_data.tar.bz2
	cd $(BASEDIR)/Data/Genome && tar xjf $(BASEDIR)/tmp/mm9_data.tar.bz2

###########################################################################
# download data and map reads to primers/genome
###########################################################################
mapping:
	@ cd Makefiles && $(MAKE) -f Makefile.mapping

###########################################################################
# Install HiCLib libraries
###########################################################################
hiclib:
	@cd Makefiles && $(MAKE) -f Makefile.hiclib

###########################################################################
# Perform analysis with HiFive
###########################################################################
hifive:
	@ cd Makefiles && $(MAKE) -f Makefile.hifive

###########################################################################
# Perform analysis with HiFive
###########################################################################
hicpipe:
	@ cd Makefiles && $(MAKE) -f Makefile.hicpipe

###########################################################################
# Perform analysis with HiCNorm
###########################################################################
hicnorm:
	@ cd Makefiles && $(MAKE) -f Makefile.hicnorm

###########################################################################
# Perform analysis with Matrix Balancing
###########################################################################
matrixbalancing:
	@ cd Makefiles && $(MAKE) -f Makefile.mb

.SECONDARY:
.PHONY: all