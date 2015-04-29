
#####################################################################################################
# Utility gmake canned recipe
#####################################################################################################

# checks return code of last shell run
define ASSERT
if [ $$? -ne 0 ]; then exit 1; fi
endef

# user defined parameters are in the cfg file
include $(CFG)

######################################################################################################
# Parameters: experimental bias modelling
######################################################################################################

# file describing normaliztion model
MFN?=models/map_len_gc.mdl

# estimate biases using (trans/far_cis/close_cis/both) contacts
FILTER?=far_cis

# threshold between close and far cis definitions, relevant if FILTER equals far_cis or close_cis
CIS_THRESHOLD?=500000

######################################################################################################
# Parameters: output contact maps
######################################################################################################

# maximal simultaneous jobs running on sge cluster
CLUSTER_MAXJOBS?=24

######################################################################################################
# Parameters: R paths
######################################################################################################

# path to the Rscript binary
RSCRIPT?=Rscript

######################################################################################################
# Main pipeline: 
# - Check input files
# - Build regular and smooth contact matrices
######################################################################################################

all: check_user_params init
# print parameters
	@echo Input fragment ends table: $(FENDS_TABLE)
	@echo Input mat file: $(MAT_FILE)
	@echo Output directory: $(ODIR)
	@echo Dataset title: $(DATASET)
	@echo Star activity threshold \(bp\): $(SEGMENT_LEN_THRESHOLD)
	@echo Correction model file: $(MFN)
	@echo Correction bias type: $(FILTER)
	@echo close/far cis threshold \(relevant only for cis biases\): $(CIS_THRESHOLD)
	@if [ "${CLUSTER}"=="1" ]; then echo Distributing over Sun Grid Engine; else echo Running sequentially on local machine; fi
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@ $(MAKE) $(ODIR)/$(DATASET).model

check_user_params:
# check user parameters
	@if [ -z ${FENDS_TABLE} ]; then echo FENDS_TABLE not defined, see README; exit 1; fi
	@if [ -z ${MAT_FILE} ]; then echo MAT_FILE not defined, see README; exit 1; fi
	@if [ -z ${DATASET} ]; then echo DATASET not defined, see README; exit 1; fi
	@if [ -z ${ODIR} ]; then echo ODIR not defined, see README; exit 1; fi
# check input files
	@if [ ! -e ${FENDS_TABLE} ]; then echo Input fends table ${FENDS_TABLE} does not exist; exit 1; fi
	@if [ ! -e ${MAT_FILE} ]; then echo Input fends table ${MAT_FILE} does not exist; exit 1; fi

# remove output directory
clean: check_user_params
	rm -rf $(ODIR)

#####################################################################################################
# Binaries
#####################################################################################################

bin/smooth_matrix: src/smooth_matrix.cpp
	@mkdir -p bin
	g++ $^ -O2 -o $@ -Wall
bin/model_integrate: src/model_integrate.cpp
	@mkdir -p bin
	g++ $^ -O2 -o $@ -Wall
init: bin/smooth_matrix bin/model_integrate

.DELETE_ON_ERROR:
.SECONDARY:

######################################################################################################
# PIPELINE
######################################################################################################

######################################################################################################
# Step 0: dataset initializtion
#
# create link from map_bin.f and map.bin_ranges into results dir
# extract limited fends table from full fends table
#
######################################################################################################

$(ODIR)/%.fends $(ODIR)/%.mat: $(MAT_FILE)
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 0: initialization, $*

	mkdir -p $(ODIR)
	mkdir -p tmp
	mkdir -p log
	ln -sf $(CURDIR)/const_correction/map_bin.f $(ODIR)/$*_map_bin.f
	ln -sf $(CURDIR)/const_correction/map.bin_ranges $(ODIR)/$*_map.bin_ranges
	ln -sf $(MAT_FILE) $(ODIR)/$*.mat
	head -n 1 ${FENDS_TABLE} > $(ODIR)/$*.fends
	cat ${FENDS_TABLE} | awk '$$5 == 1' >> $(ODIR)/$*.fends

######################################################################################################
# Step 2: prepare nm file
# create the nm file from the mat and fend tables
######################################################################################################

$(ODIR)/%.nm: $(ODIR)/%.fends $(ODIR)/%.mat
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 2: compute nm table, $*
	srun $(RSCRIPT) R/model_preprocess.r $(ODIR)/$* fends $(MFN) 3 $@ $(FILTER) $(CIS_THRESHOLD) $(CLUSTER) $(CLUSTER_MAXJOBS) $(RSCRIPT)

######################################################################################################
# Step 3: learn model
# Create model biases (.f files) from the nm file
######################################################################################################

$(ODIR)/%.model: $(ODIR)/%.nm
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 3: learn model, $*
	srun $(RSCRIPT) R/learn_model.r $(ODIR)/$* $* $(MFN)
	touch $@
