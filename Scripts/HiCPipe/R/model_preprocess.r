
options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <prefix> <fends suffix> <model file> <bin round> <output file> <filter> <cis threshold> <cluster> <max jobs on cluster> <Rscript>\n", script.name))
  q(status=1) 
}

source("R/bin_fields.r")
library(R.utils)
npath=function(fn)
{
  split = strsplit(fn, split = "/")[[1]]
  if (length(split) > 1)
  {    
    sfn = split[length(split)]
    return (paste(normalizePath(getParent(fn)), sfn, sep="/"))
  } else {
    return (fn)
  }
  
}

prefix = npath(args[1])
fend.suffix = args[2]
model.fn = npath(args[3])
round = as.numeric(args[4])
ofn = npath(args[5])
filter = args[6]
cis.threshold = as.numeric(args[7])
cluster = args[8]
max.njobs = args[9]
rscript.binary = args[10]

fends.fn = paste(prefix, fend.suffix, sep=".")
mtable = read.delim(model.fn, stringsAsFactors=F)
fields = mtable[, "field"]
rfields = mtable[mtable$type == "optimize", "raw_field"]
qns = mtable[mtable$type == "optimize", "size"]

# 1) bin
cat("Step I: binning fields\n")
cat(sprintf("binning fields: %s\n", paste(rfields, collapse=",")))
cat(sprintf("command: bin.fields(fends.fn=\"%s\", ofn.prefix=\"%s\", fields=c(%s), qns=c(%s), round=%f)\n",
            fends.fn, prefix, paste(rfields, collapse=","), paste(qns, collapse=","), round))
bin.fields(fends.fn=fends.fn, ofn.prefix=prefix, fields=rfields, qns=qns, round=round)

# 2) observed counts
cat("Step II: observed counts\n")
command = paste("./lscripts/compute_n.pl", prefix, "none", filter, cis.threshold, "F", paste(fields, collapse=" "))
cat(sprintf("%s\n", command))
if (system(command) != 0)
  stop("error in compute_n.pl")

# 3) expected counts
cat("Step III: expected counts\n")
command = paste("srun ", rscript.binary, " R/total_expected_counts_wrapper.r", prefix, model.fn, filter, cis.threshold, cluster, max.njobs)
cat(sprintf("%s\n", command))
if (system(command) != 0)
  stop("error in total_expected_counts_wrapper.r")

# 4) unite results
cat("Step IV: unite results\n")
bfields = sapply(fields, function(x) paste(x, 1:2, sep=""))
command = paste("./lscripts/append_table.pl ", prefix, ".total_counts ",
                prefix, ".counts ", ofn, " count 0 F ", paste(bfields, collapse=" "), sep="")
cat(sprintf("%s\n", command))
if (system(command) != 0)
  stop("error in append_table.pl")

cat("model_preprocess.r done\n")

q(status=0) 
