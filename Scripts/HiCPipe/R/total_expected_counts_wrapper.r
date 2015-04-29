
options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <input prefix> <model file> <filter> <cis.threshold> <use cluster> <max jobs on cluster>\n",
              script.name))
  q(status=1) 
}

ifn.prefix = args[1]
model.ifn = args[2]
filter = args[3]
cis.threshold = as.numeric(args[4])
cluster = (args[5] == "1")
max.njobs = as.numeric(args[6])

mtable = read.delim(model.ifn)
mfields = mtable$field
maxvals = mtable$size

if (cluster) {
  cat("Using Sun Grid Engine cluster\n")
} else {
  cat("Not using Sun Grid Engine cluster, running sequentially on local machine\n")
}

source("R/model_predict.r")
compute.total.counts(prefix=ifn.prefix, cluster=cluster, max.njobs=max.njobs, ofields=mfields, max.vals=maxvals, filter=filter, cis.threshold=cis.threshold)

q(status=0) 
