
options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <input prefix> <fend suffix> <model file> <filter> <cis.threshold> <ofn> <use cluster> <max jobs on cluster> <ofields.x.count> <ofields.x> <ofields.y.count> <ofields.y>\n",
              script.name))
  cat(sprintf("      Note that <model file> can be set to 'none'\n"))
  q(status=1) 
}

ifn.prefix = args[1]
fend.suffix = args[2]
model.ifn = args[3]
filter = args[4]
cis.threshold = as.numeric(args[5])
ofn = args[6]
cluster = (args[7] == "1")
max.njobs = as.numeric(args[8])

# get x fields
ofields.x.count = as.integer(args[9])
index = 10+ofields.x.count
ofields.x = args[10:(index-1)]
# get y fields
ofields.y.count = as.integer(args[index])
ofields.y = args[(index+1):length(args)]
if (ofields.y.count != length(ofields.y))
  stop("error in ofields.y.count")

cat(sprintf("ofields.x: %s\n", paste(ofields.x, collapse=",")))
cat(sprintf("ofields.y: %s\n", paste(ofields.y, collapse=",")))

if (cluster) {
  cat("Using Sun Grid Engine cluster\n")
} else {
  cat("Not using Sun Grid Engine cluster, running sequentially on local machine\n")
}

source("R/model_predict.r")


cat(sprintf("compute.expected.counts(ifn.prefix=\"%s\", fend.suffix=\"%s\", model.ifn=\"%s\", filter=\"%s\", cis.threshold=%d, ofn=\"%s\", ofields.x=c(\"%s\"), ofields.y=c(\"%s\"))\n",
            ifn.prefix, fend.suffix, model.ifn, filter, cis.threshold, ofn, paste(ofields.x, collapse=","), paste(ofields.y, collapse=",")))

compute.expected.counts(ifn.prefix=ifn.prefix, cluster=cluster, max.njobs=max.njobs, fend.suffix=fend.suffix, model.ifn=model.ifn,
                        filter=filter, cis.threshold=cis.threshold, ofn=ofn, ofields.x=ofields.x, ofields.y=ofields.y)

q(status=0) 
