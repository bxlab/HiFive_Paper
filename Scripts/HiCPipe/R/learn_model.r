
options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <ifn prefix> <dataset> <model file> [only seed]\n", script.name))
  q(status=1) 
}

ifn.prefix = args[1]
dataset = args[2]
model.fn = args[3]
only.seed = ifelse (length(args) >= 4, argv[4], F)
cat(sprintf("learn.model(ifn.prefix=\"%s\", dataset=\"%s\", model.fn=\"%s\", only.seed=%s)\n",
            ifn.prefix, dataset, model.fn, only.seed))

source("R/model.r")
learn.model(ifn.prefix=ifn.prefix, dataset=dataset, model.fn=model.fn, clean=T, only.seed=only.seed)
q(status=0) 
