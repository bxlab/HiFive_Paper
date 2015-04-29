#########################################################################################################
# Internal functions
#########################################################################################################

get.short.fn=function(fn)
{  
  grep = gregexpr("/", fn)[[1]]
  substr(fn, grep[length(grep)]+1, nchar(fn))
}

get.ofield.args=function(ranges, fields, suffix)
{
  from = ranges[paste(fields, suffix, "from", sep="_")]
  to = ranges[paste(fields, suffix, "to", sep="_")]
  paste(length(fields), paste(fields, from, to, collapse=" "))
}

model.predict=function(ofield.ranges, lib.dir="",
                       fends.fn="h1_cl_gc_1000.binned",
                       tmp.dir="tmp/", log.dir="log/",
                       prior=1, mfields=NULL, mfields.maxvals=NULL, mfields.fns=NULL,
                       filter="trans", cis.threshold=100000,
                       ofields.x=c("f1", "f2"), ofields.y=c("f1", "f2"))
{
  if (!is.null(mfields))
    model.args = paste(length(mfields), paste(mfields, mfields.maxvals, mfields.fns, collapse=" "))
  else
    model.args = 0
  ofield.x.args = get.ofield.args(ofield.ranges, ofields.x, "x")
  ofield.y.args = get.ofield.args(ofield.ranges, ofields.y, "y")
  
  index = ofield.ranges["index"]
  from.fend = ofield.ranges["from.fend"]
  to.fend = ofield.ranges["to.fend"]
  
  binary = paste(lib.dir, "/bin/model_integrate", sep="")
  
  ofn=paste(tmp.dir, get.short.fn(fends.fn), ".", index, sep="")
  log=paste(log.dir, get.short.fn(fends.fn), ".", index, sep="")
  command = paste(binary, fends.fn, ofn, from.fend, to.fend, filter, cis.threshold, prior,
                  model.args, ofield.x.args, ofield.y.args, ">>", log, "2>&1", sep=" ")
  cat(sprintf("%s\n", command))
  
  system(paste("echo hostname: `hostname` > ", log, sep=""))
  system(paste("echo command: ", command, " >> ", log, sep=""))
  
  # remove ofn before run
  system(paste("rm -rf", ofn))

  sleep = round(runif(1, 1, 10))
  system(paste("echo sleeping ",sleep," before command ... >> ", log, sep=""))
  system(paste("sleep", sleep))
  system(paste("echo sleep done, running command >> ", log, sep=""))
  
  ec = system(command)
  if (ec != 0)
  {  
    cat(sprintf("command: %s\n", command))
    stop(paste("Command failed, error code:", ec))
  }

  return (ofn)
}

expand.ranges=function(fields=c("f1", "f2") , from.vals=c(1,1), to.vals=c(10,10), num.splits=4)
{
  dim = length(fields)
  num.splits.dim = num.splits^(1/dim)
  field.nparts = rep(1, length(fields))
  for (i in 1:dim)
    field.nparts[i] = min(to.vals[i] - from.vals[i] + 1, ceiling(num.splits.dim))
  
  params = list()
  for (i in 1:dim)
  {
    levels = seq(from.vals[i], to.vals[i])
    if (field.nparts[i] > 1) {
      lookup.table = data.frame(t(sapply(split(levels ,cut(levels, field.nparts[i])), range)))
      rownames(lookup.table) = NULL
      colnames(lookup.table) = c("from", "to")
    } else {
      lookup.table = data.frame(from=from.vals[i], to=to.vals[i])
    }
    params[[fields[i]]] = lookup.table
  }
  index.list = rev(lapply(params, function(x) 1:(dim(x)[1])))
  egrid = expand.grid(index.list, stringsAsFactors=F)

  result = NULL
  for (i in 1:dim)
  {
    ind = egrid[,fields[i]]
    lookup.table = params[[fields[i]]]
    t = lookup.table[ind,]
    colnames(t) = paste(fields[i], c("from", "to"), sep="_")
    rownames(t) = NULL
    result = if(i==1) t else cbind(result, t)
  }

  result
}

model.predict.split=function(lib.dir="", split.by.output.bins=F,
                             fends.fn="h1_cl_gc_1000.binned",
                             ofn="h1_cl_gc_1000.expected_counts",
                             cluster=T, filter="trans", cis.threshold=100000,
                             max.njobs=400, prior=1, fends=NULL,
                             mfields=NULL, mfields.maxvals=NULL, mfields.fns=NULL,
                             ofields.x=c("frag_len_bin","frag_gc_bin"), from.vals.x=c(1,1), to.vals.x=c(5,5),
                             ofields.y=c("frag_len_bin","frag_gc_bin"), from.vals.y=c(1,1), to.vals.y=c(5,5))
{
  cat(sprintf("fends file: %s\n", fends.fn))
  cat(sprintf("ofn: %s\n", ofn))
  cat(sprintf("prior: %f\n", prior))
  cat(sprintf("filter: %s, cis threshold: %d\n", filter, cis.threshold))
  if (!is.null(mfields))
    cat(sprintf("model fields: %s (%s)\n", paste(mfields, collapse=","), paste(mfields.maxvals, collapse=",")))
  else
    cat(sprintf("no model fields\n", paste(mfields, collapse=",")))
  cat(sprintf("output fields X: %s (%s)\n", paste(ofields.x, collapse=","),
              paste(from.vals.x, to.vals.x, sep="-", collapse=",")))
  cat(sprintf("output fields Y: %s (%s)\n", paste(ofields.y, collapse=","),
              paste(from.vals.y, to.vals.y, sep="-", collapse=",")))
  cat(sprintf("Contact filter: %s, cis threshold: %d\n", filter, cis.threshold))
  
  ofields = c(paste(ofields.x, "x", sep="_"), paste(ofields.y, "y", sep="_"))
  from.vals = c(from.vals.x, from.vals.y)
  to.vals = c(to.vals.x, to.vals.y)
  nbins = prod(to.vals - from.vals + 1)

  num.splits = max.njobs
  
  # split run according to output space or fends space
  if (is.null(fends) || split.by.output.bins) {
    cat(sprintf("Splitting jobs according to output bin space\n"))
    # explode ranges into table
    ranges = expand.ranges(fields=ofields, from.vals=from.vals, to.vals=to.vals, num.splits=num.splits)
    ranges$from.fend = rep(0, dim(ranges)[1])
    ranges$to.fend = rep(0, dim(ranges)[1])
    concat = T
  } else {
    cat(sprintf("Splitting jobs according to fend space\n"))
    fend.breaks = sapply(split(fends,cut(fends, num.splits)), range)
    sranges = expand.ranges(fields=ofields, from.vals=from.vals, to.vals=to.vals, num.splits=1)
    ranges = NULL
    for (i in 1:dim(fend.breaks)[2])
      ranges = rbind(ranges, cbind(sranges, fend.breaks[1,i], fend.breaks[2,i]))
    names(ranges)[dim(ranges)[2]-1] = "from.fend"
    names(ranges)[dim(ranges)[2]] = "to.fend"
    concat = F
  }

  # add running index field
  cat('DEBUGGING\n')
  cat(sprintf("Fends: %i, Split: %i\n",as.integer(is.null(fends)),as.integer(split.by.output.bins)))
  cat(sprintf("Range names: %s\n",paste(names(ranges),collapse=' ')))
  ranges$index = 1:dim(ranges)[1]
  #cat(sprintf("%s\n",paste(paste(paste(ranges$map_bin_x_from,ranges$map_bin_x_to,sep='-'),paste(ranges$frag_gc_x_from,ranges$frag_gc_x_to,sep='-'),paste(ranges$frag_len_x_from,ranges$frag_len_x_to,sep='-'),sep=' '),collapse='\n')))
  cat(sprintf("%s\n",paste(paste(paste(ranges$from.fend,ranges$to.fend,sep='-'),collapse=','))))
  
  ntasks = dim(ranges)[1]
  njobs = min(ntasks, max.njobs)
  cat(sprintf("Number of bins: %f\n", nbins))
  cat(sprintf("Number of tasks: %d\n", ntasks))
  cat(sprintf("Average square bins per tasks: %f\n", round(nbins/ntasks)))
  cat(sprintf("Number of running jobs: %d\n", njobs))

  tmp.dir = paste(lib.dir, "/tmp/", sep="")
  log.dir = paste(lib.dir, "/log/", sep="")
  
  system(sprintf("mkdir -p %s", tmp.dir))
  system(sprintf("mkdir -p %s", log.dir))
  
  system(paste("rm -rf ", log.dir, get.short.fn(fends.fn), "*", sep=""))
  system(paste("rm -rf ", tmp.dir, get.short.fn(fends.fn), "*", sep=""))
  system(paste("rm -rf ", ofn, sep=""))

  library(Rsge)
  sge.setDefaultOptions()
  sge.options(sge.save.global=F)
  sge.prefix = paste("tmp/Rsge.", get.short.fn(fends.fn), ".", sep="")
  cat(sprintf("Rsge temp files: %s*\n", sge.prefix))
  sge.options(sge.file.prefix=sge.prefix)

  cat(sprintf("Prefix: %s\n",sge.getOption('sge.file.prefix')))
  cat(sprintf("Trace: %s\n",as.character(sge.getOption('sge.trace'))))
  cat(sprintf("Debug: %s\n",as.character(sge.getOption('sge.debug'))))
  cat(sprintf("Qsub: %s\n",as.character(sge.getOption('sge.qsub'))))
  cat(sprintf("Qsub Options: %s\n",as.character(sge.getOption('sge.qsub.options'))))
  cat(sprintf("User options: %s\n",as.character(sge.getOption('sge.user.options'))))
  cat(sprintf("Qsub blocking: %s\n",as.character(sge.getOption('sge.qsub.blocking'))))
  cat(sprintf("Script: %s\n",as.character(sge.getOption('sge.script'))))
  
  result = sge.parRapply(ranges, model.predict, lib.dir=lib.dir,
                         njobs=njobs, join.method=c, cluster=cluster,
                         fends.fn=fends.fn, log.dir=log.dir,
                         prior=prior, mfields=mfields, mfields.maxvals=mfields.maxvals, mfields.fns=mfields.fns,
                         filter=filter, cis.threshold=cis.threshold, ofields.x=ofields.x, ofields.y=ofields.y, 
                         function.savelist=c("get.short.fn", "get.ofield.args"))
  if (class(result) == "list" && class(result[[1]]) == "try-error")
    return (-1)

  cat(sprintf("Merging results into %s\n", ofn))
  if (concat) {
    for (i in seq_along(result)) {
      if (!file.exists(result[i]))
        stop(paste("file does not exist:", result[i]))

      if (i==1)
        ec = system(sprintf("cat %s > %s", result[i], ofn))
      else
        ec = system(sprintf("cat %s | %s/lscripts/remove_header.pl >> %s", result[i], lib.dir, ofn))
      if (ec != 0)
        stop("Error merging results")
    }
  } else {
    ofields = c(paste(ofields.x, "1", sep=""), paste(ofields.y, "2", sep=""))
    table = NULL
    cat(sprintf("merging %d files\n", length(result)))
    for (i in seq_along(result)) {
      if (!file.exists(result[i]))
        stop(paste("file does not exist:", result[i]))

      ttable = read.delim(result[i], stringsAsFactors=F)
      if (i == 1) {
        table = ttable
        next
      }
      if (any(ttable[,ofields] != table[,ofields]))
        stop("result does not have the same field order")
      table$value = table$value + ttable$value
    }
    write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")
  }
  
  system(paste("rm -rf ", tmp.dir, get.short.fn(fends.fn), "*", sep=""))
}

omit.dups=function(table, nfields)
{
  diff = table[,1:nfields+nfields] - table[,1:nfields]

  if (nfields == 1)
    return (table[diff>=0,])
  
  ind = rep(F, dim(table)[1])
  eq = rep(T, dim(table)[1])
  for (i in 1:nfields)
  {
    if (i<nfields)
      ind = ind | (eq & (diff[,i] > 0))
    else
      ind = ind | (eq & (diff[,i] >= 0 ))
    eq = eq & (diff[,i] == 0)
  }
  table[ind,]
}

#########################################################################################################
# Wrapper functions
#########################################################################################################

# used to setup the model - get total possible counts for each bin
compute.total.counts=function(prefix="h1_cl_gc_200", cluster=F,
                              ofields=c("frag_gc_bin", "map_bin"),
                              max.njobs=400,
                              max.vals=c(20,5), filter="trans", cis.threshold=100000)
{
  fends.fn = paste(prefix, ".binned", sep="")
  table = read.delim(fends.fn, stringsAsFactors=F)
  fends = table$fend
  
  ofn = paste(prefix, ".total_counts", sep="")

  if (cluster) {
    has.cluster = (system("which qsub") == 0)
    if (!has.cluster)
      stop("must run on host which supports sun grid engine (qsub)")
  }
  
  model.predict.split(lib.dir=getwd(), split.by.output.bins=T, max.njobs=max.njobs, 
                      fends.fn=fends.fn, ofn=ofn, cluster=cluster, fends=fends, 
                      filter=filter, cis.threshold=cis.threshold,
                      ofields.x=ofields, from.vals.x=rep(1, length(ofields)), to.vals.x=max.vals,
                      ofields.y=ofields, from.vals.y=rep(1, length(ofields)), to.vals.y=max.vals)
  
  table = read.delim(ofn, stringsAsFactors=F)
  result = omit.dups(table, length(ofields))
  names(result)[dim(result)[2]] = "total"
  write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

# used to check model performance in a user-defined xy plane
compute.expected.counts=function(ifn.prefix="h1_cl_gc_1000_1M", cluster=T,
                                 fend.suffix="cbinned", max.njobs=400,
                                 filter="trans", cis.threshold=100000,
                                 model.ifn="map_len_gc.model",
                                 ofn="h1_cl_gc_1000_1M.e_contact",
                                 ofields.x=c("coord_bin"), ofields.y=c("coord_bin"))
{
  if (cluster) {
    has.cluster = (system("which qsub") == 0)
    if (!has.cluster)
      stop("must run on host which supports sun grid engine (qsub)")
  }
  
  fends.ifn = paste(ifn.prefix, fend.suffix, sep=".")
  prior.ifn = paste(ifn.prefix, "prior", sep=".")

  # get all bin ranges 
  table = read.delim(fends.ifn, stringsAsFactors=F)
  fends = table$fend
  for (field in c(ofields.x, ofields.y))
    if (!is.element(field, names(table)))
      stop(paste("field", field, "not found"))
  from.vals.x = numeric(length(ofields.x))
  to.vals.x = numeric(length(ofields.x))
  from.vals.y = numeric(length(ofields.y))
  to.vals.y = numeric(length(ofields.y))
  for (i in seq_along(ofields.x))
  {  
    from.vals.x[i] = range(table[,ofields.x[i]])[1]
    to.vals.x[i] = range(table[,ofields.x[i]])[2]
  }
  for (i in seq_along(ofields.y))
  {  
    from.vals.y[i] = range(table[,ofields.y[i]])[1]
    to.vals.y[i] = range(table[,ofields.y[i]])[2]
  }

  # prior
  prior = read.delim(prior.ifn, header=F)
  prior = prior[1,1]

  # model parameters
  mtable = read.delim(model.ifn)
  if (dim(mtable)[1] > 0)
  {  
    mfields = mtable$field
    mfields.maxvals = mtable$size
    mfields.fns = paste(paste(ifn.prefix, mtable$field, sep="_"), ".f", sep="")
  }
  else {
    mfields = NULL
    mfields.maxvals = NULL
    mfields.fns = NULL
  }
  model.predict.split(lib.dir=getwd(), split.by.output.bins=T, fends.fn=fends.ifn, ofn=ofn, cluster=cluster, fends=fends, 
                      filter=filter, cis.threshold=cis.threshold, prior=prior, max.njobs=max.njobs, 
                      mfields=mfields, mfields.maxvals=mfields.maxvals, mfields.fns=mfields.fns,
                      ofields.x=ofields.x, from.vals.x=from.vals.x, to.vals.x=to.vals.x,
                      ofields.y=ofields.y, from.vals.y=from.vals.y, to.vals.y=to.vals.y)
}
