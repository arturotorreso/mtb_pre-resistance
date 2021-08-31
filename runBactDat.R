#!/usr/bin/env Rscript

library(ape)
library(BactDating)

## Collect arguments

args = commandArgs(trailingOnly=TRUE)

if ('-r' %in% args || '--random' %in% args) {
  randomize = TRUE
  args = args[!args %in% c('-r','--random')]
} else {
  randomize = FALSE
}

if ('-f' %in% args || '--fixMu' %in% args) {
  updateMu = FALSE
  args = args[!args %in% c('-f','--fixMu')]
} else {
  updateMu = TRUE
}


## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

# Help section
if("--help"  %in% args || "-h"  %in% args || "-?" %in% args) {
  cat("
      About: Wrapper for BactDating
      Usage: ./runBactDat [OPTIONS] ml.tree
      ./runBactDat -d outdir/ ml.tree
      Options:
      -i, --input <file>          - ML tree. Dates separated by _ at the end of name.
      -r, --random                - Randomize dates.
      -d, --dir <path>            - Directory for the output. Add slash to create new directory.
                                    Default: current working directory
      -p, --prefix <str>          - Prefix for the output files. Default: input name + BactDating
      -t, --sample_times <file>   - File with dates.
                                    If absent, it will try to get them from the name (last field separated by '_').
                                    File format: no header, two columns tab separated: id date 
      
      -l, --length <int>          - Length of the alignment used for the tree.
      -m, --model <str>           - Model to use. Options: poisson, negbin, strictgamma, relaxedgamma, [mixedgamma]
      -n, --nIter <int>           - Number of MCMC iterations to perform [10e4]
      -M, --initMu <int>          - Initial rate of substitutions per genome (not per site),
                                    or zero to use root-to-tip estimate [NA]
      -f, --fixMu                 - If present, don't update the mutation rate [FALSE]
      -u, --updateRoot            - TRUE/FALSE to update the root (default: FALSE for rooted trees, TRUE for unrooted trees)
      -b, --minbralen             - Minimum branch length for the input tree (in number of substitutions)
                                    Default: Minimum branch length of tree
      -S  --initSigma             - Initial std on per-branch substitution rate
      -A  --initAlpha             - Initial coalescent time unit
      -s, --showProgress          - TRUE/FALSE to show progress bar (default: FALSE)
      -T, --thinInt               - Thinning interval (default: ceiling(nIter/1000))
      -h, -?, --help              - This help message.
      \n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg value or -a value)

args = paste(args[c(TRUE, FALSE)],args[c(FALSE, TRUE)])

parseArgs <- function(x) strsplit(sub("-*", "", x), " ")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))

argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#Set the working directory or create a new one

wd = ifelse ("d" %in% names(argsL) || "dir" %in% names(argsL),
             argsL[names(argsL) %in% c('d','dir')][[1]],
             getwd())

if (strsplit(wd, "")[[1]][length(strsplit(wd,"")[[1]])] == '/'){
  if (!dir.exists(file.path(wd))) {dir.create(wd)}
} else {
  wd = paste(wd,'/',sep="")
}

input_tree = argsL[names(argsL) %in% c('i','input')][[1]]
alignment_length = as.numeric(argsL[names(argsL) %in% c('l','length')][[1]])
prefix = ifelse ("p" %in% names(argsL) || "prefix" %in% names(argsL),
                 argsL[names(argsL) %in% c('p','prefix')][[1]],
                 paste(input_tree,'_BactDating', sep=''))

dates_file = ifelse ("t" %in% names(argsL) || "sample_times" %in% names(argsL),
                     argsL[names(argsL) %in% c('t','sample_times')][[1]],
                     FALSE)


model = ifelse ("m" %in% names(argsL) || "model" %in% names(argsL),
                argsL[names(argsL) %in% c('m','model')][[1]],
                'mixedgamma')

nIter = ifelse ("n" %in% names(argsL) || "nIter" %in% names(argsL),
                argsL[names(argsL) %in% c('n','nIter')][[1]],
                10e4)
nIter = as.integer(nIter)

initMu = ifelse ("M" %in% names(argsL) || "initMu" %in% names(argsL),
                as.numeric(argsL[names(argsL) %in% c('M','initMu')][[1]]),
                NA)

updateRoot = ifelse ("u" %in% names(argsL) || "updateRoot" %in% names(argsL),
                 argsL[names(argsL) %in% c('u','updateRoot')][[1]],
                 NA)

minbralen = ifelse ("b" %in% names(argsL) || "minbralen" %in% names(argsL),
                    argsL[names(argsL) %in% c('b','minbralen')][[1]],
                    0.1)

initAlpha = ifelse ("A" %in% names(argsL) || "initAlpha" %in% names(argsL),
                    argsL[names(argsL) %in% c('A','initAlpha')][[1]],
                    NA)

initSigma = ifelse ("S" %in% names(argsL) || "initSigma" %in% names(argsL),
                    argsL[names(argsL) %in% c('S','initSigma')][[1]],
                    NA)

showProgress = ifelse ("s" %in% names(argsL) || "showProgress" %in% names(argsL),
                     argsL[names(argsL) %in% c('s','showProgress')][[1]],
                     FALSE)

thinInt = ifelse ("T" %in% names(argsL) || "thinInt" %in% names(argsL),
                argsL[names(argsL) %in% c('T','thinInt')][[1]],
                ceiling(nIter/1000))

thinInt = as.numeric(thinInt)
minbralen = as.numeric(minbralen)
initMu = as.numeric(initMu)
if (!is.na(initAlpha)){initAlpha = as.numeric(initAlpha)}
if (!is.na(initSigma)){initSigma = as.numeric(initSigma)}

tracefile = paste(wd, prefix, '_traces.txt', sep='')

#####################
####### CODE #######
####################

tree = read.tree(input_tree)

# Check whether to update root or not
if (is.na(updateRoot)) {
  updateRoot = !is.rooted(tree)
}

# Check if length is in snps/site or total number of snps.
if (sum(tree$edge.length) < 2) {
  tree$edge.length=tree$edge.length*alignment_length
}


if (!isFALSE(dates_file)){
  dates = read.table(dates_file, header = F, stringsAsFactors = F)
  dates = setNames(dates$V2, dates$V1)
} else {
  dates = sapply(tree$tip.label, function (x){as.integer(strsplit(x, '_')[[1]][length(strsplit(x, '_')[[1]])])})
}


if (randomize) {
  print ('This is gonna be random')
  dates_names = names(dates)
  dates_random = sample(unname(dates), replace=F)
  names(dates_random) = dates_names
  dates = dates_random
}

if (isFALSE(minbralen)){
  minbralen = min(tree$edge.length)
}

# Run BactDating

print (paste(
  'Running BactDating with ',
  nIter,
  ' iterations, with the ',
  model,
  ' model. Starting mutation rate: ',
  initMu,
  '. Minimum branch length: ',
  minbralen, 
  '. Update the mutation rate: ',
  updateMu,
  '. Starting Sigma: ',
  initSigma, 
  '. Starting Alpha: ',
  initAlpha,
  '. Update root: ',
  updateRoot,
  '. Recording samples every ',
  thinInt,
  ' iterations.',
  sep=''
))



res=bactdate(tree = tree,
             date = dates,
             nbIts = nIter,
             model = model,
             initMu = initMu,
             initSigma = initSigma,
             initAlpha = initAlpha,
             updateMu = updateMu,
             updateRoot = updateRoot,
             minbralen = minbralen,
             showProgress = showProgress)


# OUTPUT
save(res, file = paste(wd, prefix, '_res.RData', sep=''))
print (res)
