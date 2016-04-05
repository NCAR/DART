#! /opt/R/bin/Rscript --vanilla

## Avoid shell arithimetic:
## you can pass arguments as text calculations you'd like R to perform.

args<-commandArgs(TRUE)

if(length(args)!=4) {
  print('usage: ./randomSeq distName length param1 param2')
  q()
}

doEvalParseText <- function(x) eval(parse(text=x))

if(length(args)==4) {
  distName <- args[1]
  length <- as.numeric(doEvalParseText(args[2]))
  param1 <- as.numeric(doEvalParseText(args[3]))
  param2 <- as.numeric(doEvalParseText(args[4]))
}

distFunction <- 
  switch( tolower(distName),
          'uniform'  = runif,
          'normal'   = rnorm,
          'gaussian' = rnorm,
          { error('No such distribution yet provided for.')
            q()
          }
         )

cat(do.call(distFunction, list(length, param1, param2)), fill=1)
