library(rwrfhydro)
gregorianOrigin <- "1601-01-01 00:00:00 UTC"

## compare a saved restart file against
##   1) one in the current ensemble advance.
##   2) what in the posterior? file

#setwd("")

ensMember <- '0005'

cat(
"Verifying RESTARTS for the single ensemble member: ", ensMember, "

Checking restart file information from 4 locations:
1) current ensemble restarts in main dir (end)
2) OUTPUT/model_integration...
3) Prior_Diag.nc
4) Posterior_Diag.nc

")

rst.forDate <- GetNcdfFile(paste0('restart.hydro.',ensMember,'.nc'), q=TRUE)
dateFormat='%Y-%m-%d_%H:%M:%S'
dateFormatNoSec='%Y-%m-%d_%H:%M'
dateFormatCompactNoSec='%Y%m%d%H%M'
rst.date <- as.POSIXct(attr(rst.forDate,'global')$Restart_Time,
                       format=dateFormat, tz='UTC')

priorFile <- 'Prior_Diag.nc'
prior <- GetNcdfFile(priorFile, q=TRUE)
prior.date <- as.POSIXct(tail(prior$time,1) * (24*3600),
                         origin=gregorianOrigin, tz='UTC')

posteriorFile <- 'Posterior_Diag.nc'
posterior <- GetNcdfFile(posteriorFile, q=TRUE)
posterior.date <- as.POSIXct(tail(posterior$time,1) * (24*3600),
                         origin=gregorianOrigin, tz='UTC')

if(rst.date != prior.date) stop('rst.date != prior.date')
if(posterior.date != prior.date) stop('posterior.date != prior.date')

## HYDRO_RST files
fileNameOutput <-
    paste0("OUTPUT/model_integration.",
           format(rst.date-3600, format=dateFormatCompactNoSec),
           '-',
           format(rst.date, format=dateFormatCompactNoSec),
           '.',ensMember,'/',
           'HYDRO_RST.',
           format(rst.date, format=dateFormatNoSec),
           '_DOMAIN1.', ensMember, '.nc')
fileNameCurrent <-
    paste0('restart.hydro.',ensMember,'.nc')
rst.output <- GetNcdfFile(fileNameOutput, q=TRUE, var='qlink1')
rst.current <- GetNcdfFile(fileNameCurrent, var='qlink1', q=TRUE)
cat('HYDRO_RST: qlink1 output & current match:',
    identical(rst.output$qlink1, rst.current$qlink1),'\n')

rst.prior <- prior$qlink1[,as.integer(ensMember)+2, length(prior$time)]
cat('HYDRO_RST: qlink1 output & prior match:',
    identical(as.numeric(rst.output$qlink1) , rst.prior),'\n')

## This should NOT match. The OUTPUT files are the priors
rst.posterior <- posterior$qlink1[,as.integer(ensMember)+2, length(posterior$time)]
cat('HYDRO_RST: qlink1 output & posterior do not match:',
    !identical(as.numeric(rst.output$qlink1) , rst.posterior),'\n')

## assim only files
fileNameOutput <-
    paste0("OUTPUT/model_integration.",
           format(rst.date-3600, format=dateFormatCompactNoSec),
           '-',
           format(rst.date, format=dateFormatCompactNoSec),
           '.',ensMember,'/restart.assimOnly.nc')
fileNameCurrent <-
    paste0('restart.assimOnly.',ensMember,'.nc')
rst.output <- GetNcdfFile(fileNameOutput, q=TRUE)
rst.current <- GetNcdfFile(fileNameCurrent, q=TRUE)
cat('assimOnly output and current match:', identical(rst.output, rst.current),'\n')

## SERIOUSLY MESSED UP, the index order of the copy dimension is reversed
nCopies <- length(prior$copy)
for(vv in c("qBucketMult", "qSfcLatRunoffMult")) {
    rst.prior <- prior[[vv]][nCopies-1-as.integer(ensMember), length(prior$time)]
    cat(paste0('assimOnly: ',vv,' output & prior match:'),
        identical(as.numeric(rst.output[[vv]]) , rst.prior),'\n')

    ## This should NOT match. The OUTPUT files are the priors
    rst.posterior <- posterior$qlink1[,nCopies-1-as.integer(ensMember), length(posterior$time)]
    cat(paste0('assimOnly: ',vv,' output & prior do not match:'),
        !identical(as.numeric(rst.output[[vv]]) , rst.posterior),'\n')
}




