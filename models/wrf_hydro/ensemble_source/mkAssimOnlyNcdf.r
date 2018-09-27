## the mkNcdf routine is in my path (src looks there)
## please ask for it if you want it, it's really useful.
src('mkNcdf.r')
varList = list()

varList[[1]] <- list( name='precipMult',
                      longname='Precipitation Multiplier',
                      units='-',
                      precision = 'double',
                      missing = -9999,
                      dimensionList =
                           list(
                                scalar=list(name='scalar',values=1,
                                  units='-', unlimited=FALSE,
                                  create_dimvar=FALSE)
                                ),
                      data = 1:1 )
globalAttList <- list()
globalAttList[[1]] <- list(name='Restart_Time',value='2012-07-05_00:00:00',
                           precision="text" )

## hell, make an ensemble here... for now
varList[[1]]$data <- .8
dum <- mkNcdf( varList, globalAttList, 'restart.assimOnly.0001.nc' )
varList[[1]]$data <- 1
dum <- mkNcdf( varList, globalAttList, 'restart.assimOnly.0002.nc' )
varList[[1]]$data <- 1.3
dum <- mkNcdf( varList, globalAttList, 'restart.assimOnly.0003.nc' )
