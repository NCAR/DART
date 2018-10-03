PlotEnsTimeSeries <- function(pinfo) {
## PlotEnsTimeSeries: Plots time series of ensemble members, mean and truth.
#
#
# PlotEnsTimeSeries is intended to be called by 'plot_ens_time_series'
# The only input argument is a structure with model-dependent
# components.
#
# USAGE: PlotEnsTimeSeries( pinfo );
#
# STRUCTURE COMPONENTS FOR low-order models
# truth_file      name of netCDF DART file with copy tagged 'true state'
# diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
#                 and 'ensemble spread'
# var             name of netCDF variable of interest
# var_inds        indices of variables of interest

# Example 1
  if(FALSE){
    library(rwrfhydro)
    dartPath <-
      strsplit(grep('wrf_hydro_dart', readLines('~/.wrf_hydro_tools'), value=TRUE),'=')[[1]][2]
    source(paste0(dartPath,'/R/PlotEnsTimeSeries.R'))

    pinfo <- list()
    pinfo$truth_file  = 'True_State.nc'
    pinfo$diagn_file  = 'Prior_Diag.nc'
    pinfo$model       = 'wrf_hydro'      ## not really necessary currently

    pinfo$var         = 'qlink1'
    ##pinfo$var_inds  = which(trimws(ncdump("../DOMAIN/RouteLink.nc",'gages', q=TRUE))!='')
    pinfo$var_inds    = c(130,136,157)
    pinfo$var_ind_name= trimws(ncdump("../DOMAIN/RouteLink.nc",'gages', q=TRUE)[pinfo$var_inds])
    PlotEnsTimeSeries(pinfo)
       
    pinfo$var         = 'qBucketMult'
    pinfo$var_inds    = NA
    pinfo$var_ind_name= NA
    PlotEnsTimeSeries(pinfo)

    pinfo$var         = 'qSfcLatRunoffMult'
    PlotEnsTimeSeries(pinfo)

  }
  ## Code starts here -----------------------------------------------------------------

  if(!file.exists(pinfo$diagn_file)) stop(paste0("No such file: ", pinfo$diagn_file))
  pinfo$diagn_file_abs <- if(substr(pinfo$diagn_file,1,1)=='/') pinfo$diagn_file else paste0(getwd(),'/',pinfo$diagn_file)
  
  ## Query the metadata to determine which "copy" is the ensemble mean.
  ncMeta <- ncdump(pinfo$diagn_file, q=TRUE)
  varunits <- ncMeta$var[[pinfo$var]]$units
  nVarDims <- ncMeta$var[[pinfo$var]]$ndims
  varDims <- unlist(lapply(ncMeta$var[[pinfo$var]]$dim, function(vv) vv$name))
  nMembers <- ncMeta$dim$copy$len
  copyNames <- ncdump(pinfo$diagn_file,'CopyMetaData', q=TRUE)
  ens_mean_index <- which( trimws(copyNames) ==  'ensemble mean')
  whEnsMems <- grep('ensemble member',copyNames)
  gregorianOrigin <- "1601-01-01 00:00:00 UTC"
  time <- as.POSIXct( ncdump(pinfo$diagn_file,'time', q=TRUE) * 24*3600,
                      origin=gregorianOrigin, tz='UTC')

  ## JLM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ## WOW, going to need to check that this is correct via performing the calculation...
  # because of issues with ncdf4
  ## JLM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ## If the truth is known, great.
  ## JLM, deal with truth when we have/fake it
  #have_truth <- FALSE
  #if(file.exists(pinfo$truth_file)) {
  #  have_truth  = TRUE
  #  truth_index = get_copy_index(pinfo$truth_file, 'true state')
  #}

  ## This works for scalar multipliers (e.g. qBucketMult, qSfcLatRunoffMult)
  if(all(varDims==c("scalar", "copy", "time"))) {
    ens_mean <-
      ncdump(pinfo$diagn_file, pinfo$var, q=TRUE)[ens_mean_index, ,drop=FALSE]
    ens_members <-
      ncdump(pinfo$diagn_file, pinfo$var, q=TRUE)[whEnsMems, ,drop=FALSE]
    if(!(abs(mean(ens_members[,1])-ens_mean[1,1])) < 1e-7)
      warning("The ensemble mean does not match a spot check of its value.\n",
              "Please investigate, may be a dimension issue with ncdf4 package.")

    dimnames(ens_mean) <- list()
    dimnames(ens_mean)[[2]] <- format(time,'%Y-%m-%d %H:%M:%S')
    dimnames(ens_mean)[[1]] <- trimws(copyNames[ens_mean_index])
    meanDf <- plyr::adply(ens_mean, c(1,2))
    meanDf <- plyr::rename(meanDf, c(X1='member', X2='time', V1=pinfo$var))

    dimnames(ens_mean) <- list()
    dimnames(ens_members)[[1]] <- trimws(copyNames[whEnsMems])
    dimnames(ens_members)[[2]] <- format(time,'%Y-%m-%d %H:%M:%S')
    ensDf <- plyr::adply(ens_members, c(1,2))
    ensDf <- plyr::rename(ensDf, c(X1='member', X2='time', V1=pinfo$var))

    ensDf <- rbind(ensDf, meanDf)
    ensDf$time <- as.POSIXct(ensDf$time, tz='UTC')

  }

  ## This works for qlink1
  if(all(varDims==c("links", "copy", "time"))) {
    ens_mean <-
      ncdump(pinfo$diagn_file, pinfo$var, q=TRUE)[pinfo$var_inds, ens_mean_index, ,drop=FALSE]
    ens_members <-
      ncdump(pinfo$diagn_file, pinfo$var, q=TRUE)[pinfo$var_inds, whEnsMems, ,drop=FALSE]
    if(!(abs(mean(ens_members[1,,1])-ens_mean[1,1,1])) < 1e-7)
      warning("The ensemble mean does not match a spot check of its value.\n",
              "Please investigate, may be a dimension issue with ncdf4 package.")

    dimnames(ens_mean) <- list()
    dimnames(ens_mean)[[1]] <- pinfo$var_ind_name
    dimnames(ens_mean)[[2]] <- trimws(copyNames[ens_mean_index])
    dimnames(ens_mean)[[3]] <- format(time,'%Y-%m-%d %H:%M:%S')
    meanDf <- plyr::adply(ens_mean, c(1,2,3))
    meanDf <- plyr::rename(meanDf, c(X1='gage', X2='member', X3='time', V1=pinfo$var))

    dimnames(ens_members) <- list()
    dimnames(ens_members)[[1]] <- pinfo$var_ind_name
    dimnames(ens_members)[[2]] <- trimws(copyNames[whEnsMems])
    dimnames(ens_members)[[3]] <- format(time,'%Y-%m-%d %H:%M:%S')
    ensDf <- plyr::adply(ens_members, c(1,2,3))
    ensDf <- plyr::rename(ensDf, c(X1='gage', X2='member', X3='time', V1=pinfo$var))

    ensDf <- rbind(ensDf, meanDf)
    ensDf$time <- as.POSIXct(ensDf$time, tz='UTC')

  }



  ## if ( have_truth )
  ## obs?
  
  ggP <-
    ggplot2::ggplot(ensDf) +
    ggplot2::geom_line(ggplot2::aes_string(x='time', y=pinfo$var, color='member')) +
    ggplot2::ggtitle(pinfo$diagn_file_abs) +
    ggplot2::theme_bw()

  if('gage' %in% names(ensDf)) ggP <- ggP + ggplot2::facet_wrap(~gage, ncol=1, scale='free_y')


  print(ggP)
  return(invisible(list(gg=ggP, data=ensDf, pinfo=pinfo)))

}
