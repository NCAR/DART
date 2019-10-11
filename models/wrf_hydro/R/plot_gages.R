library(rwrfhydro)
library(data.table)
#library(ggplotly)
library(ggplot2)

# Config
## Set the path to the current run.
runDir <- "/glade/scratch/jamesmcc/wrfhydro_dart/sixmile/runs/test1/member_000"
runTag <- ''
# Add path to observations... obs_seq?

## -------------------------------------------------------
setwd(runDir)

rlFile <- grep('route_link_f',readLines("hydro.namelist"), ignore.case=TRUE, value=TRUE)
rlFile <- strsplit(rlFile,'(=|\'|\")')[[1]]
rlFile <- rlFile[length(rlFile)]

rl <- GetNcdfFile(rlFile, q=TRUE, var='time', exclude=TRUE)
theLinks <- rl$link[which(trimws(rl$gages) != '')]
theGages <- rl$gages[which(trimws(rl$gages) != '')]
link2gage <- theGages
names(link2gage) <- theLinks

exclVarsChrt <- c('reference_time','time')

GetAllChrtout <- function(path, runId) {
   # Identify the files
   modelChrtoutFiles <- list.files(path=path, pattern='CHRTOUT', full.names=TRUE)
   # Get the data using plyr. This can easily be parallelized on most systems but not shown here.
   modelData <-  
     plyr::ldply(rwrfhydro::NamedList(modelChrtoutFiles),
                 rwrfhydro::GetNcdfFile, variables=exclVarsChrt, exclude=TRUE, quiet=TRUE)
   # Work with it as a data.table.
   modelData <- data.table::as.data.table(modelData)
   # Calculate the times using POSIXct
   modelData[, `:=`(POSIXct=as.POSIXct(basename(.id), format='%Y%m%d%H', tz='UTC'))]
   # Return is the last expression
   modelData
 }

chrtData <- GetAllChrtout(path='.', runId='modeled')

chrtData.gages <- chrtData[ feature_id %in% theLinks, ]
chrtData.gages$site_no <-
     trimws(link2gage[as.character(chrtData.gages$feature_id)])
setnames(chrtData.gages, 'streamflow', 'Discharge (cms)')

# Simple way to format the y-axis labels
num2Decimals <- function(x, nPlaces=2) round(x*10^nPlaces)/(10^nPlaces)

SiteLabeller <- function(site) {
  labels <- paste0(attributes(obsDischarge)$siteInfo$station_nm, 
                   ' (',attributes(obsDischarge)$siteInfo$site_no,')')
  names(labels) <- attributes(obsDischarge)$siteInfo$site_no
  labels[site]
}

#obsModAssimColors <- c(observed='black', modeled='royalblue3', assim='tomato3')

#ggplotly(
gg <-
  ggplot(chrtData.gages) +
  geom_point(aes(x=POSIXct, y=`Discharge (cms)`), size=.5) + 
  facet_wrap(~site_no, ncol=1) + #, 
  #labeller=labeller(site_no=SiteLabeller)) +
  scale_y_continuous(trans='log', labels=num2Decimals, name='Discharge (cms)') +
  scale_x_datetime(name='2013') + 
  #scale_color_manual(name='', values=obsModAssimColors) +
  theme_bw(base_size=13) +
  theme(axis.text.y = element_text(angle = 30, size = 10)) 
