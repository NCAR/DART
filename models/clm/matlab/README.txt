# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


Contents: Instructions to help visualize the Prior/Posterior/True_State netCDF files
using clm_get_var.m and clm_plot_var.m scripts.  These scripts reconstitute vector
formatted netCDF files (e.g. CLM restart, DART stage files) such that they
can be viewed in lon/lat format

*) to check if I preserved the metadata of the CLM variables (given their 'sparse' representation),
   I used clm_get_var.m  to read a variable from the DART Posterior/Prior/True_State netCDF files:

   fname      = 'clm_restart.nc';
   varname    = 'frac_sno';
   copystring = 'ensemble member 1';
   levelindex = 1;
   timeindex  = 1;
   x = clm_get_var(fname,varname,copystring,levelindex,timeindex);

   % clm_get_var uses the 'landfrac' variable as the ocean mask.
   % anything with a landfrac of zero is assumed to be ocean and is
   % set to Matlab's "missing value" NaN.

   h = imagesc(x.lon, x.lat, x.datmat);
   set(h,'AlphaData',~isnan(x.datmat));  % sets NaNs to the background color
   set(gca,'YDir','normal');
   h = title(x.varname);
   set(h,'Interpreter','none');
   worldmap;
   colorbar;



*) to check the localization ... it is useful to create a single observation with the
   standard DART tools: create_obs_sequence, create_fixed_network, perfect_model_obs ...
   Given that observation, severely restrict the area impacted by the observation by
   setting the localization to a small number ...  &assim_tools_nml:cutoff = 0.05
   is just a few gridpoints at the 1-degree resolution of CLM. 
   This should make the half-width (about) 0.05radians * 40000km/(2pi radians) = 318+ km

   To ensure the observation is impacting the appropriate part of the state vector,
   it is easy to difference the Prior and Posterior diagnostic files. The problem is
   that the metadata is then zero (it is the same in both files, naturally), so it
   is convenient to differ the files and then add back the metadata from one of the
   parents:

   Example 1: (CLM 4.5 and prior versions)
   ncdiff analysis.2000-01-06-00000.nc preassim.2000-01-06-00000.nc Innov.nc
   ncks -A -v lon,lat,levgrnd,area,landfrac,cols1d_ixy,cols1d_jxy,pfts1d_ixy,pfts1d_jxy,cols1d_wtxy,pfts1d_wtxy preassim.2000-01-06-00000.nc  Innov.nc

   Example 2: (CLM5.0 and later versions) 
   ncdiff clm_analysis_member_0001_d01.2011-01-02-00000.nc clm_preassim_member_0001_d01.2011-01-02-00000.nc Innov.nc
   ncks -A -v lon,lat,levgrnd,levsoi,levdcmp,levlak,area,landfrac,land1d_ixy,land1d_jxy,land1d_wtxy,land1d_ityplun,cols1d_ixy,cols1d_jxy,cols1d_lon,cols1d_lat,cols1d_ityplun,pfts1d_ixy,pfts1d_jxy,cols1d_wtxy,pfts1d_wtxy,pfts1d_lon,pfts1d_lat,pfts1d_ityplun clm_analysis_member_0001_d01.2011-01-02-00000.nc  Innov.nc

   % Now - visualize as before, with one added twist. 
   % To differentiate between ocean and a difference of zero,
   % plot everything with a difference of zero as a light gray color,
   % while still plotting ocean colors (NaNs) as the background color
   % (normally, white).

   fname      = 'Innov.nc';
   varname    = 'frac_sno';
   copystring = 'ensemble member 1';
   levelindex = 1;
   timeindex  = 1;
   x = clm_get_var(fname,varname,copystring,levelindex,timeindex);

   h = imagesc(x.lon, x.lat, x.datmat);
   set(h,'AlphaData',~isnan(x.datmat));
   set(gca,'YDir','normal');
   h = title(sprintf('Innovations in %s',x.varname));
   set(h,'Interpreter','none');
   worldmap;
   colorbar;
   bob = colormap;
   bob(1,:) = 0.9;
   colormap(bob);

   I looked at the obs_seq.final to manually grab the location:

   loc3d = [1.527163095495038        0.5125752760427027         0.000000000000000     -1 ]*180/pi
   (which results in)
   loc3d = 87.5000   29.3684         0  -57.2958
   hold on;
   plot(87.5, 29.3684,'o','MarkerSize',20)
   plot(87.5, 29.3684,'x','MarkerSize',20);
   print -dpng CLM_cutoff_0.05.png             


*) It will be useful to perform these single-observation localization checks at 
   various parts of the DART state vector ... at different levels, for different
   state variables ... 


*) The variables that use levtot or levsno are problematic. Each levsno depth is
different from one ensemble member to the next, so there is no single coordinate variable to
write out. I'd have to read "zsno(column,levsno)" for each ensemble member (but not part 
of the state vector?) and do something with it. For the time being, I am just going to set 
all the snow depth values (the midpoints of the snow layers, or the surface of the snow layers)
to a dumb value that is obviously not right. This only pertains to the 'coordinate' values,
the state values/prognostic variable values will be correct.

TJH Thu Aug 25 14:42:28 MDT 2011

I have to scan the incoming values from the CLM restart files and replace them with MISSING
values (or perhaps the value that causes the machine to STOP). T_LAKE is defined (column,levlak)
even though 90+% of the columns are not lakes! At present, the values in those columns is undefined,
there is an index array that specifies the portion of the columns with valid values. I don't want
to get into it at that kind of detail ... too many decisions to track, too many variables to support.

