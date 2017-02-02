%
% Script to explore the values of the "lake" columns. 
% Since lakes use a different formulation for snow, the SNLSNO variable is insufficient
% to define if the corresponding values are indeterminate or not.
%
% DART needs to know when to replace all the indeterminate values with DART missing.
% I want to use something more robust than using snlsno == 0 to identify all columns
% that have no snow.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% add the 'easy' netcdf support

if (exist('nc_varget') ~= 2)
   addpath /contrib/matlab
   ncstartup
end

format LONGENG

% this was the file that spawned the whole issue
% snowdp = nc_varget('camforcing_6hrSolar.clm2.r.2000-12-31-00000.nc','SNOWDP');

% fname   = '/glade/proj3/cseg/inputdata/ccsm4_init/I2000CN_f09_g16_c100305/0001-01-01/I2000CN_f09_g16_c100305.clm2.r.0001-01-01-00000.nc';
fname   = 'camforcing_6hrSolar.clm2.r.2000-12-31-00000.nc';

if (exist('lakes','var') ~=2)
   ityplun          = nc_varget(fname,'cols1d_ityplun');  %(column) "land unit types" ... lakes == 3
   SNLSNO           = nc_varget(fname,'SNLSNO');          %(column) "number of snow layers"
   SNOWDP           = nc_varget(fname,'SNOWDP');          %(column) "snow depth"
   WA               = nc_varget(fname,'WA');              %(column) "water in the unconstrained aquifer"
   WT               = nc_varget(fname,'WT');              %(column) "total water storage"
   ZWT              = nc_varget(fname,'ZWT');             %(column) "water table depth"
   frac_sno         = nc_varget(fname,'frac_sno');        %(column) "fraction of ground covered by snow (0 to 1)"
   T_GRND           = nc_varget(fname,'T_GRND');          %(column) "ground temperature"
   T_GRND_U         = nc_varget(fname,'T_GRND_U');        %(column) "urban ground temperature"
   T_GRND_R         = nc_varget(fname,'T_GRND_R');        %(column) "rural ground temperature"
   DZSNO            = nc_varget(fname,'DZSNO');           %(column, levsno) "snow layer thickness"
    ZSNO            = nc_varget(fname, 'ZSNO');           %(column, levsno) "snow layer depth"
   ZISNO            = nc_varget(fname,'ZISNO');           %(column, levsno) "snow interface depth"
   T_SOISNO         = nc_varget(fname,'T_SOISNO');        %(column, levtot) "soil-snow temperature"
   T_LAKE           = nc_varget(fname,'T_LAKE');          %(column, levlak) "lake temperature"
   snw_rds          = nc_varget(fname,'snw_rds');         %(column, levsno) "effective radius" "um"
   mss_bcpho        = nc_varget(fname,'mss_bcpho');       %(column, levsno) "hydrophobic black carbon mass" "kg m-2"
   mss_bcphi        = nc_varget(fname,'mss_bcphi');       %(column, levsno) long_name = "hydrophilic black carbon mass" "kg m-2"
   mss_ocpho        = nc_varget(fname,'mss_ocpho');       %(column, levsno) "hydrophobic organic carbon mass" "kg m-2"
   mss_ocphi        = nc_varget(fname,'mss_ocphi');       %(column, levsno) "hydrophilic organic carbon mass" "kg m-2"
   mss_dst1         = nc_varget(fname,'mss_dst1');        %(column, levsno) "dust species 1 mass" "kg m-2"
   mss_dst2         = nc_varget(fname,'mss_dst2');        %(column, levsno) "dust species 2 mass" "kg m-2"
   mss_dst3         = nc_varget(fname,'mss_dst3');        %(column, levsno) "dust species 3 mass" "kg m-2"
   mss_dst4         = nc_varget(fname,'mss_dst4');        %(column, levsno) "dust species 4 mass" "kg m-2"
   flx_absdv        = nc_varget(fname,'flx_absdv');       %(column, levsno1) "flux absorption factors (direct, VIS)" "fraction"
   flx_absdn        = nc_varget(fname,'flx_absdn');       %(column, levsno1) "flux absorption factors (direct, NIR)" "fraction"
   flx_absiv        = nc_varget(fname,'flx_absiv');       %(column, levsno1) "flux absorption factors (diffuse, VIS)" "fraction"
   flx_absin        = nc_varget(fname,'flx_absin');       %(column, levsno1) "flux absorption factors (diffuse, NIR)"
   qflx_snofrz_lyr  = nc_varget(fname,'qflx_snofrz_lyr'); %(column, levsno) "ice freezing rate" "kg m-2 s-1"

   lakes = find(ityplun == 3); % These are the colums that are lakes.
   lakeswithsnow   = find((ityplun == 3) & (SNOWDP >  0));
   lakeswithnosnow = find((ityplun == 3) & (SNOWDP <= 0));
   fprintf('\nThe number of lake snow depths >  0 is %d\n',length(lakeswithsnow))
   fprintf('The number of lake snow depths <= 0 is %d\n',length(lakeswithnosnow))
   fprintf('The number of lakes                 is %d\n',length(lakes))
end

%% summarize

vars1d = {'SNLSNO', 'SNOWDP', 'WA', 'WT', 'ZWT', ...
          'frac_sno', 'T_GRND', 'T_GRND_U', 'T_GRND_R'};

columninds = 1:length(SNLSNO);

for i = 1:length(vars1d)

   long_name = nc_attget(fname,vars1d{i},'long_name');
   units     = nc_attget(fname,vars1d{i},'units');

   subplot(3,1,1)

      eval(sprintf('vals = %s(lakes);',vars1d{i}))
      plot(columninds(lakes),vals,'*');
      h = title(sprintf('%s %s over lakes',long_name,vars1d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars1d{i},min(vals(:)),max(vals(:)))

      uniquevals = unique(vals);
      uniquevals(1:min([10 length(uniquevals)]))

   subplot(3,1,2)

      eval(sprintf('vals = %s(lakeswithsnow);',vars1d{i}))
      plot(columninds(lakeswithsnow),vals,'*');
      h = title(sprintf('%s %s over lakes WITH snow',long_name,vars1d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars1d{i},min(vals(:)),max(vals(:)))

   subplot(3,1,3)

      eval(sprintf('vals = %s(lakeswithnosnow);',vars1d{i}))
      plot(columninds(lakeswithnosnow),vals,'*');
      h = title(sprintf('%s %s over lakes with NO snow',long_name,vars1d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars1d{i},min(vals(:)),max(vals(:)))

   disp('pausing, hit any key to continue ...'), pause
end

%% exploring the 2D variables that are be snow-related

vars2d = { 'DZSNO', 'ZSNO', 'ZISNO', 'T_SOISNO', 'T_LAKE', ...
           'snw_rds', 'mss_bcpho', 'mss_bcphi', 'mss_ocpho', 'mss_ocphi', ...
           'mss_dst1', 'mss_dst2', 'mss_dst3', 'mss_dst4', 'flx_absdv', ...
           'flx_absdn', 'flx_absiv', 'flx_absin', 'qflx_snofrz_lyr'};

for i = 1:length(vars2d)

   long_name = nc_attget(fname,vars2d{i},'long_name');
   units     = nc_attget(fname,vars2d{i},'units');

   subplot(3,1,1)

      eval(sprintf('vals = %s(lakes,:);',vars2d{i}))
      plot(columninds(lakes),vals,'*');
      h = title(sprintf('%s %s over lakes',long_name,vars2d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars2d{i},min(vals(:)),max(vals(:)))

      uniquevals = unique(vals);
      uniquevals(1:min([10 length(uniquevals)]))

   subplot(3,1,2)

      eval(sprintf('vals = %s(lakeswithsnow,:);',vars2d{i}))
      plot(columninds(lakeswithsnow),vals,'*');
      h = title(sprintf('%s %s over lakes with snow',long_name,vars2d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars2d{i},min(vals(:)),max(vals(:)))

   subplot(3,1,3)

      eval(sprintf('vals = %s(lakeswithnosnow,:);',vars2d{i}))
      plot(columninds(lakeswithnosnow),vals,'*');
      h = title(sprintf('%s %s over lakes with NO snow',long_name,vars2d{i}));
      set(h,'Interpreter','none')
      ylabel(units)
      xlabel('column indices'); axis([1 length(SNLSNO) -Inf Inf])
      fprintf('\nRange of %s is %d to %d\n',vars2d{i},min(vals(:)),max(vals(:)))

   disp('pausing, hit any key to continue ...'), pause

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
