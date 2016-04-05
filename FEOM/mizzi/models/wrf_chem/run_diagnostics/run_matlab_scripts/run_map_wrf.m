fname_diff         = '/glade/p/acd/mizzi/DART_TEST_AVE/MOPCOMB_Exp_3_RtDA_60M_p30p30_sp4/Cat_Diff.nc';
fname_post         = '/glade/p/acd/mizzi/DART_TEST_AVE/MOPCOMB_Exp_3_RtDA_60M_p30p30_sp4/Cat_Posterior_Diag.nc';
fname_prior        = '/glade/p/acd/mizzi/DART_TEST_AVE/MOPCOMB_Exp_3_RtDA_60M_p30p30_sp4/Cat_Prior_Diag.nc';
%
varname       = 'co_d01';
copystring    = 'ensemble mean';
levelindx     = 2;
timeindx      = 1;
map_wrf(fname_diff, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc diff_10_06;
%
timeindx      = 3;
map_wrf(fname_diff, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc diff_10_18;
%
timeindx      = 1;
map_wrf(fname_post, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc post_10_06;
%
timeindx      = 3;
map_wrf(fname_post, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc post_10_18;
%
timeindx      = 1;
map_wrf(fname_prior, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc prior_10_06;
%
timeindx      = 3;
map_wrf(fname_prior, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc prior_10_18;
%
timeindx      = 4;
map_wrf(fname_prior, varname, levelindx, timeindx, copystring);
worldmap;
axis off;
print -dpsc prior_11_00;
