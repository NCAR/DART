function hop_test_check(file0, file1, file2, varname)
% DART hop_test_check 
%
% USAGE:
% hop_test_check(file0, file1, file2);
%
% file0    is the filename of the initial state 
% file1    is the filename of the 'long hop' case
% file2    is the filename of the 'short hop' case
%
% The difference between file2 and file1 will be scaled by the 
% time tendency - the amount the field changed between file0 and file1.
%
% Example:
%
% file0 = '/gpfs/ptmp/thoar/restarts/CAM/caminput_1.nc';
% file1 = '/glade/scratch/thoar/archive/hop_24/rest/2008-11-01-00000/hop_24.cam_0001.i.2008-11-01-00000.nc';
% file2 = '/glade/scratch/thoar/archive/hop_12/rest/2008-11-01-00000/hop_12.cam_0001.i.2008-11-01-00000.nc';
% hop_test_check(file0, file1, file2)
%
% Example:
%
% file3 = '/gpfs/ptmp/thoar/restarts/CLM/clminput_1.nc';
% file4 = '/glade/scratch/thoar/archive/hop_24/rest/2008-11-01-00000/hop_24.clm2_0001.r.2008-11-01-00000.nc';
% file5 = '/glade/scratch/thoar/archive/hop_12/rest/2008-11-01-00000/hop_12.clm2_0001.r.2008-11-01-00000.nc';
% hop_test_check(file3, file4, file5)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( (exist(file0,'file') ~= 2) || ...
     (exist(file1,'file') ~= 2) || ...
     (exist(file2,'file') ~= 2) )
   fprintf('One of the input files does not exist.\n')
   fprintf('%s, \n',file0)
   fprintf('%s, or\n',file1)
   error('%s',file2)
end

%% extract all the variable names from both files and find the intersection.
vars1 = sort(FindProgVars(file1));
vars2 = sort(FindProgVars(file2));
vars  = sort(FindCommonVars(vars1,vars2));

%% loop over the variables and extract them from both files.
%  take the difference and convert it to a percentage of 'file1-file0' ... 
% Since Matlab automatically squeezes out singleton dimensions, and 'time'
% is usually a singleton in these files, I need to squeeze the singleton
% dimension out of the variable dimension strings, too.

bob = jet128; colormap(bob);

if (nargin > 3)
    vars = {varname};
    pausecmd = 'fprintf(''pausing at level %d ... hit a key to continue\n'',ilevel); pause';
else
    pausecmd = 'fprintf(''           level %d ... \n'',ilevel); pause(1.0)';
end

for i = 1:length(vars)

   varinfo       = nc_getvarinfo(file0,vars{i});
   nonsingletons = (varinfo.Size > 1);
   if (varinfo.Nctype == 2)
       % Character string variables need not be checked.
       fprintf('Skipping   %s\n',vars{i})
       continue
   else
       fprintf('Comparing  %s\n',vars{i})
   end
   hop.varname   = vars{i};
   hop.dimsizes  = varinfo.Size(nonsingletons);
   hop.dimnames  = varinfo.Dimension(nonsingletons);
   hop.levdim    = find(strcmp(hop.dimnames,'lev'));

   twohop        = nc_varget(file2,vars{i});
   onehop        = nc_varget(file1,vars{i});
   start         = nc_varget(file0,vars{i});
   hop.units     = GetAttribute(file0,vars{i},'units');
   hop.tendency  = onehop - start;
   hop.change    = twohop - onehop;

   switch BestPlotType(hop)
      case 'bylevel'
         for ilevel = 1:hop.dimsizes(hop.levdim)
            myplot(hop, ilevel)
            eval(pausecmd)
         end  
      case 'horizontalslab'
         my2dplot(hop)
         disp('          pausing - hit any key to continue ...')
         pause
      case 'lineplot'
      otherwise
   end
 
end


function x = BestPlotType(hopobj)

   if ( ~isempty(hopobj.levdim) )
       x = 'bylevel';
   elseif (length(hopobj.dimsizes) == 1)
       x = 'lineplot';
   elseif (isempty(hopobj.levdim))
       x = 'horizontalslab';
   else
       x = 'dunno';
   end


function my2dplot(hopobj)
%% Make some plots
%

slabmin = min(hopobj.change(:));
slabmax = max(hopobj.change(:));

orgmin = min(hopobj.tendency(:));
orgmax = max(hopobj.tendency(:));
datmax = max(abs([orgmin orgmax]));

if orgmin == orgmax
    clim = [-1 1];
else
    clim = [-datmax datmax];
end
histedges = FindHistogramEdges(clim,101);

sbpos = [0.10 0.06 0.80 0.28; 
         0.10 0.43 0.80 0.28;
         0.10 0.80 0.80 0.15]; 

% Need to know what dimensions we are plotting

figure(1); clf; orient tall; 

subplot('position',sbpos(3,:))
   histc(hopobj.change(:),histedges)
   str1 = sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)', ...
                   slabmin, hopobj.units, slabmax, hopobj.units);
   title(str1,'Interpreter','none')
   xlabel(hopobj.units)

subplot('position',sbpos(2,:))
   imagesc(hopobj.change,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   str2 = sprintf('%s difference from hopping',hopobj.varname);
   title(str2,'Interpreter','none')
   axis image
   ylabel(sprintf('%s (index)',hopobj.dimnames{1}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',hopobj.units)
   
subplot('position',sbpos(1,:))
   imagesc(hopobj.tendency,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   str3 = sprintf('%s tendency',hopobj.varname);
   title(str3,'Interpreter','none')
   axis image
   ylabel(sprintf('%s (index)',hopobj.dimnames{1}))
   xlabel(sprintf('%s (index)',hopobj.dimnames{2}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',hopobj.units)


function myplot(hopobj,levelindx)
%% Make some plots
%

slab    = squeeze(hopobj.change(levelindx,:,:));
slabmin = min(slab(:));
slabmax = max(slab(:));

if slabmin == slabmax, return; end

orgslab = squeeze(hopobj.tendency(levelindx,:,:));
orgmin  = min(orgslab(:));
orgmax  = max(orgslab(:));
datmax  = max(abs([orgmin orgmax]));
clim    = [-datmax datmax];
histedges = FindHistogramEdges(clim,101);

sbpos = [0.10 0.06 0.80 0.28; 
         0.10 0.43 0.80 0.28;
         0.10 0.80 0.80 0.15]; 

figure(1); clf; orient tall; 
subplot('position',sbpos(3,:))
   histc(slab(:),histedges)
   str1 = sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)', ...
                    slabmin, hopobj.units, slabmax, hopobj.units);
   title(str1,'Interpreter','none')
   xlabel(hopobj.units)

subplot('position',sbpos(2,:))
   imagesc(slab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   str2 = sprintf('%s level %d difference from hopping',hopobj.varname,levelindx);
   title(str2,'Interpreter','none')
   axis image
   ylabel(sprintf('%s (index)',hopobj.dimnames{2}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',hopobj.units)
   
subplot('position',sbpos(1,:))
   imagesc(orgslab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   str3 = sprintf('%s level %d tendency',hopobj.varname,levelindx);
   title(str3,'Interpreter','none')
   axis image
   ylabel(sprintf('%s (index)',hopobj.dimnames{2}))
   xlabel(sprintf('%s (index)',hopobj.dimnames{3}))
   h = colorbar('vert');
   set(get(h,'YLabel'),'String',hopobj.units)



function vars = FindProgVars(fname)
%% Determine the multi-dimensional variables that are not coordinate variables.
% These are probably the variables of interest.

fileinfo  = nc_info(fname);
nvars     = length(fileinfo.Dataset);
isPROGvar = zeros(nvars,1);

% We are not interested in the 1D variables, so just skip them.
for i = 1:nvars
   if ( length(fileinfo.Dataset(i).Size) > 1 ), isPROGvar(i) = 1; end
end

if (sum(isPROGvar) == 0)
   error('No multidimensional variables in %s',fname)
end

% coerce just the names into a cell array 

varind = 0;
for i = 1:nvars
   if (isPROGvar(i) > 0)
      varind = varind + 1;
      vars{varind} = fileinfo.Dataset(i).Name;
   end
end



function vars = FindCommonVars(vars1,vars2)
% 
k = 0;

for i = 1:length(vars1)
   if ( any(strcmp(vars1{i},vars2)) )
      k = k + 1;
      vars{k} = vars1{i};
   end
end


function x = jet128()
x = jet(128);
x(64:65,:) = 1;

function attvalue = GetAttribute(fname,varname,attname)
varinfo = nc_getvarinfo(fname,varname);
attvalue = [];
for i = 1:length(varinfo.Attribute)
   if (strcmp(varinfo.Attribute(i).Name,attname))
      attvalue = varinfo.Attribute(i).Value;
      break
   end
end


function x = FindHistogramEdges(minmax,nedges)

x = minmax(1) + [1:nedges]*(minmax(2)-minmax(1))/(nedges-1);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
