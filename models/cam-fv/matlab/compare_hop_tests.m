function compare_hop_tests(name1, file1, name2, file2, varname)
% DART compare_hop_tests 
%
% cd /glade/scratch/thoar/
% ncdiff cam_hop_proposed_redo/run/cam_hop_proposed_redo.cam.i.2004-01-11-21600.nc \
%                  cam_hop_proposed/run/cam_hop_proposed.cam.i.2004-01-11-21600.nc \
%        FT_minus_onehop.nc
%
% ncdiff                       cam_hop_FF/run/cam_hop_FF.cam.i.2004-01-11-21600.nc \
%                  cam_hop_proposed/run/cam_hop_proposed.cam.i.2004-01-11-21600.nc \
%        FF_minus_onehop.nc
%
% USAGE:
%
% name1 = 'FT - onehop'
% file1 = '/glade/scratch/thoar/camfv_hop/rs_minus_onehop.nc';
% name2 = 'FF - onehop'
% file2 = '/glade/scratch/thoar/camfv_hop/FF_minus_onehop.nc';
%
% compare_hop_tests(name1, file1, name2, file2)
%
% -or-
%
% compare_hop_tests(name1, file1, name2, file2,{'PS','T','US','VS','Q','CLDICE','CLDLIQ'})


%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( (exist(file1,'file') ~= 2) || (exist(file2,'file') ~= 2) )
   fprintf('%s, \n',file1)
   fprintf('%s, or\n',file2)
   error('does not exist.')
end

%% extract all the variable names from both files and find the intersection.
vars1 = sort(FindProgVars(file1));
vars2 = sort(FindProgVars(file2));
vars  = sort(FindCommonVars(vars1,vars2));

if (nargin > 4)
    if (iscell(varname))
       vars = varname;
    else
       vars = {varname};
    end
    pausecmd = 'fprintf(''pausing at level %d ... hit a key to continue\n'',levelindx); pause';
else
    pausecmd = 'fprintf(''           level %d ... \n'',levelindx); pause(1.0)';
end

%% loop over the variables and extract them from both files.
% Since Matlab automatically squeezes out singleton dimensions, and 'time'
% is usually a singleton in these files, I need to squeeze the singleton
% dimension out of the variable dimension strings, too.

bob = jet128; colormap(bob);

for i = 1:length(vars)

   varname = vars{i};

   switch (varname(1))
      case{'n','p','s','d'}
         fprintf('Skipping %s\n',varname)

      otherwise
         varinfo       = nc_getvarinfo(file1,vars{i});
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
         hop.units     = GetAttribute(file1,vars{i},'units');
      
         hop.name1     = name1;
         hop.data1     = nc_varget(file1,vars{i});
         hop.name2     = name2;
         hop.data2     = nc_varget(file2,vars{i});
         hop.change    = abs(hop.data2) - abs(hop.data1);
      
         switch BestPlotType(hop)
            case 'bylevel'
               for ilevel = max(1,15):hop.dimsizes(hop.levdim)
                  myplot(hop, ilevel, pausecmd)
               end  
            case 'horizontalslab'
               my2dplot(hop)
            case 'lineplot'
            otherwise
         end
   end
end


function myplot(hopobj,levelindx,pausecmd)
%% Make some plots

slab    = squeeze(hopobj.change(levelindx,:,:));
slabmin = min(slab(:));
slabmax = max(slab(:));

if slabmin == slabmax, return; end

datmax    = max(abs([slabmin slabmax]));
clim      = [-datmax datmax];
histedges = FindHistogramEdges(clim,101);

test1     = squeeze(hopobj.data1(levelindx,:,:));
test2     = squeeze(hopobj.data2(levelindx,:,:));
slab2min  = min(min(test1(:)), min(test2(:)));
slab2max  = max(max(test1(:)), max(test2(:)));
datmax2   = max(abs([slab2min slab2max]));
clim2     = [-datmax2 datmax2];

sbpos = [0.10 0.050 0.80 0.225; 
         0.10 0.325 0.80 0.225;
         0.10 0.600 0.80 0.225;
         0.10 0.900 0.80 0.050];

figure(1); clf; orient tall; 
subplot('position',sbpos(4,:))
   [n,bin] = histc(slab(:),histedges);
   bar(histedges, n, 'histc');
   str1 = sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)', ...
                    slabmin, hopobj.units, slabmax, hopobj.units);
   title(str1,'Interpreter','none')
   xlabel(sprintf('"%s" positive values means "restart" is better',hopobj.units))

subplot('position',sbpos(3,:))
   imagesc(test1,clim2);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel(sprintf('%s level %d %s',hopobj.varname,levelindx,hopobj.name1));
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

subplot('position',sbpos(2,:))
   imagesc(test2,clim2);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel(sprintf('%s level %d %s',hopobj.varname,levelindx,hopobj.name2));
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

subplot('position',sbpos(1,:))
   imagesc(slab,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel('nosurfrest - surfrest');
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

   eval(pausecmd)


function my2dplot(hopobj)
%% Make some plots
%

slabmin = min(hopobj.change(:));
slabmax = max(hopobj.change(:));

if slabmin == slabmax
    fprintf('no difference in %s\n',hopobj.varname)
    return
end

datmax    = max(abs([slabmin slabmax]));
clim      = [-datmax datmax];
histedges = FindHistogramEdges(clim,101);

test1     = hopobj.data1(:,:);
test2     = hopobj.data2(:,:);
slab2min  = min(min(test1(:)), min(test2(:)));
slab2max  = max(max(test1(:)), max(test2(:)));
datmax2   = max(abs([slab2min slab2max]));
clim2     = [-datmax2 datmax2];

sbpos = [0.10 0.050 0.80 0.225; 
         0.10 0.325 0.80 0.225;
         0.10 0.600 0.80 0.225;
         0.10 0.900 0.80 0.050];

% Need to know what dimensions we are plotting

figure(1); clf; orient tall; 

subplot('position',sbpos(4,:))
   [n,bin] = histc(hopobj.change(:),histedges);
   bar(histedges, n, 'histc');
   str1 = sprintf('(min %0.5g %s) hopping difference histogram (max %0.5g %s)', ...
                   slabmin, hopobj.units, slabmax, hopobj.units);
   title(str1,'Interpreter','none')
   xlabel(sprintf('"%s" positive values means "restart" is better',hopobj.units))

subplot('position',sbpos(3,:))
   imagesc(test1,clim2);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel(sprintf('%s %s',hopobj.varname,hopobj.name1));
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

subplot('position',sbpos(2,:))
   imagesc(test2,clim2);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel(sprintf('%s %s',hopobj.varname,hopobj.name2));
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

subplot('position',sbpos(1,:))
   imagesc(hopobj.change,clim);
   set(gca,'YDir','normal','TickDir','out','XMinorTick','on','FontSize',14)
   ylabel('nosurfrest - surfrest');
   axis image
   h = colorbar('vert'); set(get(h,'YLabel'),'String',hopobj.units)

disp('pausing - hit any key to continue ...'); pause
   

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


function x = FindHistogramEdges(minmax,nedges)

x = minmax(1) + [1:nedges]*(minmax(2)-minmax(1))/(nedges-1);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
