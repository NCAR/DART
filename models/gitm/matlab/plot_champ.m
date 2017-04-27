function plot_champ(dart_install_dir, truth_run_dir, pbs_file_t, middle_run_dir, pbs_file_m, dart_run_dir, pbs_file_d)
%% This file plots data from
% DART (preassim.nc, analysis.nc), 
% truth simulation and middle simulation from data files 3DALL_*.b00*), and 
% GITM satellite .dat files created by herot.engin.umich.edu/bigdisk1/bin/makerun_wr.pl
%
% WHAT are TRUTH and MIDDLE simulations? 
% Imagine a hypothetical situation:
% you have an ensemble (with just 2 members) with the only thing differing
% for these members is f107 (129 and 131 respectively). What I do is I run
% assimilation on these members, but also run 2 more simulations: truth
% (which uses f107 as was measured by NOAA on that day - say 175.3) and a
% MIDDLE SIMULATION (which uses f107 of 130 - a middle (or mean) of all
% ensemble members' f107s and IS NOT subject to assimilation!
% Hence it is called the middle simulation - a
% simulation representative of the middle of ensemble if it WASNOT
% assimilated upon (also sometimes called GITM without DART).
%
% INPUTS: 
%
% -dart_install_dir - where is DART installed on this computer? (used to read the champ(_with_rho.dat) and grace files in the 
%   models/gitm/GITM2/run/ directory; f107.txt in models/gitm/GITM2/srcData 
%   directory; and use the text_GITM converter in observations/text_GITM/work)
%   My paths-naming convention is that when I store path in a variable, I do not terminate it with
%   a slash (so that when I append something to it, I need to put an extra
%   slash). Just a convention. 
%
% -truth_run_dir - where are the GITM files for the truth run (data/3DALL_*.b00* and data/cham_*.sat and pbs_file)
%
% -pbs_file_t - name of the run script used to simulate truth (inside truth_run_dir, or can be ../../pbs_file.sh)
%
% -middle_run_gitm - what middle-gitm run data files you want to plot (don't have to be in ${wd}) (as in 3DALL_t*.b0001) #the middle run is different from the true one because of f107 it uses - I want to see what "middle of ensembe" would do if I just let it run without adjusting it. For example, if truth is 150 and ensemble is centered about 130, then I want to run three simulations: the truth (without assimilation, f=150), the middle (without assimilation, f=130) and the ensemble with assimilation (f is adjusted on the fly). right now middle_run_gitm is set to be the same as truth_run_dir, but that's only to check if all the files are in the right places - change this when you do the actual runs.
%
% -pbs_file_m - name of the pbs file relative to $middle_run_gitm, which containes grid information (like nblockslon and nlons) and awt
%
% -dart_run_dir - where are the netcdf files for the ensemble itself (P*_Diag.nc and pbs_file)
%
% -pbs_file_d - name of the run script used to simulate the ensemble (inside dart_run_dir, or can be ../../pbs_file.sh)
%
% 
% OUTPUTS
%  files in pdf format in the matlab directory
%
% Examples
%
% Alexey Morozov alexeymor at google mail 08/15/2012

% DART $Id$
% CREDIT: Alexey Morozov


%%% Variable-names convention (examples of time, coordinate, and density)
%if a variable ends in T, this variable comes from GITM truth simulation,
% ie LonT is the array of longitude coordinates coming from GITM binary data
% files (here usually 3DALL_*.b00*).
% tt is time gitm - ie time coming from the b00* files (ie times when you have truth data
% available)
%if a variable ends in M, it comes from gitm Middle simulation (read line 6)
%if a variable ends in D, it comes from DART netcdf (*.nc) files
%if a variable ends in C, it comes from GITM satellite champ.dat file
%if a variable ends in G, it comes from GITM satellite grace.dat file
%else, it is a variable that is a combination or extension of the above

%Times (usually MINUTES, sometimes converted to hours or days...):
% tt              min, GITM truth time (GITM writes out data files at prespecified intervals (see UAM.in:#SAVEPLOTS), so usually this is something like 30:30:1440) hopefully tt=tm (middle time - see line 6
% td              min, DART time (DART skips an assimilation step if there were no observations in that window, so usually this is a sparse version of GITM time)
% tc              min, CHAMP time (CHAMP took measurements on avarage every 47 seconds with min of 10 seconds and
% tg              min, GRACE time 

%Longitudes (latitudes are very similar):
% LonT            deg, comes from *.b0001  (GITM  time, tt) %Longitude in middle simulations assumed to be equal to LonT
% LonD            deg, comes from *.nc     (DART  time, td)
% LonC          deg, comes from *.dat    (CHAMP time, tc)
% LonCT          deg, is LonC lin-interpolated to tt
% LonCD          deg, is LonC lin-interpolated to td
% lona            index of LonD closest to LonC at each time instant in tc
% wlon            index of LonD closest to Ann Arbor, MI

%Altitudes:
% AltT            m, comes from *.b0001  (GITM  time, tt)
% AltD            m, comes from *.nc     (DART  time, td)
% AltC          m, comes from *.dat    (CHAMP time, tc)
% AltCT          m, is AltC lin-interpolated to tt
% AltCD          m, is AltC lin-interpolated to td
% alta            index of AltD closest to AltC at each time instant in tc
% walt            index of AltD closest to Ann Arbor, MI

%Densities
% RhoT            truth density at gridpoints, comes from GITM_T 
% RhoM            middle density at gridpoints, comes from GITM_M
% RhoDr           density at gridpoints, comes from DART, pRior
% RhoDo           density at gridpoints, comes from DART, pOsterior
% RhoC          density at CHAMPloc-n, comes from CHAMP
% RhoCD          is RhoC lin-interpolated-in-time to td
% RhoCT          is RhoC lin-interpolated-in-time to tt
% RhoDi           is RhoDo     lin-interpolated-in-time-lon-lat-alt to tc-LonC2-LatC2-AltC2
% RhoDib          is RhoDo     lin-interpolated-in-     lon-lat-alt to    LonCD-LatCD-AltCD
% RhoTC         is RhoT      lin-interpolated-in-     lon-lat-alt to    LonCT-LatCT-AltCT
% RhoA           is RhoT(wlon,wlat,walt,:) (A for Ann Arbor)


%% Clear
p=pwd 
addpath(p) %so that matlab always has access to .m files in this dir no matter where it cd's
format compact
set(0,'DefaultAxesFontSize',12)

%% Load GITM truth data using 3DALL files
kind = 37; %kind of obs (used only if you're making perfect model files from scratch (advanced feature, not used in the initial tests, used in mat2dart))
kind_sd = 2.6e-13; %standard deviation of obs, obs uncertainty (used only with kind above when writing input files to DART)

[~,awt] = unix(['grep "awt=" < ' truth_run_dir '/' pbs_file_t ' | awk ''{print $1}''']); i=find(awt=='='); awt=str2double( awt((i+1):(end-1)));    %what assimilation window did you use, in seconds, pbs_*.sh:line 23
awt=awt/60; %convert aw to minutes
y2 = 2;    %when did you start assimilation, mod(year,2000)?
d1 = 335;  %when did you start assimilation, day
nd = 1;    %how many days did you run assimilation?
nom=nd*1440-1; %calculate upper time limit, minutes (used in reading gitm truth and middle files)

nom=30; %only for initial tests (30 minute simulations). REMOVE this line when you want to plot real data.

tt=awt:awt:nom;  %Time Gitm

%read RhoT from GITM binary (data) files
[LonT, LatT, AltT, RhoT] = ...
    data2mat3d_m (y2, d1, tt, truth_run_dir, pbs_file_t); 

% %read RhoC from GITM champ.dat file and interpolate it to tt to get RhoCT
[tc, LonC, LatC, AltC, RhoC, RhoCu ...
    LonCT, LatCT, AltCT, RhoCT, RhoCuT ] = ...
    dat2mat(y2, d1, tt, [dart_install_dir '/models/gitm/GITM2/run'], 'champ.dat'); %load champ data

%read RhoC from GITM champ.dat file and interpolate it to tt to get RhoCT
%[tc, LonC, LatC, AltC, RhoC, RhoCu ...
%    LonCT, LatCT, AltCT, RhoCT, RhoCuT ] = ...
%    dats2mat(y2, d1, nd, tt, [dart_install_dir '/observations/CHAMP/work']); %load champ data

%VTecTC=calc_vtec(tt, LonT, LatT, AltT, RhoT, LonCT, LatCT); %experimental

%interpolate RhoT to champ location to get RhoTC 
RhoTC=interpolateN( RhoT,{LonT;LatT;AltT;tt},{LonCT;LatCT;AltCT;tt}); %interpolate RhoT to interpolated champ location

%read f107 at tt-like times (really see tf1)
[f107T, tf1, f107Ta]=f107_p([dart_install_dir '/models/gitm/GITM2/srcData'],2000+y2,nan,d1,nd,0); %tom is the true f107 as measured on Earth on that day


%read RhoG from GITM grace.dat file and interpolate it to tt to get RhoGT
[tg, LonG, LatG, AltG, RhoG, RhoGu ...
    LonGT, LatGT, AltGT, RhoGT, RhoGuT ] = ...
    dat2mat(y2, d1, tt, [dart_install_dir '/models/gitm/GITM2/run'], 'grace.dat'); %load grace data

%interpolate RhoT to grace location to get RhoTG 
RhoTG=interpolateN( RhoT,{LonT;LatT;AltT;tt},{LonGT;LatGT;AltGT;tt}); %interpolate RhoT to interpolated champ location

% Interpolate RhoT to subsolar point (SP) 
LonTS=mod( (1440-mod(tt,1440))/1440*360+180 , 360); %subsolar point (point on Earth closest to the Sun) moves Westward at 360deg/day=15deg/hr. It is where the local noon is (ie where it is 12:00 military time, not 00:00).
LatTS=  (23+26/60)*sin(2*pi/365*((d1+tt/1440)-78) ); %subsolar point moves between the Tropic of Cancer (Northern Solstice, let's say Jun 22) and the Tropic of Capricorn (Southern Solstice, let's say Dec 21). 
%+ so it is a sinusoid with magnitude of +-23.43deg (that's where the Tropics are) and frequency of 2pi rev/365 days. The phase shift is approximately datenum(2002,3,20)-datenum(2002,1,1) 
AltTS=393983.5 *ones(size(tt)); %about 400km, one of the gridpoints in GITM in this particular setup
RhoTS=interpolateN( RhoT(:,:,:,:), ... %doesn't work too well if you're interpolating to the gridpoint (it takes the cube closest to the gridpoint, if you want perfection at gridpoints, use interpolateN_b - bilinear interp)
    {LonT;LatT;AltT;tt}, ...
    {LonTS;LatTS;AltTS;tt}); 

%VTecTS=calc_vtec(tt, LonT, LatT, AltT, RhoT, LonTS, LatTS); %experimental

[~,wlon]=min(abs(LonT-(96.2626+180) )) %lon closest to AA (WI in 9x9)
[~,wlat]=min(abs(LatT-42.2526)) %lat closest to AA (WI in 9x9)
walt=35;

LonTA=LonT(wlon)+0*tt; %interpolating Truth to gridpoint closest to Ann arbor really comes down to picking one value and making it an array (ie 42+[0 0 0 0...] = [42 42 42 42...] ), because that gridpoint doesn't move over time
LatTA=LatT(wlat)+0*tt;
AltTA=AltT(walt)+0*tt;
RhoTA=squeeze(RhoT(wlon,wlat,walt,:))'; %this turns out to be picking a 1D array from a big 4D array

%VTecTA=calc_vtec_grid(AltT, RhoT, wlon, wlat); %experimental


%% Write a DART input file for SP example (example 1 in JASTP paper, co0 in PBS naming convention)

%IMPORTANT: took me couple months to catch this: DART input file (obs_seq.out) should contain
%data at reasonable time resolution. What is reasonable? Let's compare to
%CHAMP real data. CHAMP real data comes about every 47 seconds. So if you
%have truth data every (one) minute, it is reasonable, but DON'T make the
%mistake I made - I recorded truth data every 30 minutes and expected DART
%to work more or less the same for both real and simulated data - they
%didn't - the simulated data was way too sparse compared to real and it had convergence issues. 
%
%Moral of the story: you should have truth data about every minute, so tt
%should = 1:1:nom, ta=1:1:nom, (tt can= 1:1:nom, but doesn't have to?) 
%You can still do assimilation every 30 minutes, but DART will use all that
%data at one instance (it doesn't do time interpolation, so it assumes
%whatever was in this assimilation window, happened at this assimilation
%step.

i = find( tt>0 & tt<=nom ); %at what times do you want density in this file to be written out? >0, <nom means find first nonzero time and last time that you want assimilation to run till (controlled by variable called nd)
%mat2dart(y2, d1, tt(i), [dart_install_dir '/observations/text_GITM/work'], ...
%    [truth_run_dir '/obs_seq_01to03_SP_lin_2p6_2d.out'], ...
%    LonTS(i), LatTS(i), AltTS(i), RhoTS(i), kind_sd, kind); 

%% Write a GITM input file for SP 

i = find( tt>0 & tt<=nom ); %at what times do you want density in this file to be written out? >0, <nom means find first nonzero time and last time that you want assimilation to run till (controlled by variable called nd)
%mat2dat(y2, d1, tt(i), [truth_run_dir '/subp.dat'], ...
%    LonTS(i), LatTS(i), AltTS(i), RhoTS(i), kind_sd); 


%% Write a DART input file for Ann Arbor example (example 2? in JASTP paper, co2 in PBS naming convention)

i = find( tt>0 & tt<=nom ); %at what times do you want density in this file to be written out? >0, <nom means find first nonzero time and last time that you want assimilation to run till (controlled by variable called nd)
%mat2dart(y2, d1, tt(i), [dart_install_dir '/observations/text_GITM/work'], ...
%    [truth_run_dir '/obs_seq_01to03_AA_lin_2p6_2d.out'], ...
%    LonTA(i), LatTA(i), AltTA(i), RhoTA(i), kind_sd, kind);


%% Write a GITM input file for Ann Arbor 

i = find( tt>0 & tt<=nom ); %at what times do you want density in this file to be written out? >0, <nom means find first nonzero time and last time that you want assimilation to run till (controlled by variable called nd)
%mat2dat(y2, d1, tt(i), [truth_run_dir '/anna.dat'], ...
%    LonTA(i), LatTA(i), AltTA(i), RhoTA(i), kind_sd);

%% Write a DART input file for Champ Simulated example (example 3? in JASTP paper, co4 in PBS naming convention)

i = find( tt>0 & tt<=nom ); %at what times do you want density in this file to be written out? >0, <nom means find first nonzero time and last time that you want assimilation to run till (controlled by variable called nd)
%mat2dart(y2, d1, tt(i), [dart_install_dir '/observations/text_GITM/work'], ...
%    [truth_run_dir '/obs_seq_01to03_CS_lin_2p6_2d.out'], ...
%    LonTC(i), LatTC(i), AltTC(i), RhoTC(i), kind_sd, kind);

%for comparison plot, see next section

%% Write a DART input file for Champ Real example (example 4? and5 in JASTP paper, co6 and co8 in PBS naming convention)

i = find( tc>0 & tc<=nom );
%mat2dart(y2, d1, tc(i), [dart_install_dir '/observations/text_GITM/work'], ...
%    [truth_run_dir '/obs_seq_01to05_CR_lin_3x_2d.out'], ...
%    LonC(i), LatC(i), AltC(i), RhoC(i), 3*RhoCu(i), kind)

% plot(tt,RhoTC,'b',t,r,'r',tc(i),RhoC(i),'g') % compare the normal and interpolated RhoS

%% Load GITM middle data from 3DALL

[~,awm] = unix(['grep "awt=" < ' middle_run_dir '/' pbs_file_m ' | awk ''{print $1}''']); i=find(awm=='='); awm=str2double( awm((i+1):(end-1)));    %what assimilation window did you use, in seconds, pbs_*.sh:line 23
awm=awm/60; %convert aw to minutes
tm=awm:awm:nom;  %Time Middle (I hope you ran T and M for the same amound of time, so I don't have different nom's here and there, but i will never know...

%load binary files
[~,~,~, RhoM] = ... %hopefully lons and lats are the same in M as in T, please don't run truth and middle at different resolutions as I assume LonM=LonT!
    data2mat3d_m (y2, d1, tm, middle_run_dir, pbs_file_m); %read the actual values
% cd(middle_run_dir);load('rhom_half');

%read RhoC from GITM champ.dat file and interpolate it to tm (which might be different from tt) to get RhoCM
[~, ~, ~, ~, ~, ~, ... %LonC and other *C were already read in truth_run
    LonCM, LatCM, AltCM, RhoCM, RhoCuM ] = ...
    dat2mat(y2, d1, tm, [dart_install_dir '/models/gitm/GITM2/run'], 'champ.dat'); %load champ data

%interpolate RhoM to champ location to get RhoMC 
RhoMC=interpolateN( RhoM,{LonT;LatT;AltT;tm},{LonCM;LatCM;AltCM;tm}); %interpolate RhoM to interpolated champ location

%read RhoG from GITM grace.dat file and interpolate it to tt to get RhoGT
[~, ~, ~, ~, ~, ~, ... %LonG and other *G were already read in truth_run
    LonGM, LatGM, AltGM, RhoGM, RhoGuM ] = ...
    dat2mat(y2, d1, tm, [dart_install_dir '/models/gitm/GITM2/run'], 'grace.dat'); %load grace data

%interpolate RhoM to grace location to get RhoMG 
RhoMG=interpolateN( RhoM,{LonT;LatT;AltT;tm},{LonGM;LatGM;AltGM;tm}); %interpolate RhoM to interpolated champ location

%interpolate RhoM to subsolar point location (assume tm=tt, LonM=LonT...)
RhoMS=interpolateN( RhoM(:,:,:,:), ... %doesn't work too well if you're interpolating to the gridpoint (it takes the cube closest to the gridpoint, if you want perfection at gridpoints, use interpolateN_b - bilinear interp)
    {LonT;LatT;AltT;tt}, ...
    {LonTS;LatTS;AltTS;tt}); 
%interpolate RhoM to ann arbor (assume tm=tt, LonM=LonT...)
RhoMA=squeeze(RhoM(wlon,wlat,walt,:))'; %this turns out to be picking a 1D array from a big 4D array



%% Load DART Data from NC files

cd(dart_run_dir)
[~,noe] = unix(['grep "noe=" < ' pbs_file_d ' | awk ''{print $1}''']); i=find(noe=='='); noe=str2double(noe((i+1):(end-1)));    %Number Of Ensembles (ensemble members to be precise); be weary of punctuation - str2double doesn't like semicolon

[~,cfh] = unix([ 'grep "cfh=" < ' pbs_file_d ' | awk ''{print $1}''' ]); i=find(cfh=='='); cfh=str2double( cfh((i+1):(end-1)) );   %horizontal cutoff (radians) also in input.nml:line 83
cfh=cfh*180/pi; %convert cfh to degrees
[~,cfv] = unix([ 'grep "cfv=" < ' pbs_file_d ' | awk ''{print $1}''' ]); i=find(cfv=='='); cfv=str2double( cfv((i+1):(end-1)) );   %vertical cutoff (meters) also in input.nml:line 254

[~,aw] = unix(['grep "aw=" < ' pbs_file_d ' | awk ''{print $1}''']); i=find(aw=='='); aw=str2double( aw((i+1):(end-1)));    %what assimilation window did you use, in seconds, pbs_*.sh:line 23
aw=aw/60; %convert aw to minutes


%load PRIOR
ncid  = netcdf.open('preassim.nc','NOWRITE');
LonD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'LON') ) );
LatD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'LAT') ) );
AltD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'ALT') ) );
% TempEr= netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Temperature') ); %Temperature is  not read right now (reading less variables requires less memory)
% NdsEr = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'iO_3P_NDensityS') );
% IdsEr = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'ie_IDensityS') );
f107Dr= netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'f107') );
RhoDr = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Rho') );
netcdf.close(ncid)

%load POSTERIOR
ncid  = netcdf.open('analysis.nc','NOWRITE');
td =  netcdf.getVar(ncid, netcdf.inqVarID(ncid,'time'))... %td=Time DART (number of minutes since beginning of 12/01/2002)
    + datenum([1601 01 01 00 00 00]) ... %DART starts its timekeeping from 01/01/1601
    - datenum([2002 12 01 00 00 00]);    %This file counts time       from 12/01/2002
td=td*24*60; %convert DART time to minutes
max(td)      %display maximum of DART time as a safety check
% LonD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'LON') ) ); %LonD and such should be the same as un Prior, so no need to read again
% LatD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'LAT') ) );
% AltD   = double(netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'ALT') ) );
% TempEo= netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Temperature') );
% NdsEo = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'iO_3P_NDensityS') );
% IdsEo = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'ie_IDensityS') );
f107Do= netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'f107') );
RhoDo = double( netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Rho') ) );
netcdf.close(ncid)

%% Interp4 C->D->C
LonCD=interp1np(tc,LonC,td);
LatCD=interp1(tc,LatC,td);
AltCD=interp1(tc,AltC,td);
RhoCD=interp1(tc,RhoC,td);

RhoDrmC=interpolateN(squeeze( RhoDr(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonCD;LatCD;AltCD;td}); %m - mean, b - backwards, i - interpolated
RhoDomC=interpolateN(squeeze( RhoDo(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonCD;LatCD;AltCD;td});
RhoDrsC=interpolateN(squeeze( RhoDr(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonCD;LatCD;AltCD;td});
RhoDosC=interpolateN(squeeze( RhoDo(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonCD;LatCD;AltCD;td});

LonGD=interp1np(tg,LonG,td);
LatGD=interp1(tg,LatG,td);
AltGD=interp1(tg,AltG,td);
RhoGD=interp1(tg,RhoG,td);

RhoDrmG=interpolateN(squeeze( RhoDr(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonGD;LatGD;AltGD;td}); %m - mean, b - backwards, i - interpolated
RhoDomG=interpolateN(squeeze( RhoDo(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonGD;LatGD;AltGD;td});
RhoDrsG=interpolateN(squeeze( RhoDr(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonGD;LatGD;AltGD;td});
RhoDosG=interpolateN(squeeze( RhoDo(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonGD;LatGD;AltGD;td});

LonSD=interp1nr(tt,LonTS,td);
LatSD=interp1(tt,LatTS,td,'linear','extrap');
AltSD=interp1(tt,AltTS,td,'linear','extrap');
RhoSD=interp1(tt,RhoTS,td,'linear','extrap');

RhoDrmS=interpolateN(squeeze( RhoDr(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonSD;LatSD;AltSD;td}); %m - mean, b - backwards, i - interpolated, sp -SubsolarPoint
RhoDomS=interpolateN(squeeze( RhoDo(:,:,:,1,:) ),{LonD;LatD;AltD;td},{LonSD;LatSD;AltSD;td});
RhoDrsS=interpolateN(squeeze( RhoDr(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonSD;LatSD;AltSD;td});
RhoDosS=interpolateN(squeeze( RhoDo(:,:,:,2,:) ),{LonD;LatD;AltD;td},{LonSD;LatSD;AltSD;td});


%% Check if data makes sense


lona=0*tc; %preallocate
lata=0*tc;
alta=0*tc;
for i=1:length(tc)
    [~, lona(i)] = min(abs(LonD(:) - LonC(i)));
    [~, lata(i)] = min(abs(LatD(:) - LatC(i)));
    [~, alta(i)] = min(abs(AltD(:) - AltC(i)));
end

%display mins and maxs
disp(' ')
disp(['longitude indeces champ traveled through [ ' num2str([min(lona) max(lona)]) ' ]'])
disp(['latitude  indeces champ traveled through [ ' num2str([min(lata) max(lata)]) ' ]'])
disp(['altitude  indeces champ traveled through [ ' num2str([min(alta) max(alta)]) ' ]'])
cd(dart_run_dir)     %come back to parent directory



%display more checks
disp(' ')
disp('check that units (orders of magntude) are the same Column-wise and values are the same on the first 2 lines')
disp([td(1) LonD(wlon) LatD(wlat) AltD(walt)/1000]) %check that units (orders of magntude) are the same Column-wise and values are the same on the first 2 lines
% disp([LonT(wlon)     LatT(wlat) AltT(walt)/1000])
disp([tt(1) LonCT(1)   LatCT(1)   AltCT(1)/1000])
disp([tc(1) LonC(1) LatC(1) AltC(1)/1000]) %tc(5515) is closest to tt(0) for my set of truth data - might differ from yours

cd(dart_run_dir)




%% switch to hours
tt=tt/60;
tm=tm/60;
td=td/60;
tc=tc/60;
tg=tg/60;
tf1=tf1/60;


%% figure 7,8,9 in gitm/GITM2/srcDoc/thermo.pdf (for other figures, search this file for "thermo.pdf")
plot_cr6 %currently segfaults for some reason - no clue why
		       disp('past cr6')
%% CHAMP traj map

figure(5); clf
set(gcf,'position',[150 150 500 400]);
load('topo.mat','topo','topomap1')
contour(0:359,-90:89,topo,[0 0],'k')


hold on


cma=7*10^-12;
cmi=2*10^-12;



%%%% grid cells boundaries
% for i=1:(length(LonD)-1)
%     line((LonD(i)+LonD(i+1))/2*[1 1], [-90 90],'color',[.5 .5 .5])
% end
% for i=1:(length(LatD)-1)
%     line([0 360], (LatD(i)+LatD(i+1))/2*[1 1],'color',[.5 .5 .5])
% end

%%% grid cell centers 
% s=4; %MarkerSize
% [x,y]=meshgrid(LonD,LatD);
% qt=RhoT(:,:,walt,1)'; %transpose is needed only for
% qt=reshape(qt,prod(size(qt)),1); %columns (values at fixed lon, varying lat) are put into one big column (column-wise first, =qt(:), but that looks more cryptic) 
% % cma=max(qt(:));
% % cmi=min(qt(:));
% % cma=10*10^-12;
% % cmi=1*10^-12;
% set(gca,'clim',[cmi cma])
% ylabel(colorbar,'Mass Density [kg m^{-3}]');
% ca=jet(100);
% cc=ca(round( (qt-cmi)./(cma-cmi)*99 )+1,:);
% 
% % scatter(x(:),y(:),s,cc,'filled');
% for i=1:(length(LonD))
% for j=1:(length(LatD))
%     plot(LonD(i),LatD(j),'o', ...
%          'MarkerSize',s, ...
%         'MarkerFaceColor',cc(j+(i-1)*length(LatD),:), ...
%         'MarkerEdgeColor','none')
%     
% end
% end

%%% champ sample data 
t=(1:240)+80; %about two CHAMP's orbits, starting after some time 5515
s=6; %MarkerSize
% cma=max(RhoC(t) );
% cmi=min(RhoC(t) );
% cma=10*10^-12;
% cmi=1*10^-12;
set(gca,'clim',[cmi cma])
ylabel(colorbar,'Mass Density [kg m^{-3}]');
ca=jet(100);
cc=ca(max(min( round( (RhoC(t)-cmi)./(cma-cmi)*99 )+1, 100), 1),:);

h=text(LonC(t(1)),LatC(t(1))+5,'Start','color','b','fontsize',12,'fontweight','b');
text(LonC(t(end))-18,LatC(t(end))-5,'Finish','color','b','fontsize',12,'fontweight','b')
% scatter(LonC(t),LatC(t),30,cc,'filled')
for i=1:(length(t)-1)
    if LonC(t(i))>280 && LonC(t(i+1))<80
        plot([LonC(t(i)) LonC(t(i+1))+360],LatC(t([i i+1])),'-o', ...
            'color',cc(i,:), 'MarkerSize',s, ...
            'MarkerFaceColor',cc(i,:), ...
            'MarkerEdgeColor','none')
        plot([LonC(t(i))-360 LonC(t(i+1))],LatC(t([i i+1])),'-o', ...
            'color',cc(i,:), 'MarkerSize',s, ...
            'MarkerFaceColor',cc(i,:), ...
            'MarkerEdgeColor','none')
    else
        plot(LonC(t([i i+1])),LatC(t([i i+1])),'-o', ...
            'color',cc(i,:), 'MarkerSize',s, ...
            'MarkerFaceColor',cc(i,:), ...
            'MarkerEdgeColor','none')
    end
end




hold off

xlabel('Longitude [deg]');
set(gca,'XTick',0:60:360)
ylabel('Latitude [deg]');
set(gca,'YTick',-90:30:90)

xlim([0 360]);
ylim([-90 90]);
print(gcf,'-dpdf','champ_traj')


%% Figure that used to be in thermo.pdf paper (maps with mass density on top)

figure(5); clf
set(gcf,'position',[150 150 500 400]);
load('topo.mat','topo','topomap1')
contour(0:359,-90:89,topo,[0 0],'k')


hold on


cma=7*10^-12;
cmi=1*10^-12;



%%%% grid cells boundaries
% for i=1:(length(LonD)-1)
%     line((LonD(i)+LonD(i+1))/2*[1 1], [-90 90],'color',[.5 .5 .5])
% end
% for i=1:(length(LatD)-1)
%     line([0 360], (LatD(i)+LatD(i+1))/2*[1 1],'color',[.5 .5 .5])
% end

%%% grid cell centers 
s=4; %MarkerSize
[x,y]=meshgrid(LonD,LatD);
qt=RhoT(:,:,walt,1)'; %transpose is needed only for
qt=reshape(qt,prod(size(qt)),1); %columns (values at fixed lon, varying lat) are put into one big column (column-wise first, =qt(:), but that looks more cryptic) 
% cma=max(qt(:));
% cmi=min(qt(:));
% cma=10*10^-12;
% cmi=1*10^-12;
set(gca,'clim',[cmi cma])
ylabel(colorbar,'Mass Density [kg m^{-3}]');
ca=jet(100);
cc=ca(max(min( round( (qt-cmi)./(cma-cmi)*99 )+1, 100), 1),:);

% scatter(x(:),y(:),s,cc,'filled');
for i=1:(length(LonD))
for j=1:(length(LatD))
    plot(LonD(i),LatD(j),'o', ...
         'MarkerSize',s, ...
        'MarkerFaceColor',cc(j+(i-1)*length(LatD),:), ...
        'MarkerEdgeColor','none')
    
end
end

%%% champ sample data 
% t=5515+(1:240)+80; %about two CHAMP's orbits, starting after some time
% s=6; %MarkerSize
% % cma=max(RhoC(t) );
% % cmi=min(RhoC(t) );
% % cma=10*10^-12;
% % cmi=1*10^-12;
% set(gca,'clim',[cmi cma])
% ylabel(colorbar,'Mass Density [kg m^{-3}]');
% ca=jet(100);
% cc=ca(round( (RhoC(t)-cmi)./(cma-cmi)*99 )+1,:);
% 
% h=text(LonC(t(1)),LatC(t(1))+5,'Start','color','b','fontsize',12,'fontweight','b');
% text(LonC(t(end))-18,LatC(t(end))-5,'Finish','color','b','fontsize',12,'fontweight','b')
% % scatter(LonC(t),LatC(t),30,cc,'filled')
% for i=1:(length(t)-1)
%     if LonC(t(i))>280 && LonC(t(i+1))<80
%         plot([LonC(t(i)) LonC(t(i+1))+360],LatC(t([i i+1])),'-o', ...
%             'color',cc(i,:), 'MarkerSize',s, ...
%             'MarkerFaceColor',cc(i,:), ...
%             'MarkerEdgeColor','none')
%         plot([LonC(t(i))-360 LonC(t(i+1))],LatC(t([i i+1])),'-o', ...
%             'color',cc(i,:), 'MarkerSize',s, ...
%             'MarkerFaceColor',cc(i,:), ...
%             'MarkerEdgeColor','none')
%     else
%         plot(LonC(t([i i+1])),LatC(t([i i+1])),'-o', ...
%             'color',cc(i,:), 'MarkerSize',s, ...
%             'MarkerFaceColor',cc(i,:), ...
%             'MarkerEdgeColor','none')
%     end
% end




hold off

xlabel('(a)  Longitude [deg]');
set(gca,'XTick',0:60:360)
ylabel('Latitude [deg]');
set(gca,'YTick',-90:30:90)

xlim([0 360]);
ylim([-90 90]);

print(gcf,'-dpdf','gitm_h')

%% Figure that used to be in thermo.pdf paper (GITM's altitude grid)

figure(5); clf

set(gcf,'position',[150 150 500 400]);


s=4; %MarkerSize
[x,y]=meshgrid(LonT,AltT/1000);
qt=squeeze(log10(RhoT(:,1,:,1)))'; %transpose is needed only for
qt=reshape(qt,prod(size(qt)),1); %columns (values at fixed lon, varying lat) are put into one big column (column-wise first, =qt(:), but that looks more cryptic) 
% cma=max(qt(:));
% cmi=min(qt(:));
cma=-6;
cmi=-14;
set(gca,'clim',[cmi cma])
ylabel(colorbar,'Log_{10} (Mass Density [kg m^{-3}])');
ca=jet(100);
cc=ca(max(min( round( (qt-cmi)./(cma-cmi)*99 )+1, 100), 1),:);
hold on
% scatter(x(:),y(:),s,cc,'filled');
for i=1:(length(LonT))
for j=1:(length(AltT))
    plot(LonT(i),AltT(j)/1000,'o', ...
         'MarkerSize',s, ...
        'MarkerFaceColor',cc(j+(i-1)*length(AltT),:), ...
        'MarkerEdgeColor','none')
    
end
end
hold off



% set(gca,'XTick',LonT)
a=[AltT(1) AltT(13:4:end)']/1000;
set(gca,'YTick',a)
set(gca,'YTickLabel',round(a))
set(gca,'YDir','normal')
set(gca,'clim',[cmi cma])

xlim([0 360])
ylim([min(AltT) max(AltT)]/1000)
xlabel('(b)  Longitude [deg]')
set(gca,'XTick',0:60:360)
ylabel('Altitude (km)')

box on
print(gcf,'-dpdf','gitm_v')

%% Figure 1b in thermo.pdf %for 1a see sun_rot.m
figure(1)
set(gcf,'position',[521   242   420   273]);
% set(gcf,'position',[150 150 500 400]);

% subplota(1,2,1,2,.1,0)
% subplot(1,2,2)

time_plot(tt,RhoTC,'GITM output', ...
    [],[],'GITM without EAKF', ...
    tc,RhoC,'CHAMP measurement', ...
    tc,RhoCu,'CHAMP measurement +/-SD', ...
    [],[],50,[],'CHAMP measurement localized', ...
    [],[],[],'EAKF ensemble mean', ...
    [],[],[],'EAKF ensemble mean +/- SD', ...
    [],[],[], ...
    [],[4 1],-1.5:.01:3,'na', ...
    'NE','none',[],0.5,10^12)
box on

ylim([0.9 8.5])
xlim([1 4.2])
ylabel('Mass density  x10^{-12} [kg m^{-3}]')
xlabel('(b)  Hours since 00UT 01/12/2002')

set(gca,'XTick',0:0.5:48)


print(gcf,'-depsc','p_diff')


%% Figure 6 in thermo.pdf: Map Gaspari Cohn demonstrated via lon-lat

figure(9); clf
set(gcf,'position',[150 150 500 250]);
set(gca,'position',[.11 .15 0.8 0.8])
% subplota(1,1,1,1,0,0);
load('topo.mat','topo','topomap1')

[n,t]=meshgrid(LonD,LatD);
[tts,ns,tss,as,rs,gs]=loc2(nan, LonD(wlon), LatD(wlat), AltD(walt), ...
    0*n(:), n(:), t(:), AltD(walt)+0*n(:), 0*n(:), ...
    cfh, cfv);



% h=scatter(ns,tss,20,[ones(size(gs)) 1-gs 1-gs ],'filled'); %color
h=scatter(ns,tss,20,[1-gs 1-gs 1-gs ],'filled'); %bw

hold on
contour(0:359,-89:90,topo,[0 0],'k')
hold off

xlabel('Longitude [deg]');
set(gca,'XTick',0:60:360)
ylabel('Latitude [deg]');
set(gca,'YTick',-90:30:90)

xlim([0 360]);
ylim([-90 90]);
box on

 set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 1.825 8 4]);  
 print(gcf,'-depsc',['loc_h'])

%% Map Gaspari Cohn demonstrated via lon-alt

figure(10); clf
set(gcf,'position',[150 150 500 250]);
set(gca,'position',[.11 .15 0.8 0.8])
% subplota(1,1,1,1,0,0);
load('topo.mat','topo','topomap1')

[n,a]=meshgrid(LonD,AltD);
[tts,ns,tss,as,rs,gs]=loc2(nan, LonD(wlon), LatD(wlat), AltD(walt), ...
    0*n(:), n(:), LatD(wlat)+0*n(:), a(:), 0*n(:), ...
    cfh, cfv);



h=scatter(ns,as/1000,20,[ones(size(gs)) 1-gs 1-gs ],'filled');

% hold on
% contour(0:359,-89:90,topo,[0 0],'k')
% hold off

xlabel('Longitude [deg]');
set(gca,'XTick',0:60:360)
ylabel('Altitude (km)');
set(gca,'YTick',AltD(1:4:end)/1000)
set(gca,'YTickLabel',round(AltD(1:4:end)/1000))

% set(gca,'YTick',-90:30:90)

ylim([min(AltD) max(AltD)]/1000)
xlim([0 360]);
box on
grid on

set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 1.825 8 4]);  
print(gcf,'-depsc',['loc_v'])


%% Map Gaspari Cohn demonstrated via lon-GC

figure(11); clf
set(gcf,'position',[150 150 500 250]);
set(gca,'position',[.11 .15 0.8 0.8])
% subplota(1,1,1,1,0,0);
load('topo.mat','topo','topomap1')

n=1:.1:359;
a=AltD(walt)+0*n;
[tts,ns,tss,as,rs,gs]=loc2(nan, LonD(wlon), LatD(wlat), AltD(walt), ...
    0*n(:), n(:), LatD(wlat)+0*n(:), a(:), 0*n(:), ...
    cfh, cfv);



h=scatter(ns,gs,20,[ones(size(gs)) 1-gs 1-gs ],'filled');

% hold on
% contour(0:359,-89:90,topo,[0 0],'k')
% hold off

xlabel('Longitude [deg]');
set(gca,'XTick',0:60:360)
ylabel('Localization Function ');
% set(gca,'YTick',-90:30:90)

 ylim([-.1 1.1])
xlim([0 360]);
box on
grid on
set(gcf,'PaperType', 'A5');  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize'))); set(gcf,'PaperPosition',[.125 1.825 8 4]);  
print(gcf,'-depsc',['loc_f'])

%% Map Gaspari Cohn demonstrated via ground tracks

figure(9); clf

% qt = RhoT(:,:,alt2,:);
% qr = RhoDr(:,:,alt2,1,:);
% yl =[  min( [min(qt(:)) min(qr(:))] ) max( [max(qt(:)) max(qr(:))] )  ];


subplota(1,1,1,1,0,0);
load('topo.mat','topo','topomap1')


dn=LonD(lon1)-LonC;
dt=LatD(lat1)-LatC;
da=AltD(walt)-AltC;
gc=GC(cfh,dn);
gc=GC(cfh,dn).*GC(cfh,dt);
gc=GC(cfh,dn).*GC(cfh,dt).*GC(cfv,da);

[gci,j] = sort(gc,1,'ascend');
tci=tc(j);
lnci=LonC(j);
ltci=LatC(j);
laci=AltC(j);
rci=RhoC(j);
tti=gci>.01; %test
gs=gci(tti); %select obs
tts=tci(tti); %select time
lns=lnci(tti); %select rho
lts=ltci(tti); %select rho
las=laci(tti); %select rho
rs=rci(tti); %select rho

h=scatter(lns,lts,50,[ones(size(gs)) 1-gs 1-gs ],'filled');

hold on
contour(0:359,-89:90,topo,[0 0],'k')
hold off

xlim([0 360])
ylim([-90 90])
set(gca,'XTick',[])
set(gca,'YTick',[])


%% obs_diag

al = plot_rmse_xxx_evolution_ALEX('obs_diag_output.nc','spread');
ale=al{5};
te=ale(:,1);
te=(te-te(1))*24*60;

close
close
a=findall(gcf,'type','axes');
set(a(end),'YScale','log');
print(gcf,'-dpng','obs_diag')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
