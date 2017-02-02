function [LongitudeT,LatitudeT,AltitudeT,RhoT]=data2mat3d_m(y2, d1, tg, directory, pbs_file)
% GITM files (3DALL_*.b000*) reader

% DART $Id$
% CREDIT: Alexey Morozov


p=pwd; %save current path
cd(directory) %go to where the pbs_file for this run is

% GITM parameters (read automatically from pbs_file
[~,nLons] = unix(['grep "ncn=" < ' pbs_file ' | awk ''{print $1}''']); i=find(nLons=='='); nLons=str2double( nLons((i+1):(end-1)) );  % number of longitudes PER block (from ModSize.f90)
[~,nLats] = unix(['grep "nct=" < ' pbs_file ' | awk ''{print $1}''']); i=find(nLats=='='); nLats=str2double( nLats((i+1):(end-1)) ); % number of latitudes  PER block
[~,nAlts] = unix(['grep "nca=" < ' pbs_file ' | awk ''{print $1}''']); i=find(nAlts=='='); nAlts=str2double( nAlts((i+1):(end-1)) ); % number of altitudes  total (since there is no vertical blocking)
[~,nBlocksLon] = unix(['grep "nbn=" < ' pbs_file ' | awk ''{print $1}''']); i=find(nBlocksLon=='='); nBlocksLon=str2double( nBlocksLon((i+1):(end-1)) ); % number of blocks longitude-wise (from UAM.in)
[~,nBlocksLat] = unix(['grep "nbt=" < ' pbs_file ' | awk ''{print $1}''']); i=find(nBlocksLat=='='); nBlocksLat=str2double( nBlocksLat((i+1):(end-1)) ); % number of blocks latitude -wise

cd([directory '/data']) %go to where the data itself is (usually just in 'data' folder here)

%preallocate
a=zeros( (nLons*nBlocksLon+4)*(nLats*nBlocksLat+4)*(nAlts+4), 39, length(tg)); %preallocate memory for giant read-in array

for k=1:length(tg) %iterate over number of minutes
  [yy mo dd hh mm ss]=datevec(datenum(2000+y2,0,0)+d1+tg(k)/1440);
    
    
    fname=['3DALL_t' num2str(mod(yy,2000),'%02.0f') num2str(mo,'%02.0f') num2str(dd,'%02.0f') ...
			     '_' num2str(hh,'%02.0f') num2str(mm,'%02.0f') num2str(ss,'%02.0f') ...
					 '.bin'];
			     
			     disp(fname) %enable for debug
			     id=fopen(fname,'r');

    nvars=fread(id,1,'int32')/8; %number of scalar variables  to be read (var at that time, at that location is a single number I call "SCALAR")
    version=fread(id,nvars,'float64'); 
    nvars=fread(id,1,'int32')/8
        
    nvars=fread(id,1,'int32')/4; 
    dims=fread(id,nvars,'int32'); 
    nvars=fread(id,1,'int32')/4

    nvars=fread(id,1,'int32')/4; 
    NVARS=fread(id,nvars,'int32'); 
    nvars=fread(id,1,'int32')/4

for iVar=1:NVARS
    nvars=fread(id,1,'int32') 
    var_names(:,iVar)=fread(id,nvars,'char');            %write(iOutputUnit_) Variables(iVar)
    nvars=fread(id,1,'int32')
end

var_names=var_names'; %transpose the char array so that string is in row

    nvars=fread(id,1,'int32')/4; 
    time=fread(id,nvars,'int32'); 
    nvars=fread(id,1,'int32')/4

for iVar=1:NVARS
    nvars=fread(id,1,'int32')/8; 
    a(:,iVar,k)=fread(id,nvars,'float64');            %write(iOutputUnit_) Variables(iVar)
    nvars=fread(id,1,'int32')/8
end

        fclose(id);        
        
        

end

%% 2. Unpacking Part

%The names are hardcoded from my copy of some 3DALL_... .header file
%really should use var_names instead, but that has exclamation marks and brackets
varnames={'Longitude'
    'Latitude'
    'Altitude'
    'Rho'
    'OU3NP'
    'OD2N'
    'ND2N'
    'NU4NS'
    'NO'
    'NU2ND'
    'NU2NP'
    'H'
    'He'
    'CO2'
    'OU1ND'
    'Temperature'
    'VDnNeast'
    'VDnNnorth'
    'VDnNup'
    'VDnNupOU3NP'
    'VDnNupOD2N'
    'VDnNupND2N'
    'VDnNupNU4NS'
    'VDnNupNO'
    'O_4SP_UpN'
    'OD2UpN'
    'ND2UpN'
    'NUpN'
    'NOUpN'
    'OU2NDUpN'
    'OU2NPUpN'
    'HUpN'
    'HeUpN'
    'e'
    'eTemperature'
    'iTemperature'
    'VDiNeast'
    'VDiNnorth'
    'VDiNup'};

for i=1:length(varnames) %really should use var_names instead, but that has exclamation marks and brackets
    eval([varnames{i} ' = reshape(a(:,i,:),(nLons*nBlocksLon+4),(nLats*nBlocksLat+4),(nAlts+4), length(tg));'])
    eval([varnames{i} 'T = ' varnames{i} '(3:end-2, 3:end-2, 3:end-2, :);'])
end
LongitudeT=squeeze(LongitudeT(:,1,1,1)); %these don't change over time and other dims
LatitudeT=squeeze(LatitudeT(1,:,1,1));
AltitudeT=squeeze(AltitudeT(1,1,:,1));
clear a

for i=1:length(varnames)
    eval(['clear ' varnames{i}]) %clean up the left over temporary variables
end
%% MINIMUMS AND MAXIMUMS

disp('MINIMA AND MAXIMA')
for i=1:length(varnames)
    disp([varnames{i} 'T ', num2str(min(eval([varnames{i} 'T(:)']))), ' ', num2str(max(eval([varnames{i} 'T(:)']))) ]);
end

%     disp(['RhoT ', num2str(min(RhoT(:))), ' ', num2str(max(RhoT(:))) ]);

cd(p) %come back to where we were before calling this function

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
