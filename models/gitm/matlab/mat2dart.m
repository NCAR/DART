function mat2dart(y2, d1, t, converter_dir, obs_file, Lon, Lat, Alt, Rho, RhoU, kind)
%mat2dart(y2, d1, t, converter_dir, obs_file, Lon, Lat, Alt, Rho, RhoU, kind)
%t is the time at which you have Lon, Lat, Alt, Rho and RhoU available (so
%length(t) must= length(Rho) = length(Lon) etc. RhoU can be scalar, in
%which case it is just made into array (RhoU*ones(length(1),1))
%REQUIRES the DART/observations/text_GITM converter to be compiled!!!

% DART $Id$
% CREDIT: Alexey Morozov

sp=ones(length(t),1);

if isscalar(RhoU) %if RhoU is given as 1 number (scalar), assume it's constant for the whole time
    RhoU=RhoU*sp;
end


[yy mo dd hh mm ss]=datevec(datenum(2000+y2,0,0)+d1+t(:)'/1440); %t(:)' makes it a row. Hopefully you didn't do something funny like give me a matrix for t and I accidentally flattened it!


p=pwd; %save current path
cd(converter_dir) %go to where the txt to dart converter is

%% writing intermideary text file
idt=fopen('text.txt','w'); %intermideary text file

text_M=[num2str(kind*sp) 32*sp ... 
    num2str(Lat(:),'%9.5f') 32*sp ...   %what (:) does is make something (matrix, or row, or column) into a column. So no matter what you gave me, I'll make it into a column
    num2str(Lon(:),'%10.5f') 32*sp ...
    num2str(Alt(:),'%8.1f') 32*sp ...
    num2str(yy') 32*sp ...
    num2str(mo','%02.0f') 32*sp ...
    num2str(dd','%02.0f') 32*sp ...
    num2str(hh','%02.0f') 32*sp ...
    num2str(mm','%02.0f') 32*sp ...
    num2str(ss','%02.0f') 32*sp ...
    num2str(Rho(:),'%12.6e') 32*sp ...
    num2str(RhoU(:)) 32*sp ...
    10*sp];

fwrite(idt,text_M');

%% converting the text file into DART obs file
unix('rm obs_seq.out'); %remove old dart files (if you don't remove, new stuff will get APPENDED to the old)

disp('just so you know, here is what the end of the intermediary text file looks like')
unix('tail text.txt');
disp('and that''s how many lines it has')
unix('wc -l < text.txt');


s=unix('./text_to_obs'); %run the converter
if s~=0 
    unix('ls'); %show me what files are in this directory - LiSt the directory
    error(['Converter did not run. ' 10 'A) is this the correct directory at all?' 10 'B) Is it compiled (./quickbuild.sh) (ie Is text_to_obs present in what is shown above?)?'])
end

disp('just so you know, here is what the end of DART file looks like')
unix('tail -n 11 obs_seq.out'); %just so you know, here is what the end of DART file looks like

unix(['mv obs_seq.out ' obs_file]); %rename and possibly move (if obs_file has at least one '/' in it) the input to the dart file

cd(p)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
