function Rho=mat2dat(y2, d1, t, dat_file, Lon, Lat, Alt, Rho, RhoU)

%mat2dat - no interpolation or extrapolation done here - just write it as is
%t is the time at which you have Lon, Lat, Alt, Rho and RhoU available (so
%length(t) must= length(Rho) = length(Lon) etc. RhoU can be scalar, in
%which case it is just made into array (RhoU*ones(length(1),1))

% DART $Id$
% CREDIT: Alexey Morozov

sp=ones(length(t),1);

if isscalar(RhoU) %if RhoU is given as 1 number (scalar), assume it's constant for the whole time
    RhoU=RhoU*sp;
end

[yy mo dd hh mm ss]=datevec(datenum(2000+y2,0,0)+d1+t(:)'/1440); %t(:)' makes it a row. Hopefully you didn't do something funny like give me a matrix for t and I accidentally flattened it!

%% write 
idt=fopen(dat_file,'w'); 
fwrite(idt,[' '  10]);

fwrite(idt,['File made on : ' datestr(now) ' ' 10]);
fwrite(idt,['Matlab code : mat2dat.m, pwd: ' pwd  10]);

fwrite(idt,[' '  10]);
fwrite(idt,['#START '  10]);

text_M=[num2str(yy') 32*sp ...
    num2str(mo','%02.0f') 32*sp ...
    num2str(dd','%02.0f') 32*sp ...
    num2str(hh','%02.0f') 32*sp ...
    num2str(mm','%02.0f') 32*sp ...
    num2str(ss','%02.0f') 32*sp ...
    '0'*sp 32*sp ...
    num2str(Lon(:),'%10.5f') 32*sp ...
    num2str(Lat(:),'%9.5f') 32*sp ...   %what (:) does, is make something (matrix, or row, or column) into a column. So no mater what you gave me, I'll make it into a column
    num2str(Alt(:)/1000,'%7.3f') 32*sp ...
    num2str(Rho(:),'%12.6e') 32*sp ...
    num2str(RhoU(:)) 32*sp ...
    10*sp];

fwrite(idt,text_M');

disp(['just so you know, here is what the end of the ' dat_file ' file looks like'])
unix(['tail ' dat_file]);
disp('and that''s how many lines it has')
unix(['wc -l < ' dat_file]);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
