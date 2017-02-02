function [] = wrslice(file,arr,k,varargin)
%% wrslice(filename,arr,n,k) writes an array of shape [nx,ny,...]
% to a direct access binary file (float or double precisision) named
% by the string 'filename'. The file may contain multi-dimensional
% data. 'n' is the record number to be written where the record length
% is the size of the array 'arr'.
%
% eg.   wrslice('t.xyt',T,4);
%
% rdsclice(filename,A,k,type) writes an array of type 'type'.
% where type can be one of 'real*4' or 'real*8'. The default is 'real*8'.
%
% DART $Id$

% Default word-length
WORDLENGTH=8;
rtype='real*8';
ieee='b';

% Check optional arguments
args=char(varargin);
while (size(args,1) > 0)
 if deblank(args(1,:)) == 'real*4'
  WORDLENGTH=4;
  rtype='real*4';
 elseif deblank(args(1,:)) == 'real*8'
  WORDLENGTH=8;
  rtype='real*8';
 elseif deblank(args(1,:)) == 'n' | deblank(args(1,:)) == 'native'
  ieee='n';
 elseif deblank(args(1,:)) == 'l' | deblank(args(1,:)) == 'ieee-le'
  ieee='l';
 elseif deblank(args(1,:)) == 'b' | deblank(args(1,:)) == 'ieee-be'
  ieee='b';
 elseif deblank(args(1,:)) == 'c' | deblank(args(1,:)) == 'cray'
  ieee='c';
 elseif deblank(args(1,:)) == 'a' | deblank(args(1,:)) == 'ieee-le.l64'
  ieee='a';
 elseif deblank(args(1,:)) == 's' | deblank(args(1,:)) == 'ieee-be.l64'
  ieee='s';
 else
  sprintf(['Optional argument ' args(1,:) ' is unknown'])
  return
 end
 args=args(2:end,:);
end


N=size(arr);
nnn=prod(N);

[fid mess]=fopen(file,'r+',ieee);
if fid == -1
 [fid mess]=fopen(file,'w',ieee);
 if fid == -1
  sprintf('Error while opening file:\n%s',mess)
  return
 end
end
st=fseek(fid,nnn*(k-1)*WORDLENGTH,'bof');
if st ~= 0
 mess=ferror(fid);
 sprintf('There was an error while positioning the file pointer:\n%s',mess)
 return
end
count=fwrite(fid,arr,rtype);
if count ~= nnn
 sprintf('Not enough data was available to be read: off EOF?')
 return
end
st=fclose(fid);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
