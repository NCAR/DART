function MakeManyObs

% DART $Id$
% CREDIT: Alexey Morozov


obcount = 1;

%----------------------------------------------------------------------
%% Create obs exactly on heights (verttype == 3)

ALT = [100000, 101663.3, 103337.4, 105035.3, 106773.5, 108573.3, 110464.5, ...
       112492.9, 114729.3, 117279.1, 120304, 124025.5, 128512, 133760.8, ...
       139833.2, 146773.5, 154606.7, 163337.8, 172952.9, 183422, 194703.2, ...
       206747.6, 219503.3, 232918.5, 246943.8, 261533, 276642.7, 292232.6, ...
       308263.8, 324695.1, 341495.2, 358628.3, 376059.6, 393755.4, 411684.1, ...
       429816.1, 448123.9, 466582.8, 485170.7, 503868.2, 522658.4, 541526.7, ...
       560460.4, 579449, 598483.4, 617556, 636660.5, 655791.5, 674944.6, 694116.2] ;

N = length(ALT);

obcount = write_obs(obcount, 0.0, 3)
for height = 1:N
   obcount = write_obs(obcount, ALT(height), 3)
end
obcount = write_obs(obcount, 700000.0, 3)

%----------------------------------------------------------------------
%% Create lots of obs evenly spaced in heights (verttype == 3)

ALT = [0:1000:700000];

for height = 1:length(ALT)
   obcount = write_obs(obcount, ALT(height), 3)
end

%----------------------------------------------------------------------
%% Create obs on many levels (verttype == 1)

for height = 0:0.05:(N+1)
   obcount = write_obs(obcount, height, 1)
end


%----------------------------------------------------------------------
% helper function
%----------------------------------------------------------------------

function obcount = write_obs(iobs, vert, verttype )

fname = sprintf('obs_seq.%04d.in',iobs);

fid = fopen(fname,'wt');
fprintf(fid,' obs_sequence\n');
fprintf(fid,'obs_kind_definitions\n');
fprintf(fid,'           1\n');
fprintf(fid,'          38 GND_GPS_VTEC \n');
fprintf(fid,'  num_copies:            0  num_qc:            0\n');
fprintf(fid,'  num_obs:            1  max_num_obs:            1\n');
fprintf(fid,'  first:            1  last:            1\n');
fprintf(fid,' OBS            1\n');
fprintf(fid,'          -1          -1          -1\n');
fprintf(fid,'obdef\n');
fprintf(fid,'loc3d\n');
fprintf(fid,sprintf('     2.321287905152458        -1.431169986635350         %f      %d\n', vert, verttype));
fprintf(fid,'kind\n');
fprintf(fid,'          38\n');
fprintf(fid,'    90     146796\n');
fprintf(fid,'   27.0400000000000\n');
fclose(fid);
obcount = iobs + 1

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
