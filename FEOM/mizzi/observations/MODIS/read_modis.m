% reads MODIS HDF files, process the data and dump the data as ascii file
% the output of this code is the input of modis_ascii_to_obs in DART
%
   clear all;
   iday_start=29;
   iday_end=30;
   ifile_start=1;
%
% where's the data?
   rootdir='/glade/p/work/mizzi/TRUNK/DART_CHEM/observations/MODIS';
   rundir='/glade/scratch/mizzi/MODIS_ASCII';
   date_dir=2008060100:100:2008063000;
   indx=1;
%   for indx=1:30
   for indx=iday_start:iday_end
      cd(rundir);
      date_str=int2str(date_dir(indx));
      datadir=strcat('/glade/p/acd/mizzi/AVE_TEST_DATA/MODIS/',date_str);
%
% directory listing
      list=dir(datadir);
%
% cut other files not HDF
      list=list(3:end);
%      list(1).name
%      list(2).name
%      list(3).name
%
% file listing
      filelist={list.name};
%
% number of files
      nfiles=length(filelist);
%
% for now force it, since i do not know how to read a local attribute   
      scale_factor=0.001;
%
% time conversion (MODIS time is at TAI time seconds since 1993-1-1 00:00)
      time_factor = 12*365+3*366; %? how many days up to current year(2008) since 1993
      sec2day     = 1/(24*60*60);
      month_array = [31,29,31,30,31,30,31,31,30,31,30,31]; % 2008 is a leap year
%
% returns the cumulative sum along different dimensions of an array.
      cumdays = cumsum(month_array);
      year = 2008;
      otype = 1;
      vert  = 0;
      file_output=strcat('/glade/scratch/mizzi/MODIS_ASCII/modis_ascii_',date_str,'.input')
      fid = fopen(file_output,'w');
%
% MAIN LOOP
      file_record=strcat('/glade/scratch/mizzi/MODIS_ASCII/matlab_record.txt')
      fif = fopen(file_record,'w');
      ifile_end=nfiles;
      for ifile=ifile_start:ifile_end
         cd(datadir);
%    
% Convert cell array of matrices to single matrix
         filename=cell2mat(filelist(ifile));
         fprintf(fif,'Day %d, File %d, NFiles %d, Filename %s\n',indx,ifile,nfiles,filename);
%    
% filelist(1) is  'MYD04_L2.A2010235.0000.051.2010240042123.hdf'    
% after  filename=cell2mat(filelist(1)) 
% filename is      MYD04_L2.A2010235.0000.051.2010240042123.hdf
% up to this point, it only involve filename operation, nothing to do
% with actual file content
%    
         fileinfo = hdfinfo(filename);
%
% S = hdfinfo(filename) returns a structure S whose fields contain 
% information about the contents of an HDF4 or HDF-EOS file.
% filename is a string that specifies the name of the HDF4 file
%     fileinfo = 
%     Filename: [1x81 char]
%     Attributes: [1x8 struct]
%     Vgroup: [1x1 struct]
%
% disp(X) displays an array, without printing the array name. 
% If X contains a text string, the string is displayed.        
         disp(filename); 
%    
% variables of interest
         var_int={'Scan_Start_Time', 'Latitude', 'Longitude', 'Cloud_Mask_QA', ...
                  'Quality_Assurance_Land', 'Quality_Assurance_Ocean', ...
                  'Optical_Depth_Land_And_Ocean', ...
                  'Deep_Blue_Aerosol_Optical_Depth_550_Land'};
%
% var_int(1) will be 'Scan_Start_Time' , cell2mat(var_int(1)) will be Scan_Start_Time
         time=double(hdfread(filename,cell2mat(var_int(1))));
%
% time is 203 X 135 array 
         lat=double(hdfread(filename,cell2mat(var_int(2))));
         lon=double(hdfread(filename,cell2mat(var_int(3))));
%cloud_mask_qa
         cldmask=hdfread(filename,cell2mat(var_int(4)));
         qa_lnd=hdfread(filename,cell2mat(var_int(5)));
         qa_ocn=hdfread(filename,cell2mat(var_int(6)));  
         tau=double(hdfread(filename,cell2mat(var_int(7))))*scale_factor;
% Optical_Depth_Land_And_Ocean
         taud=double(hdfread(filename,cell2mat(var_int(8))))*scale_factor;
% Deep_Blue_Aerosol_Optical_Depth_550_Land 
%
% Cloud Mask Status Flag  0: Undetermined  1: determined
% qa_mask(Cloud Mask Cloudiness Flag) : 0: 0-30% cloudy, 1: 30-60%, 2: 60-90%, 3: >90%  
% surf_flag: 0: ocean, 1: coast, 2: desert, 3: land
% deep blue aerosol usefulness flag dbuf: 1 useful
% deep blue aerosol confidence flag dbcf: 0 no conf, 1: marginal, 2: good, 3: very good
% deep blue aerosol type dbat: 0 mixed, 1: dust, 2:smoke, 3:sulfate
% deep blue aerosol retrieving cond dbrc: 0: optimal,1: white sand
% 2:cloudy, 3 t(550nm)>5.0
%   
% loop for each array
         [m,n]=size(time);
         for i=1:m
            for j=1:n            
               bin  = dec2bin(cldmask(i,j),8);
%
% Convert decimal to binary number in string
%              str = dec2bin(d,n) produces a binary representation with at least n bits.
%
% Cloud Mask Status Flag, I think equal to str2double(bin(8))
               mask = bin2dec( num2str(bin(8)) );
%          
% The num2str function converts numbers to their string representations.
% This function is useful for labeling and titling plots with numeric values.
% Cloud Mask Cloudiness Flag
               qa_mask = bin2dec( num2str(bin([6:7])) );
%
%Surface Type Flag            
               surf_flag = bin2dec( num2str(bin([2:3])) );
               bin = dec2bin(qa_lnd(i,j,5), 8);
               dbrc = bin2dec( num2str(bin([2:3])) );
               dbat = bin2dec( num2str(bin([4:5])) );
               dbcf = bin2dec( num2str(bin([6:7])) );
               dbuf = bin2dec( num2str(bin([8])) );
               aod=tau(i,j);
               aodlon=lon(i,j);
               aodlat=lat(i,j);
               aodtim=time(i,j);
%            
% qc (may vary with choice)
%               if ( aod>0.055 & dbuf == 1 & dbcf >=1 & aodtim > 0 ) % & (dbrc==0 | dbrc==2) ) 
               if ( aod>0.055 && aod<1.5 && aodtim > 0 && qa_mask==0)
%
% convert time to year month day hour minute second UTC
% need to plus one, because if , ex. 0.5, but it should be 1.5 days within the month 
                  tmp = ( time(i,j) - time_factor/sec2day )*sec2day + 1;
                  tmp2 = find(cumdays >= tmp); 
                  month=tmp2(1,1);
                  if (month==1)
                     days=tmp;
                  else
                     days= ( tmp-cumdays(month-1) );
                  end
%
% B = floor(A) rounds the elements of A to the nearest integers less than or equal to A. 
% For complex A, the imaginary and real parts are rounded independently.
                  day = floor(days); 
                  hours= (days-day)*24; 
                  hour=floor(hours);
                  minutes=(hours-hour)*60;
                  minute=floor(minutes);
                  seconds=(minutes-minute)*60;
                  second=floor(seconds);
%
% seconds,minutes,hours,days are all decimals.     
% subjective estimate error
                  if (aod< 0.2)
                     aoderr = 0.1;
                  elseif (aod<1.4 && aod >=0.2)
                     aoderr = 0.05 + 0.20*aod;
                  elseif (aod>=1.4)
                     aoderr = 0.2+0.4*aod;
                  end
%            
% print stuff
                  fprintf(fid,'%2d %10.4f %10.4f %5d %5d %5d %5d %5d %5d %5d %8.4f %8.4f\n', ...
                  otype, aodlat, aodlon, vert, year, month, day, hour, minute, second, aod, aoderr);
%               
               end		%end if
            end
         end              
      end
      fclose(fid);
   end
