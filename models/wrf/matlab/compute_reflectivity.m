function ref = compute_reflectivity( qr, qg, qs, rho, temp )
%% compute_reflectivity
%
% Computes reflectivity (dBz).
%
% Inputs: qr = rain water mixing ratio (3d)
%         qg =    graupel mixing ratio (3d)
%         qs =       snow mixing ratio (3d)
%	 rho = density,      at mass pts
%       temp = temperature,  at mass pts 
% Output:
%	 ref = reflectivity, at mass pts

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

dief  = 0.224;
n0r   = 8.0e6;
n0g   = 4.0e4;
n0s   = 3.0e6;
rho_r = 1000.0;
rho_g = 917.0;
rho_s = 100.0;

[Nk Nj Ni] = size(temp) ; 

ref = zeros(Nk, Nj, Ni) ;

ar = 7.2e20/(((pi*rho_r)^1.75)*(n0r^0.75)) ;
ag_dry = dief*((rho_g/rho_r)^2)*7.2e20 / (((pi*rho_g)^1.75)*(n0g^0.75)) ;

%% This is appropriate for 10-cm radar.
ag_wet = (7.2e20/(((pi*rho_g)^1.75)*(n0g^0.75)))^0.95 ;
as_wet = 7.2e20/(((pi*rho_s)^1.75)*(n0s^0.75)) ;
as_dry = dief*((rho_s/rho_r)^2)*as_wet ;

%% RAIN
precip = rho .* qr ;

for ii = 1:Ni
for jj = 1:Nj
for kk = 1:Nk
   if( precip(kk,jj,ii) > 0.0 )
       ref(kk,jj,ii) = ref(kk,jj,ii) + ar * (precip(kk,jj,ii)^1.75);
   end
end
end
end

%% HAIL / GRAUPEL
precip = rho .* qg ;

for ii = 1:Ni
for jj = 1:Nj
for kk = 1:Nk
   if( precip(kk,jj,ii) > 0.0 )
      if ( temp(kk,jj,ii) < 273.15 )
         ref(kk,jj,ii) = ref(kk,jj,ii) + ag_dry * (precip(kk,jj,ii)^1.75) ;
      else
         ref(kk,jj,ii) = ref(kk,jj,ii) + ag_wet * (precip(kk,jj,ii)^1.6625) ;
      end
   end
end
end
end

%% SNOW
precip = rho .* qs ;

for ii = 1:Ni
for jj = 1:Nj
for kk = 1:Nk
   if( precip(kk,jj,ii) > 0.0 )
      if ( temp(kk,jj,ii) < 273.15 )
         ref(kk,jj,ii) = ref(kk,jj,ii) + as_dry * (precip(kk,jj,ii)^1.75) ;
      else
         ref(kk,jj,ii) = ref(kk,jj,ii) + as_wet * (precip(kk,jj,ii)^1.75) ;
      end
   end
end
end
end

for ii = 1:Ni
for jj = 1:Nj
for kk = 1:Nk
   if( ref(kk,jj,ii) > 0.0 )
      ref(kk,jj,ii) = 10.0 * log10(ref(kk,jj,ii)) ;
   else
      ref(kk,jj,ii) = NaN ;
   end
end
end
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
