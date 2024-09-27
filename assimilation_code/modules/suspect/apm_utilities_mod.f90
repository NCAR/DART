module apm_utilities_mod
   implicit none
   private
   public :: get_dist

   contains
   
real function get_dist(lat1,lat2,lon1,lon2)
! returns distance in km
   implicit none
   real:: lat1,lat2,lon1,lon2
   real:: lon_dif,rtemp
   real:: pi,ang2rad,r_earth
   real:: coef_a,coef_c
   pi=4.*atan(1.0)
   ang2rad=pi/180.
   r_earth=6371.393
! Haversine Code
   coef_a=sin((lat2-lat1)/2.*ang2rad) * sin((lat2-lat1)/2.*ang2rad) + & 
   cos(lat1*ang2rad)*cos(lat2*ang2rad) * sin((lon2-lon1)/2.*ang2rad) * &
   sin((lon2-lon1)/2.*ang2rad)
   coef_c=2.*atan2(sqrt(coef_a),sqrt(1.-coef_a))
   get_dist=abs(coef_c*r_earth)
end function get_dist

end module apm_utilities_mod

   
