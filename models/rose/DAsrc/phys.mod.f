c----------------------------------------------------------------------
      module phys 
c----------------------------------------------------------------------
             
      use params

      implicit none

      save

c  common blocks for heating, cooling and gravity wave param.
      real :: tcltop, co2mr(nzz), rair, gravit
      real, dimension (nz,nx,ny) :: fxh, fyh, fdzzh,
     +                              fgr, fcgr, falph, fdzz,
     +                              solht, fsponge, q_resid
      real, dimension (nztrop,ny) :: tpnmc, unmc, tnmc
      real, dimension (nztrop) :: troptref, tropsref
      real, dimension (nz) :: dtzz, zkmin, xml2, wmole, dtglob, dtadv

      real, dimension (nz,ny) :: cl_o3, cl_o3t, cl_o2, cl_o2t, uclimo
      real, dimension (nzz,ny) :: cl_o3zz, cl_h2o, cl_o3p, cl_co2
      real, dimension (nzz,nx,ny) :: qir 
      real, dimension (nz,nbcon) :: dmolch, wdif

      real, dimension (nz,nx,ny) :: ch_heat, E1, E2, E3, E4, E5, E6, E7 

      end module phys
