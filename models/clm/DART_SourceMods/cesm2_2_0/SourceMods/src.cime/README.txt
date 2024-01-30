
cime/src/drivers/mct/main/seq_rest_mod.F90 no longer needs modifications.

! Jim Edwards 6/7/17: "It turns out that a diagnostic field in the coupler was changed 
! (from cesm1_5) in a non-backward compatible way.  But I have a workaround.  On your 
! first run comment out lines 258-276 of file seq_rest_mod.F90 that will get you 
! started and after the first restart you can remove that modification."

