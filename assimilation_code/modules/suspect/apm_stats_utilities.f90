module apm_stats_utilities
  
   implicit none
   private

   public :: median_code, &
             mode_code, &
             vertical_column

   contains

subroutine median_code(fld_median,fld,nfld)
   implicit none
   integer                 :: ipt,iptt,nfld
   real                    :: fld_median,temp
   real, dimension(nfld)   :: fld
   real, dimension(1000)   :: temp_arr
!
   temp_arr(1:nfld)=fld(1:nfld)
   do ipt=2,nfld
      if(temp_arr(ipt) .lt. temp_arr(ipt-1)) then
         iptt=ipt
         do while (temp_arr(iptt) .lt. temp_arr(iptt-1))
            temp=temp_arr(iptt-1)
            temp_arr(iptt-1)=temp_arr(iptt)
            temp_arr(iptt)=temp
            iptt=iptt-1
            if(iptt .eq. 1) exit
         enddo
      endif
   enddo
   if(nfld/2*2 .eq.nfld) then
      fld_median=(temp_arr(nfld/2)+temp_arr(nfld/2+1))/2.
   else
      fld_median=temp_arr(nfld/2+1)
   endif
   return
end subroutine median_code
!
subroutine mode_code(fld_mode,fld,nfld)
   implicit none
   integer                  :: ipt,iptt,nfld
   integer                  :: cnt_max,cnt_idx
   integer, dimension(nfld) :: cnt
   real                     :: fld_mode,temp
   real                     :: bin_tol
   real, dimension(nfld)    :: fld
   real, dimension(1000)    :: temp_arr
!
   bin_tol=.1   
   temp_arr(1:nfld)=fld(1:nfld)
   do ipt=2,nfld
      if(temp_arr(ipt) .lt. temp_arr(ipt-1)) then
         iptt=ipt
         do while (temp_arr(iptt) .lt. temp_arr(iptt-1))
            temp=temp_arr(iptt-1)
            temp_arr(iptt-1)=temp_arr(iptt)
            temp_arr(iptt)=temp
            iptt=iptt-1
            if(iptt .eq. 1) exit
         enddo
      endif
   enddo
!
   cnt(:)=1
   do ipt=1,nfld
      do iptt=1,nfld
         if(ipt.eq.iptt) cycle
         if(temp_arr(ipt).eq.temp_arr(iptt)) cnt(ipt)=cnt(ipt)+1
      enddo
   enddo   
   cnt_max=cnt(1)
   cnt_idx=1
   do ipt=2,nfld
      if(cnt(ipt).gt.cnt_max) then
         cnt_max=cnt(ipt)
         cnt_idx=ipt
      endif
   enddo
!
   fld_mode=0.    
   if(cnt_max.gt.1) then
      fld_mode=temp_arr(cnt_idx)
   endif
   return
end subroutine mode_code
!
subroutine vertical_column(vert_col,prof,prs_lev,nlay)
   implicit none
   integer                 :: k,nlay
   real                    :: wt,twt
   real                    :: vert_col
   real, dimension(nlay)   :: prof
   real, dimension(nlay+1) :: prs_lev
!
   twt=0.   
   vert_col=0.
   do k=1,nlay
      wt=prs_lev(k)-prs_lev(k+1)
      twt=twt+wt
      vert_col=vert_col+wt*prof(k)
   enddo
   vert_col=vert_col/twt
   return
end subroutine vertical_column

end module apm_stats_utilities
