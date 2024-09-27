module apm_cpsr_mod
  
   implicit none
   private

   public :: cpsr_calculation, &
             mat_prd, &
             mat_tri_prd, &
             vec_to_mat, &
             diag_inv_sqrt, &
             lh_mat_vec_prd, &
             rh_vec_mat_prd, &
             mat_transpose, &
             diag_vec

   contains
  
subroutine cpsr_calculation(rnlayer_trc,nlayer,spec_cpsr,avgk_cpsr,prior_cpsr,spec_obs,prior_obs,avgk_obs,cov_obs,cnt) 
   real,parameter                               :: eps_tol=1.e-3
   integer                                      :: i,j,k,kk,cnt
   integer                                      :: info,nlayer,nlayer_trc,rnlayer_trc,qstatus
   real                                         :: sdof,sum
   integer                                      :: lwrk
   real                                         :: spec_cpsr(nlayer),prior_cpsr(nlayer),avgk_cpsr(nlayer,nlayer)
   real                                         :: spec_obs(nlayer),spec_obs_adj(nlayer),prior_obs(nlayer)
   real                                         :: avgk_obs(nlayer,nlayer),cov_obs(nlayer,nlayer)
   real                                         :: avgk_obs_adj(nlayer,nlayer),prior_obs_adj(nlayer)
   double precision,allocatable,dimension(:)    :: wrk
   double precision,allocatable,dimension(:)    :: ZV,ZW,SV_cov
   double precision,allocatable,dimension(:)    :: cp_x_r,cp_x_p
   double precision,allocatable,dimension(:)    :: rs_x_r,rs_x_p
   double precision,allocatable,dimension(:)    :: err2_rs_r
   double precision,allocatable,dimension(:,:)  :: Z,ZL,ZR,SV,U_cov,V_cov,UT_cov,VT_cov
   double precision,allocatable,dimension(:,:)  :: cp_avgk,cp_cov
   double precision,allocatable,dimension(:,:)  :: rs_avgk,rs_cov
!
   lwrk=5*nlayer
   allocate(wrk(lwrk))
   allocate(Z(nlayer,nlayer),SV_cov(nlayer),SV(nlayer,nlayer))
   allocate(U_cov(nlayer,nlayer),UT_cov(nlayer,nlayer),V_cov(nlayer,nlayer),VT_cov(nlayer,nlayer))
   allocate(cp_avgk(nlayer,nlayer),cp_cov(nlayer,nlayer),cp_x_r(nlayer),cp_x_p(nlayer))
   allocate(rs_avgk(nlayer,nlayer),rs_cov(nlayer,nlayer),rs_x_r(nlayer),rs_x_p(nlayer))       
   allocate(ZL(nlayer,nlayer),ZR(nlayer,nlayer),ZV(nlayer),ZW(nlayer))
   allocate(err2_rs_r(nlayer))
!
! COMPRESSION STEP   
! Calculate SVD of the averaging kernel (Z=U_xxx * SV_xxx * VT_xxx)
   Z(1:nlayer,1:nlayer)=dble(avgk_obs(1:nlayer,1:nlayer))
   call dgesvd('A','A',nlayer,nlayer,Z,nlayer,SV_cov,U_cov,nlayer,VT_cov,nlayer,wrk,lwrk,info)
   nlayer_trc=0
   sdof=0.
!
! Phase space truncation
   do k=1,nlayer
      if(SV_cov(k).ge.eps_tol) then
         nlayer_trc=k
         sdof=sdof+SV_cov(k)
      else
         SV_cov(k)=0
         U_cov(:,k)=0. 
         VT_cov(k,:)=0.
      endif 
   enddo
!
! Get the transpose      
   call mat_transpose(U_cov,UT_cov,nlayer,nlayer)
   call mat_transpose(VT_cov,V_cov,nlayer,nlayer)
   call vec_to_mat(SV_cov,SV,nlayer)
!
! Compress the terms in the forward operator
! Get the compressed averaging kernel (cp_avgk)
   ZL(1:nlayer,1:nlayer)=dble(avgk_obs(1:nlayer,1:nlayer))
   call mat_prd(UT_cov(1:nlayer,1:nlayer),ZL(1:nlayer,1:nlayer), &
   cp_avgk(1:nlayer,1:nlayer),nlayer,nlayer,nlayer,nlayer)
!   
! Get the compressed observation error covariance (cp_cov) 
   ZL(1:nlayer,1:nlayer)=dble(cov_obs(1:nlayer,1:nlayer))
   call mat_tri_prd(UT_cov(1:nlayer,1:nlayer),ZL(1:nlayer,1:nlayer),U_cov(1:nlayer,1:nlayer), &
   cp_cov(1:nlayer,1:nlayer),nlayer,nlayer,nlayer,nlayer,nlayer,nlayer)
!   
! Get the adjusted avgering kerrnel (avgk_obs_adj) A_adj = (I-A)  
   do k=1,nlayer
      do kk=1,nlayer
         avgk_obs_adj(k,kk)=-1.*avgk_obs(k,kk)
      enddo
      avgk_obs_adj(k,k)=avgk_obs_adj(k,k)+1.
   enddo
!
! Get the prior term (prior_obs_adj): x_p_adj = (I-A) x_p
   call lh_mat_vec_prd(dble(avgk_obs_adj(1:nlayer,1:nlayer)),dble(prior_obs(1:nlayer)), &
   ZW(1:nlayer),nlayer)
   prior_obs_adj(1:nlayer)=real(ZW(1:nlayer))
!
! Get the quasi-optimal retrieval (QOR) term (spec_obs_adj): x_r_qor = x_r - (I-A) x_p
   spec_obs_adj(1:nlayer)=spec_obs(1:nlayer)-prior_obs_adj(1:nlayer)
   ZV(1:nlayer)=dble(spec_obs_adj(1:nlayer))
   call lh_mat_vec_prd(UT_cov(1:nlayer,1:nlayer),ZV(1:nlayer),cp_x_r(1:nlayer),nlayer)
!
! Get the compressed prior (cp_x_p)   
   ZV(1:nlayer)=dble(prior_obs_adj(1:nlayer))
   call lh_mat_vec_prd(UT_cov(1:nlayer,1:nlayer),ZV(1:nlayer),cp_x_p(1:nlayer),nlayer)
!
! ROTATION STEP
! Calculate SVD of compressed error covariance (cp_cov) (Z=U_xxx * SV_xxx * VT_xxx)
   Z(1:nlayer,1:nlayer)=cp_cov(1:nlayer,1:nlayer)
   call dgesvd('A','A',nlayer,nlayer,Z,nlayer,SV_cov,U_cov,nlayer,VT_cov,nlayer,wrk,lwrk,info)
   do k=nlayer_trc+1,nlayer
      SV_cov(k)=0
      U_cov(:,k)=0. 
      VT_cov(k,:)=0.
   enddo
!
! Get the transpose   
   call mat_transpose(U_cov,UT_cov,nlayer,nlayer)
   call mat_transpose(VT_cov,V_cov,nlayer,nlayer)
   call vec_to_mat(SV_cov,SV,nlayer)
!
! Check SVD of the compressed error covariance (cp_cov). It may not be full rank.
   ZL(1:nlayer,1:nlayer)=cp_cov(1:nlayer,1:nlayer)
   call mat_tri_prd(UT_cov(1:nlayer,1:nlayer),ZL(1:nlayer,1:nlayer),U_cov(1:nlayer,1:nlayer), &
   rs_cov(1:nlayer,1:nlayer),nlayer,nlayer,nlayer,nlayer,nlayer,nlayer)
   rnlayer_trc=0   
   do k=1,nlayer
      if(rs_cov(k,k).gt.0) then
         rnlayer_trc=rnlayer_trc+1
      endif
      if(rs_cov(k,k).lt.0) then
         exit
      endif   
   enddo
   do k=rnlayer_trc+1,nlayer
      SV_cov(k)=0
      U_cov(:,k)=0. 
      VT_cov(k,:)=0.
   enddo
!
! Scale the singular vectors
   do k=1,rnlayer_trc
      U_cov(:,k)=U_cov(:,k)/sqrt(SV_cov(k))
   enddo
!
! Get the transpose   
   call mat_transpose(U_cov,UT_cov,nlayer,nlayer)
   call mat_transpose(VT_cov,V_cov,nlayer,nlayer)
   call vec_to_mat(SV_cov,SV,nlayer)
!
! Rotate terms in the forward operator
! Get the compressed and rotated averaging kernel (rs_avgk)   
   ZL(1:nlayer,1:nlayer)=cp_avgk(1:nlayer,1:nlayer)
   call mat_prd(UT_cov(1:nlayer,1:nlayer),ZL(1:nlayer,1:nlayer), &
   rs_avgk(1:nlayer,1:nlayer),nlayer,nlayer,nlayer,nlayer)
!
! Get the compressed and rotated quasi-optimal retrieval (rs_x_r)   
   call lh_mat_vec_prd(UT_cov(1:nlayer,1:nlayer),cp_x_r(1:nlayer),rs_x_r(1:nlayer),nlayer)
!
! Get the compressed and rotated prior (rs_x_p)   
   call lh_mat_vec_prd(UT_cov(1:nlayer,1:nlayer),cp_x_p(1:nlayer),rs_x_p(1:nlayer),nlayer)

! Get new errors. These should all be 1.0 outside the null space and 0.0 in the null space.
   err2_rs_r(:)=0.
   do k=1,rnlayer_trc
      err2_rs_r(k)=sqrt(rs_cov(k,k))
   enddo
!
! Assign variables for return to the calling routine
   spec_cpsr(1:nlayer)=real(rs_x_r(1:nlayer))
   prior_cpsr(1:nlayer)=real(rs_x_p(1:nlayer))
   avgk_cpsr(1:nlayer,1:nlayer)=real(rs_avgk(1:nlayer,1:nlayer))
!
! Clean up and return   
   deallocate(wrk)
   deallocate(Z,SV_cov,SV)
   deallocate(U_cov,UT_cov,V_cov,VT_cov)
   deallocate(cp_avgk,cp_cov,cp_x_r,cp_x_p)       
   deallocate(rs_avgk,rs_cov,rs_x_r,rs_x_p)       
   deallocate(ZL,ZR,ZV,ZW)
end subroutine cpsr_calculation

!-------------------------------------------------------------------------------

subroutine mat_prd(A_mat,B_mat,C_mat,na,ma,nb,mb)
!
! compute dot product of two matrics
   integer :: ma,na,mb,nb,i,j,k
   double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(na,mb)
!
! check that na=mb
   if(ma .ne. nb) then
      print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
      stop
   endif
!
! initialze the product array
   C_mat(:,:)=0.
!
! calculate inner product
   do i=1,na
      do j=1,mb
         do k=1,ma
            C_mat(i,j)=C_mat(i,j)+A_mat(i,k)*B_mat(k,j) 
         enddo
      enddo
   enddo
end subroutine mat_prd

!-------------------------------------------------------------------------------

subroutine mat_tri_prd(A_mat,B_mat,C_mat,D_mat,na,ma,nb,mb,nc,mc)
!
! compute dot product of three matrics D=A*B*C
   integer :: na,ma,nb,mb,nc,mc,i,j,k
   double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(nc,mc),D_mat(na,mc)
   double precision :: Z_mat(nb,mc)
!
! check that na=mb
   if(ma .ne. nb) then
      print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
      stop
   endif
   if(mb .ne. nc) then
      print *, 'Error in matrix dimension mb (cols) must equal nc (rows) ',mb,' ',nc
      stop
   endif
!
! initialze the product array
   Z_mat(:,:)=0.
   D_mat(:,:)=0.
!
! calculate first inner product Z=B*C
   do i=1,nb
      do j=1,mc
         do k=1,mb
            Z_mat(i,j)=Z_mat(i,j)+B_mat(i,k)*C_mat(k,j) 
         enddo
      enddo
   enddo
!
! calculate second inner product D=A*Z
   do i=1,na
      do j=1,mc
         do k=1,ma
            D_mat(i,j)=D_mat(i,j)+A_mat(i,k)*Z_mat(k,j) 
         enddo
      enddo
   enddo
end subroutine mat_tri_prd

!-------------------------------------------------------------------------------

subroutine vec_to_mat(a_vec,A_mat,n)
!
! compute dot product of two matrics
   integer :: n,i
   double precision :: a_vec(n),A_mat(n,n)
!
! initialze the product array
   A_mat(:,:)=0.
!
! calculate inner product
   do i=1,n
      A_mat(i,i)=a_vec(i) 
   enddo
end subroutine vec_to_mat

!-------------------------------------------------------------------------------

subroutine diag_inv_sqrt(A_mat,n)
!
! calculate inverse square root of diagonal elements
   integer :: n,i
   double precision :: A_mat(n,n)
   do i=1,n
      if(A_mat(i,i).le.0.) then
         print *, 'Error in Subroutine vec_to_mat arg<=0 ',i,' ',A_mat(i,i)
         call abort
      endif
      A_mat(i,i)=1./sqrt(A_mat(i,i)) 
   enddo
   return
end subroutine diag_inv_sqrt

!-------------------------------------------------------------------------------

subroutine diag_sqrt(A_mat,n)
!
! calculate square root of diagonal elements
   integer :: n,i
   double precision :: A_mat(n,n)
   do i=1,n
      if(A_mat(i,i).lt.0.) then
         print *, 'Error in Subroutine vec_to_mat arg<0 ',i,' ',A_mat(i,i)
         call abort
      endif
      A_mat(i,i)=sqrt(A_mat(i,i)) 
   enddo
end subroutine diag_sqrt

!-------------------------------------------------------------------------------

subroutine lh_mat_vec_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate left hand side scaling of column vector
   integer :: n,i,j
   double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
   s_a_vec(:)=0.
!
! conduct scaling
   do i=1,n
      do j=1,n
         s_a_vec(i)=s_a_vec(i)+SCL_mat(i,j)*a_vec(j)
      enddo 
   enddo
end subroutine lh_mat_vec_prd

!-------------------------------------------------------------------------------

subroutine rh_vec_mat_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate right hand side scaling of a row vector
   integer :: n,i,j
   double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
   s_a_vec(:)=0.
!
! conduct scaling
   do i=1,n
      do j=1,n
         s_a_vec(i)=s_a_vec(i)+a_vec(j)*SCL_mat(j,i) 
      enddo
   enddo
end subroutine rh_vec_mat_prd

!-------------------------------------------------------------------------------

subroutine mat_transpose(A_mat,AT_mat,n,m)
!
! calculate matrix transpose
   integer :: n,m,i,j
   double precision :: A_mat(n,m),AT_mat(m,n)
   do i=1,n
      do j=1,m
         AT_mat(j,i)=A_mat(i,j) 
      enddo
   enddo
end subroutine mat_transpose

!-------------------------------------------------------------------------------

subroutine diag_vec(A_mat,a_vec,n)
!
! calculate square root of diagonal elements
   integer :: n,i
   double precision :: A_mat(n,n),a_vec(n)
   do i=1,n
      a_vec(i)=A_mat(i,i) 
   enddo
end subroutine diag_vec

end module apm_cpsr_mod
