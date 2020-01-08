module random_perturbation
! module contains three functions
! ran1 returns a uniform random number between 0-1
! spread returns random number between min - max
! normal returns a normal distribution

contains
    function ran1()  !returns random number between 0 - 1
        implicit none
        real*4 ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    function spread(min,max)  !returns random number between min - max
        implicit none
        real*8 spread
        real*4 min,max
        spread=(max - min) * ran1() + min
    end function spread

    function normal(mean,sigma) !returns a normal distribution
        implicit none
        real*8 normal,tmp
        real*8 mean,sigma
        integer flag
        real*8 fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0
        call init_random_seed()
            do while(rsq.ge.1.0.or.rsq.eq.0.0) ! new from for do
                r1=2.0*ran1()-1.0
                r2=2.0*ran1()-1.0
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal

    subroutine init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
    end subroutine init_random_seed

    function random_around_mean(mean,stdev)

    integer                  :: i
    integer, parameter       :: NENS=30
    real*8                   :: mean, stdev, rnum(30)
    real*8                   :: VDIF

        random_around_mean=normal(mean,stdev)
        return
    end function random_around_mean

end module random_perturbation
