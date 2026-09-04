module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine initialize(corr1,corr2)
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize

  subroutine correlation_function(corr1,corr2,CF)
    real(dp), dimension(N), intent(in) :: corr1
    real(dp), dimension(N,N), intent(in) :: corr2
    real(dp), dimension(N), intent(out) :: CF
    integer(i4) :: i1
    do i1=1,N
      CF(i1)=corr2(i1,1)-corr1(i1)*corr1(1)
    end do
  end subroutine correlation_function

  subroutine autocorrelation(T,tmax,spin)
    integer(i4), intent(in) :: tmax
    real(dp), intent(in) :: T
    integer(i4), dimension(N,N), intent(inout) :: spin
    real(dp), dimension(tmax+1) :: auto,auto_delta
    real(dp) :: E(Nmsrs+tmax), auto1(Nmsrs)
    real(dp) :: E_ave,auto1_ave,autoj(tmax+1,Nauto)
    integer(i4) :: i,j,tt
    open(70, file = 'data/autocorr.dat', status = 'replace')
    do j=1,Nauto
      do i=1,Nmsrs+tmax
        !call montecarlo(spin,T)
        call cluster(spin,T)
        E(i)=Hamilt(spin)/(N**2)
      end do
      call mean_0(E,E_ave )
      
      do tt=0,tmax
        do i=1,Nmsrs
          auto1(i)=E(i)*E(i+tt)
        end do
        call mean_0(auto1,auto1_ave)
        auto=auto1_ave-(E_ave**2)
        autoj(tt+1,j)=auto1_ave-(E_ave**2)
      end do
    end do
    do tt=0,tmax
      call mean_scalar(autoj(tt+1,:),auto(tt+1),auto_delta(tt+1))
      write(70,*) tt,auto(tt+1),auto_delta(tt+1)
    end do
    close(70)
  end subroutine autocorrelation

end module measurements
