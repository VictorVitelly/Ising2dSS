program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call cpu_time(starting)
  
  !Write thermalization history in a file and computes autocorrelation
  call thermalize(1.5_dp)

  !Measure energy, magnetization, susceptibility, heat capacity and binder cumulant in
  !an interval of temperatures, (initial temp., final temp, n. of points between them)
  call vary_temp(0.1_dp,2.5_dp,20)

  !Measure correlation function in an interval of temperatures
  !(initial temp., final temp, n. of points between them)
  !call correlate(2.4_dp,3._dp,4)
  
  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"


contains

  subroutine thermalize(T)
  real(dp), intent(in) :: T
  integer(i4) :: i
  integer(i4), allocatable :: spin(:,:)
  real(dp) :: vol
  open(10, file = 'data/historym.dat', status = 'replace')
  open(20, file = 'data/historye.dat', status = 'replace')
  vol=real(N**2,dp)
  allocate(spin(N,N))
    call cold_start(spin)
    !call hot_start(spin)
    do i=1,100*thermalization
      if(i==1 .or. mod(i,eachsweep)==0 ) then
        write(10,*) i, Magnet(spin)/vol
        write(20,*) i, Hamiltfbc(spin)/vol
      end if
      !call montecarlo2(spin,T)
      !call flip_sign(spin,i)
      call cluster(spin,T)
    end do
    !call autocorrelation(T,30,spin)
  close(10)
  close(20)
  deallocate(spin)
  end subroutine thermalize

  subroutine vary_temp(Ti,Tf,Nts)
  real(dp), intent(in) :: Ti,Tf
  integer(i4), intent(in) :: Nts
  integer(i4), dimension(N,N) :: spin
  integer(i4) :: i,i2,j,k
  real(dp), dimension(Nmsrs2) :: E,M,suscep,heat,U4
  real(dp) :: T,vol,norm,EE,MM,E_ave,E_delta,M_ave,M_delta,E2,M2,M4
  real(dp) :: suscep_ave,suscep_delta,heat_ave,heat_delta,U4_ave,U4_delta
  !real(dp) :: csx,csx2,cs(Nmsrs2),cs2(Nmsrs2),cs_ave,cs_delta,cs2_ave,cs2_delta
  open(10, file = 'data/energy.dat', status = 'replace')
  open(20, file = 'data/magnetization.dat', status = 'replace')
  open(30, file = 'data/susceptibility.dat', status = 'replace')
  open(40, file = 'data/heat.dat', status = 'replace')
  open(50, file = 'data/binder.dat', status = 'replace')
  !open(60, file = 'data/rank.dat', status = 'replace')
  !open(70, file = 'data/rank2.dat', status = 'replace')
  norm=real(Nmsrs,dp)
  vol=real(N**2,dp)
  do k=1,Nts
  call hot_start(spin)
    T=Ti+(Tf-Ti)*real(k-1,dp)/real(Nts-1)
    !write(*,*) k, T
    E(:)=0._dp
    M(:)=0._dp
    !cs(:)=0._dp
    !cs2(:)=0._dp
    do j=1,2*thermalization
      !call montecarlo2(spin,T)
      call cluster(spin,T)
    end do
    do j=1,Nmsrs2
      E2=0._dp
      M2=0._dp
      M4=0._dp
      do i=1,Nmsrs
        do i2=1,eachsweep
          !call montecarlo2(spin,T)
          call cluster(spin,T)
          !call cluster2(spin,T,csx,csx2)
        end do
        MM=Magnet(spin)
        EE=Hamiltfbc(spin)
        E(j)=E(j)+EE
        M(j)=M(j)+abs(MM)
        E2=E2+EE**2
        M2=M2+MM**2
        M4=M4+MM**4
        !cs(j)=cs(j)+csx
        !cs2(j)=cs2(j)+csx2
      end do
      E(j)=E(j)/norm
      M(j)=M(j)/norm
      E2=E2/norm
      M2=M2/norm
      M4=M4/norm
      suscep(j)=M2-M(j)**2
      heat(j)=E2-E(j)**2
      U4(j)=1._dp-M4/(3._dp*M2**2)
      !cs(j)=cs(j)/norm
      !cs2(j)=cs2(j)/norm
    end do
    call mean_scalar(E,E_ave,E_delta)
    call mean_scalar(M,M_ave,M_delta)
    call mean_scalar(suscep,suscep_ave,suscep_delta)
    call mean_scalar(heat,heat_ave,heat_delta)
    call mean_scalar(U4,U4_ave,U4_delta)
    !call mean_scalar(cs,cs_ave,cs_delta)
    !call mean_scalar(cs2,cs2_ave,cs2_delta)
    write(10,*) T, E_ave/vol, E_delta/vol
    write(20,*) T, M_ave/vol, M_delta/vol
    write(30,*) T, suscep_ave/vol, suscep_delta/vol
    write(40,*) T, heat_ave/vol, heat_delta/vol
    write(*,*) T, U4_ave, U4_delta
    !write(60,*) T, cs_ave,cs_delta
    !write(70,*) T, cs2_ave/(vol**2), cs2_delta/(vol**2)
  end do
  
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  !close(60)
  !close(70)
  end subroutine vary_temp

end program main
