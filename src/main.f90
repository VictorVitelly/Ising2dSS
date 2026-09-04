program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call init_vecs()
  call cpu_time(starting)
  
  !Write thermalization history in a file and computes autocorrelation
  !call thermalize(2.5_dp)

  !Measure energy, magnetization, susceptibility, heat capacity and binder cumulant in
  !an interval of temperatures, (initial temp., final temp, n. of points between them)
  call vary_temp(0.1_dp,0.5_dp,30)

  !Measure correlation function in an interval of temperatures
  !(initial temp., final temp, n. of points between them)
  !call correlate(2.0_dp,3._dp,16)
  
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
    !call cold_start(spin)
    !call hot_start(spin)
    spin=-1
    do i=1,10*thermalization
      if(i==1 .or. mod(i,10)==0 ) then
        write(10,*) i, Magnet(spin)/vol
        write(20,*) i, Hamilt(spin)/vol
      end if
      !call montecarlo2(spin,T)
      call montecarlo(spin,T)
      !call flip_sign(spin,i)
      !call cluster(spin,T)
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
  real(dp), dimension(Nmsrs2) :: E,M,suscep,heat,U4,AR
  real(dp) :: T,vol,norm,EE,MM,E_ave,E_delta,M_ave,M_delta,E2,M2,M4,arr,AR_ave,AR_err
  real(dp) :: suscep_ave,suscep_delta,heat_ave,heat_delta,U4_ave,U4_delta
  !real(dp) :: csx,csx2,cs(Nmsrs2),cs2(Nmsrs2),cs_ave,cs_delta,cs2_ave,cs2_delta
  real(dp), allocatable :: corr1(:),corr2(:,:),CF(:,:,:),CF_ave(:,:),CF_err(:,:)
  real(dp) :: xi2_ave,xi2_err
  open(10, file = 'data/ene.dat', status = 'replace')
  open(20, file = 'data/mag.dat', status = 'replace')
  open(30, file = 'data/sus.dat', status = 'replace')
  open(40, file = 'data/hea.dat', status = 'replace')
  open(50, file = 'data/bin.dat', status = 'replace')
  open(60, file = 'data/art.dat', status = 'replace')
  open(70, file = 'data/xi2.dat', status = 'replace')
  !open(60, file = 'data/rank.dat', status = 'replace')
  !open(70, file = 'data/rank2.dat', status = 'replace')
  
  allocate(corr2(N,N))
  allocate(corr1(N))
  allocate(CF(N,N,Nmsrs2))
  allocate(CF_ave(N,Nts))
  allocate(CF_err(N,Nts))
  
  norm=real(Nmsrs,dp)
  vol=real(N**2,dp)
  !call cold_start(spin)
  
  do k=1,Nts
  call hot_start(spin)
  !call cold_start(spin)
  !spin(:,:)=1
    T=Ti+(Tf-Ti)*real(k-1,dp)/real(Nts-1)
    write(*,*) k, 'de', Nts
    E(:)=0._dp
    M(:)=0._dp
    AR(:)=0._dp
    CF(:,:,:)=0._dp
    !cs(:)=0._dp
    !cs2(:)=0._dp
    EE=Hamilt(spin)
    do j=1,thermalization
      !call montecarlo(spin,T)
      call montecarlobien(spin,T,EE)
      !call cluster(spin,T)
    end do
    do j=1,Nmsrs2
      E2=0._dp
      M2=0._dp
      M4=0._dp
      corr1(:)=0._dp
      corr2(:,:)=0._dp
      do i=1,Nmsrs
        do i2=1,eachsweep
          !call montecarlo(spin,T)
          call montecarlobien(spin,T,EE)
          !call cluster(spin,T)
          !call cluster2(spin,T,csx,csx2)
        end do
        call metropolisbien(spin,T,EE,arr)
        MM=Magnet(spin)
        EE=Hamilt(spin)
        AR(j)=AR(j)+arr
        E(j)=E(j)+EE
        M(j)=M(j)+abs(MM)
        E2=E2+EE**2
        M2=M2+MM**2
        M4=M4+MM**4
        call correlation(spin,corr1,corr2)
        !cs(j)=cs(j)+csx
        !cs2(j)=cs2(j)+csx2
      end do
      E(j)=E(j)/norm
      M(j)=M(j)/norm
      AR(j)=AR(j)/norm
      E2=E2/norm
      M2=M2/norm
      M4=M4/norm
      suscep(j)=M2-M(j)**2
      heat(j)=E2-E(j)**2
      U4(j)=1._dp-M4/(3._dp*M2**2)
      do i=1,N
        do i2=1,N
          CF(i,i2,j)=corr2(i,i2)/(norm)-(corr1(i)*corr1(i2))/(norm**2)
          !write(*,*) corr2(1,2)/norm, corr1(1)/norm,corr1(2)/norm
        end do
      end do  
      !write(*,*) CF(:,:,j)
      !cs(j)=cs(j)/norm
      !cs2(j)=cs2(j)/norm
    end do
    call mean_scalar(E,E_ave,E_delta)
    call mean_scalar(M,M_ave,M_delta)
    call mean_scalar(suscep,suscep_ave,suscep_delta)
    call mean_scalar(heat,heat_ave,heat_delta)
    call mean_scalar(U4,U4_ave,U4_delta)
    call mean_scalar(AR,AR_ave,AR_err)
    !call mean_scalar(cs,cs_ave,cs_delta)
    !call mean_scalar(cs2,cs2_ave,cs2_delta)
    call secondmomentum2(CF,xi2_ave,xi2_err) 
    !do i=1,N 
    !  mean_scalar(CF(i,1,:),CF_ave(i,k),CF_err(i,k))
    !end do
    write(10,*) T, E_ave/(vol), E_delta/(vol)
    write(20,*) T, M_ave/(vol), M_delta/(vol)
    write(30,*) T, suscep_ave/(vol), suscep_delta/(vol)
    write(40,*) T, heat_ave/(vol), heat_delta/(vol)
    write(50,*) T, U4_ave, U4_delta
    write(60,*) T, AR_ave,AR_err
    write(70,*) T, xi2_ave, xi2_err
    !write(70,*) T, cs2_ave/(vol**2), cs2_delta/(vol**2)
  end do
  
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  close(60)
  !close(70)
  deallocate(corr2)
  end subroutine vary_temp
  
  subroutine correlate(Ti,Tf,NTs)
  real(dp), intent(in) :: Ti,Tf
  integer(i4), intent(in) :: NTs
  integer(i4) :: i,j,k,i2
  integer(i4), allocatable :: spin(:,:)
  real(dp), allocatable :: corr1(:)
  real(dp), allocatable :: corr2(:,:)
  real(dp), allocatable :: CF(:,:),CF_ave(:,:),CF_err(:,:)
  real(dp) :: T,vol,norm
  open(60, file = 'data/corrfunc.dat', status = 'replace')
    vol=real(N**2,dp)
    norm=real(Nmsrs,dp)
    allocate(spin(N,N))
    allocate(corr1(N))
    allocate(corr2(N,N))
    allocate(CF(N,Nmsrs2))
    allocate(CF_ave(N,Nts))
    allocate(CF_err(N,Nts))
    do k=1,Nts
      T=Ti+(Tf-Ti)*real(k-1,dp)/real(Nts-1)
      write(*,*) k
      call cold_start(spin)
      do j=1,thermalization
        call montecarlo(spin,T)
        !call cluster(spint,T)
      end do
      do j=1,Nmsrs2
        call initialize(corr1,corr2)
        do i=1,Nmsrs
          do i2=1,eachsweep
            call montecarlo(spin,T)
            !call cluster(spint,T)          
          end do
          call correlation(spin,corr1,corr2)
        end do
        corr1=corr1/norm
        corr2=corr2/norm
        call correlation_function(corr1,corr2,CF(:,j))
      end do
      do j=1,N
        call mean_scalar(CF(j,:),CF_ave(j,k) ,CF_err(j,k))
      end do
    end do
    do k=1,N+1
      write(60,*) abs(k-1), CF_ave(iv(k),:), CF_err(iv(k),:)
    end do
    close(60)
    deallocate(spin)
    deallocate(corr1,corr2,CF,CF_ave,CF_err)
  end subroutine correlate

end program main
