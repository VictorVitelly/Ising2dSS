module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine initialize2(corr1,corr2)
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize2

  subroutine divideN(susc1,heat1)
    real(dp), intent(inout) :: susc1,heat1
    susc1=susc1/real(Nmsrs,dp)
    heat1=heat1/real(Nmsrs,dp)
  end subroutine divideN

  subroutine correlation(spin,k,corr1,corr2)
    integer(i4), dimension(N,N), intent(in) :: spin
    integer(i4), intent(in) :: k
    !real(dp) :: M
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
    real(dp), dimension(N) :: spinvec
    integer(i4) :: i1,i2
    !M=0._dp
    spinvec=0._dp
    do i1=1,N
      do i2=1,N
        !M=M+phi(i1,i2)
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    do i1=1,N
      corr1(i1,k)=spinvec(i1)
      do i2=1,N
        corr2(i1,i2,k)=spinvec(i1)*spinvec(i2)
      end do
    end do
  end subroutine correlation

  subroutine correlation_function(corr1,corr2,CF,CFprom)
    real(dp), dimension(N,Nmsrs), intent(in) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(in) :: corr2
    real(dp), dimension(N,N), intent(out) :: CF,CFprom
    real(dp), dimension(N) :: corr1prom,corr1delta
    real(dp), dimension(N,N) :: corr2prom,corr2delta
    integer(i4) :: i1,i2
    corr1prom=0._dp
    corr2prom=0._dp
    corr1delta=0._dp
    corr2delta=0._dp
    call mean_vector(corr1,corr1prom,corr1delta)
    call mean_matrix(corr2,corr2prom,corr2delta)
    do i1=1,N
      do i2=1,N
        CF(i1,i2)=corr2prom(i1,i2)-corr1prom(i1)*corr1prom(i2)
        CFprom(i1,i2)=Sqrt((corr2delta(i1,i2))**2+(corr1prom(i1)*corr1delta(i2))**2 +(corr1prom(i2)*corr1delta(i1) )**2)
      end do
    end do
  end subroutine correlation_function

  subroutine correlation_function2(corr1,corr2,CF,CFprom)
    real(dp), dimension(N,Nmsrs), intent(in) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(in) :: corr2
    real(dp), dimension(N,N), intent(out) :: CF,CFprom
    real(dp), dimension(N,N):: jackk
    real(dp), allocatable :: corr1m(:,:)
    real(dp), allocatable :: corr2m(:,:,:),CFm(:,:,:)
    real(dp), dimension(N) :: corr1prom,corr1delta
    real(dp), dimension(N,N) :: corr2prom,corr2delta
    integer(i4) :: i1,i2,i3
    corr1prom=0._dp
    corr2prom=0._dp
    corr1delta=0._dp
    corr2delta=0._dp
    call mean_vector(corr1,corr1prom,corr1delta)
    call mean_matrix(corr2,corr2prom,corr2delta)
    do i1=1,N
      do i2=1,N
        CF(i1,i2)=corr2prom(i1,i2)-corr1prom(i1)*corr1prom(i2)
        !CFprom(i1,i2)=Sqrt((corr2delta(i1,i2))**2+(corr1prom(i1)*corr1delta(i2))**2 +(corr1prom(i2)*corr1delta(i1) )**2)
      end do
    end do
    allocate(corr1m(N,Mbins) )
    allocate(corr2m(N,N,Mbins) )
    allocate(CFm(N,N,Mbins) )
    jackk=0._dp
    do i1=1,Mbins
    corr1m(:,i1)=0._dp
    corr2m(:,:,i1)=0._dp
      do i2=1,Nmsrs
        if(i2 .le. (i1-1)*Nmsrs/Mbins) then
          corr1m(:,i1)=corr1m(:,i1)+corr1(:,i2)
          corr2m(:,:,i1)=corr2m(:,:,i1)+corr2(:,:,i2)
        else if(i2 > i1*Nmsrs/Mbins) then
          corr1m(:,i1)=corr1m(:,i1)+corr1(:,i2)
          corr2m(:,:,i1)=corr2m(:,:,i1)+corr2(:,:,i2)
        end if
      end do
    end do
    corr1m=corr1m/(real(Nmsrs,dp) -real(Nmsrs/Mbins,dp))
    corr2m=corr2m/(real(Nmsrs,dp) -real(Nmsrs/Mbins,dp))
    do i3=1,Mbins
      do i1=1,N
        do i2=1,N
          CFm(i1,i2,i3)=corr2m(i1,i2,i3)-corr1m(i1,i3)*corr1m(i2,i3)
        end do
      end do
    end do
    do i1=1,Mbins
      jackk(:,:)=jackk(:,:)+(CF(:,:)-CFm(:,:,i1) )**2
    end do
    CFprom(:,:)=Sqrt(real(Mbins-1,dp)*jackk(:,:)/real(Mbins,dp))
  end subroutine correlation_function2
  
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
