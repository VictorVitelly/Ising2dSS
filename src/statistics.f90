module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine cold_start(spin)
    integer(i4), dimension(:,:), intent(out) :: spin
    spin=-1
  end subroutine cold_start

  subroutine hot_start(spin)
    integer(i4), dimension(:,:), intent(out) :: spin
    integer(i4) :: i,j
    real(dp) :: r1
    do i=1,size(spin,dim=1)
      do j=1,size(spin,dim=2)
        call random_number(r1)
        if(r1 .le. 0.5_dp ) then
          spin(i,j)=1
        else
          spin(i,j)=-1
        end if
      end do
    end do
  end subroutine hot_start

  subroutine montecarlo(spin,T)
    integer(i4), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: T
    integer(i4) :: i1,i2
    real(dp) :: dH,r1,p
    do i1=1,size(spin,dim=1)
      do i2=1,size(spin,dim=2)
        dH=DeltaH(spin,i1,i2)
        if(dH .le. 0._dp) then
          spin(i1,i2)=-spin(i1,i2)
        else
          call random_number(r1)
          !p=(qexp(-dH/T))**q
          p=p_metropolis(T,dH)
          if(r1 < p ) then
            spin(i1,i2)=-spin(i1,i2)
          end if
        end if
      end do
    end do
  end subroutine montecarlo

  subroutine montecarlobien(spin,T,Energy)
    integer(i4), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: T
    real(dp), intent(inout) :: Energy
    integer(i4) :: i1,i2
    real(dp) :: dH,r1,p
    do i1=1,N
      do i2=1,N
        dH=DeltaH(spin,i1,i2)
        if(dH .le. 0._dp) then
          spin(i1,i2)=-spin(i1,i2)
          Energy=Energy+dH
        else
          call random_number(r1)
          !p=(qexp(-dH/(T-(1._dp-q)*Energy) ))**q
          p=(qexp(-(Energy+dH)/T)/qexp(-Energy/T))**q
          !p=min((qexp(-dH/(T-(1._dp-q)*Energy) ))**q,1._dp)
          !p=min((qexp(-(Energy+dH)/T)/qexp(-Energy/T))**q,1._dp)
            !p=(qexp(-dH/T))**q
          !write(*,*) Energy, p, qexp(-(Energy+dH)/T), qexp(-(Energy)/T)
          if(r1 < p) then
            spin(i1,i2)=-spin(i1,i2)
            Energy=Energy+dH
          end if
        end if
      end do
    end do
  end subroutine montecarlobien
  
  subroutine metropolisbien(spin,T,Energy,AR)
    integer(i4), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: T
    real(dp), intent(out) :: AR
    real(dp), intent(inout) :: Energy
    integer(i4) :: i1,i2
    real(dp) :: dH,r1,p
    AR=0._dp
    do i1=1,N
      do i2=1,N
        dH=DeltaH(spin,i1,i2)
        if(dH .le. 0._dp) then
          spin(i1,i2)=-spin(i1,i2)
          Energy=Energy+dH
          AR=AR+1._dp
        else
          call random_number(r1)
          !p=(qexp(-dH/(T-(1._dp-q)*Energy) ))**q
          p=(qexp(-(Energy+dH)/T)/qexp(-Energy/T))**q
          !p=min((qexp(-dH/(T-(1._dp-q)*Energy) ))**q,1._dp)
          !p=min((qexp(-(Energy+dH)/T)/qexp(-Energy/T))**q,1._dp)
            !p=(qexp(-dH/T))**q
          !write(*,*) Energy, p, qexp(-(Energy+dH)/T), qexp(-(Energy)/T)
          if(r1 < p) then
            spin(i1,i2)=-spin(i1,i2)
            Energy=Energy+dH
          end if
        AR=AR+p
        end if
      end do
    end do
    AR=AR/vol
  end subroutine metropolisbien

  subroutine flip_sign(spin,i)
    integer(i4), dimension(:,:), intent(inout) :: spin
    integer(i4), intent(in) :: i
    integer(i4) :: j
    if(mod(i,eachsweep)==0) then
      spin(:,:)=-spin(:,:)
    end if
  end subroutine flip_sign
  
  subroutine cluster(spin,T) 
    real(dp),intent(in) :: T 
    integer(i4), dimension(N,N),intent(inout) :: spin
    logical, dimension(N,N) :: bond_x,bond_y
    integer(i4) :: i,j,label(N,N),parent(N*N),next_label,left_label,up_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: r,p

    do i=1,N
      do j=1,N
        if(spin(i,j)==spin(mod(i,N)+1,j) ) then
          p=p_link(T)
          call random_number(r)
          bond_x(i,j)=(r<p)
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,N)+1) ) then
          p=p_link(T)
          call random_number(r)
          bond_y(i,j)=(r<p)
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,N*N
      parent(i)=i
    end do
    next_label=1
    left_label=0
    up_label=0

    do i=1,N
      do j=1,N
        left_label=0
        up_label=0
        if(i>1 .and. bond_x(i-1,j) ) then
          left_label=label(i-1,j)
        end if
        if(j>1 .and. bond_y(i,j-1) ) then
          up_label=label(i,j-1)
        end if
        if(left_label==0 .and. up_label==0) then
          label(i,j)=next_label
          next_label=next_label+1  
        else if(left_label /= 0 .and. up_label==0) then
          label(i,j)=left_label
        else if(left_label== 0 .and. up_label/=0) then
          label(i,j)=up_label
        else
          label(i,j)=min(left_label,up_label)
          call union(left_label,up_label,parent)
        end if
      end do
    end do
    
    do j=1,N
      if(bond_x(N,j) ) then
        call union(label(1,j),label(N,j),parent )
      end if
    end do
    
    do i=1,N
      if(bond_y(i,N) ) then
        call union(label(i,1),label(i,N),parent )
      end if 
    end do
    
    do i=1,N
      do j=1,N
        label(i,j)=find(label(i,j),parent)
      end do
    end do

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.

    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    do i=1,N
      do j=1,N
        if(flip_cluster(label(i,j))) then
          spin(i,j)=-spin(i,j)
        end if
      end do
    end do
    
  end subroutine cluster
  
!This cluster subroutine measures the size of the clusters
  subroutine cluster2(spin,T,rank_tot,rank2_tot) 
    real(dp),intent(in) :: T 
    integer(i4), dimension(N,N),intent(inout) :: spin
    real(dp), intent(out) :: rank_tot,rank2_tot
    logical, dimension(N,N) :: bond_x,bond_y
    integer(i4) :: i,j,label(N,N),parent(N*N),next_label,left_label,up_label
    logical, allocatable :: flip_cluster(:)
    integer(i4),allocatable :: rank(:)
    integer(i4) :: rank_ave,rank_n
    real(dp) :: r,p
    
    do i=1,N
      do j=1,N
        if(spin(i,j)==spin(mod(i,N)+1,j) ) then
          p=1._dp-exp(-2._dp/T )
          call random_number(r)
          bond_x(i,j)=(r<p)
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,N)+1) ) then
          p=1._dp-exp(-2._dp/T )
          call random_number(r)
          bond_y(i,j)=(r<p)
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,N*N
      parent(i)=i
    end do
    next_label=1
    left_label=0
    up_label=0

    do i=1,N
      do j=1,N
        left_label=0
        up_label=0
        if(i>1 .and. bond_x(i-1,j) ) then
          left_label=label(i-1,j)
        end if
        if(j>1 .and. bond_y(i,j-1) ) then
          up_label=label(i,j-1)
        end if
        if(left_label==0 .and. up_label==0) then
          label(i,j)=next_label
          next_label=next_label+1  
        else if(left_label /= 0 .and. up_label==0) then
          label(i,j)=left_label
        else if(left_label== 0 .and. up_label/=0) then
          label(i,j)=up_label
        else
          label(i,j)=min(left_label,up_label)
          call union(left_label,up_label,parent)
        end if
      end do
    end do
    
    do j=1,N
      if(bond_x(N,j) ) then
        call union(label(1,j),label(N,j),parent )
      end if
    end do
    
    do i=1,N
      if(bond_y(i,N) ) then
        call union(label(i,1),label(i,N),parent )
      end if 
    end do
    
    do i=1,N
      do j=1,N
        label(i,j)=find(label(i,j),parent)
      end do
    end do

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.
    
    allocate(rank(next_label) ) !cluster size
    rank(:)=0  !cluster size

    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    do i=1,N
      do j=1,N
        if(flip_cluster(label(i,j))) then
          spin(i,j)=-spin(i,j)
        end if
        rank(label(i,j))=rank(label(i,j))+1 !cluster size
      end do
    end do
    
    rank_ave=0 !cluster size
    rank_n=0 !cluster size
    rank2_tot=0._dp
    do i=1,next_label !cluster size
      if(rank(i).ne.0) then !cluster size
        rank_ave=rank_ave+rank(i) !cluster size
        rank_n=rank_n+1 !cluster size
        rank2_tot=rank2_tot+real(rank(i)**2,dp)
      end if !cluster size
    end do !cluster size
    rank_tot=real(rank_ave,dp)/real(rank_n,dp) !cluster size
    
    deallocate(flip_cluster,rank) !cluster size
    
  end subroutine cluster2

  
  !Error statistics

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i4) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar

  subroutine mean_vector(x,y,deltay)
    real(dp), dimension(N,Nmsrs), intent(in) :: x
    real(dp), dimension(N), intent(out) :: y,deltay
    integer(i4) :: i1
    y=0._dp
    deltay=0._dp
    do i1=1,N
      call mean_scalar(x(i1,:),y(i1),deltay(i1))
    end do
  end subroutine mean_vector

  subroutine mean_matrix(x,y,deltay)
    real(dp), dimension(N,N,Nmsrs), intent(in) :: x
    real(dp), dimension(N,N), intent(out) :: y,deltay
    integer(i4) :: i1,i2
    y=0._dp
    deltay=0._dp
    do i1=1,N
      do i2=1,N
      call mean_scalar(x(i1,i2,:),y(i1,i2),deltay(i1,i2))
      end do
    end do
  end subroutine mean_matrix

  subroutine heat_jackk(heat1,heat2,heat_ave,deltaheat)
    real(dp), dimension(:), intent(in) :: heat1, heat2
    real(dp), intent(out) :: heat_ave, deltaheat
    integer(i4) :: N,k,i
    real(dp) :: heat1t,heat2t,jackk,Ntot
    real(dp), dimension(Mbins) :: heatmean1,heatmean2,heat_avev
      N=size(heat1)
      Ntot=real(N,dp)-real(N,dp)/real(Mbins,dp)
      call mean_0(heat1,heat1t)
      call mean_0(heat2,heat2t)
      heat_ave=heat1t-heat2t**2
      heatmean1=0._dp
      heatmean2=0._dp
      do i=1,Mbins
        do k=1,N
          if(k .le. (i-1)*N/Mbins) then
            heatmean1(i)=heatmean1(i)+heat1(k)
            heatmean2(i)=heatmean2(i)+heat2(k)
          else if(k > i*N/Mbins) then
            heatmean1(i)=heatmean1(i)+heat1(k)
            heatmean2(i)=heatmean2(i)+heat2(k)
          end if
        end do
        heat_avev(i)=(heatmean1(i)/Ntot) -(heatmean2(i)/Ntot)**2
      end do
      do k=1,Mbins
        jackk=jackk+(heat_avev(k)-heat_ave )**2
      end do
      deltaheat=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp))
  end subroutine heat_jackk
  
  subroutine correlation(spin,corr1,corr2)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp), dimension(N,N), intent(inout) :: corr2
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp) :: spinvec(N)
    integer(i4) :: i1,i2
    real(dp) :: x 
    spinvec=0._dp
    x=real(N,dp)
    do i1=1,N
      do i2=1,N
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    spinvec(:)=spinvec(:)/x
    do i1=1,N
      corr1(i1)=corr1(i1)+(spinvec(i1))
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+spinvec(i1)*spinvec(i2)
      end do
    end do
  end subroutine correlation
  
  subroutine secondmomentum2(CF,xi2_ave,xi2_err)
  real(dp),dimension(N,N,Nmsrs2),intent(in) :: CF
  real(dp),intent(out) :: xi2_ave,xi2_err 
  integer(i4) :: i1,i2,i3
  real(dp) :: F1(Nmsrs2),F2(Nmsrs2),F12(Nmsrs2),F12_ave,F12_err
  F1(:)=0._dp
  F2(:)=0._dp
  do i1=1,Nmsrs2
    do i2=1,N
      do i3=1,N
        F1(i1)=F1(i1)+CF(i2,i3,i1)
        F2(i1)=F2(i1)+CF(i2,i3,i1)*COS(real(i2-1,dp)*2._dp*PI/real(N,dp)) &
        &*COS(real(i3-1,dp)*2._dp*PI/real(N,dp))
      end do
    end do
  end do
  do i1=1,Nmsrs2
    F12(i1)=F1(i1)/F2(i1)
    !write(*,*) F1(i1), F2(i1), F12(i1)
  end do
  call mean_scalar(F12,F12_ave,F12_err)
  xi2_ave=sqrt( (F12_ave -1._dp))/(2._dp*abs(SIN(PI/real(N,dp))) ) 
  xi2_err=F12_err/(4._dp*sqrt(F12_ave-1._dp)*abs(SIN(PI/real(N,dp))) )
  !write(*,*) xi2_ave,xi2_err
  end subroutine secondmomentum2

end module statistics
