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
    real(dp) :: dH,r1,p,Energy
    Energy=Hamilt(spin)
    do i1=1,size(spin,dim=1)
      do i2=1,size(spin,dim=2)
        dH=DeltaH(spin,i1,i2)
        if(dH .le. 0._dp) then
          spin(i1,i2)=-spin(i1,i2)
          Energy=Energy+dH
        else
          call random_number(r1)
          p=p_metropolis(T,dH,Energy)
          if(r1 < p ) then
            spin(i1,i2)=-spin(i1,i2)
            Energy=Energy+dH
          end if
        end if
      end do
    end do
  end subroutine montecarlo
  
  subroutine montecarlo2(spin,T)
    integer(i4), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: T
    integer(i4) :: i1,i2
    real(dp) :: dH,r1,p,Energy
    Energy=Hamilt(spin)
    do i1=1,size(spin,dim=1)
      do i2=1,size(spin,dim=2)
        dH=DeltaH(spin,i1,i2)
        call random_number(r1)
        if(r1 < p_metropolis(T,dH,Energy) ) then
          spin(i1,i2)=-spin(i1,i2)
          Energy=Energy+dH
        end if
      end do
    end do
  end subroutine montecarlo2
  
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
  
  subroutine binder(M2,M4,U,U_delta)
    real(dp), dimension(Nmsrs), intent(in) :: M2,M4
    real(dp), intent(out) :: U,U_delta
    real(dp) :: M2_ave,M4_ave,Ntot,jackk
    real(dp) :: M2_1(Mbins),M4_1(Mbins),U_ave(Mbins)
    integer(i4) :: i,k
    Ntot=real(Nmsrs,dp)-real(Nmsrs,dp)/real(Mbins,dp)
    call mean_0(M2,M2_ave)
    call mean_0(M4,M4_ave)
    U=1._dp -M4_ave/(3._dp*M2_ave**2)
    M2_1=0._dp
    M4_1=0._dp
    jackk=0._dp
    do i=1,Mbins
      do k=1,Nmsrs
        if(k .le. (i-1)*Nmsrs/Mbins) then
          M2_1(i)=M2_1(i)+M2(k)
          M4_1(i)=M4_1(i)+M4(k)
        else if(k > i*Nmsrs/Mbins) then
          M2_1(i)=M2_1(i)+M2(k)
          M4_1(i)=M4_1(i)+M4(k)        
        end if
      end do
    U_ave(i)=1._dp-M4_1(i)*Ntot/(3._dp*M2_1(i)**2)
    end do
    do k=1,Mbins
      jackk=jackk+(U_ave(k)-U )**2
    end do
    U_delta=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp))
  end subroutine binder

end module statistics
