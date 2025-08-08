module functions
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters, only : N, q
    implicit none

contains

  function iv(i)
    integer(i4), intent(in) :: i
    integer(i4) :: iv
    if(i==N+1) then
      iv=1
    else if(i==0) then
      iv=N
    else
      iv=i
    end if
  end function iv

  function Hamilt(spin)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp) :: Hamilt,neigh
    integer(i4) :: i,j
    Hamilt=0._dp
    do i=1,size(spin,dim=1)
      do j=1,size(spin,dim=2)
        !neigh=real(spin(iv(i+1),j)+spin(iv(i-1),j)+spin(i,iv(j+1))+spin(i,iv(j-1)),dp)
        !Hamilt=Hamilt-real(spin(i,j),dp)*neigh/2._dp
        neigh=real(spin(iv(i+1),j)+spin(i,iv(j+1)),dp)
        Hamilt=Hamilt-real(spin(i,j),dp)*neigh
      end do
    end do
  end function Hamilt

  function DeltaH(spin,i,j)
    integer(i4), dimension(:,:), intent(in) :: spin
    integer(i4),intent(in) :: i,j
    real(dp) :: DeltaH,neigh
    neigh=real(spin(iv(i+1),j)+spin(iv(i-1),j)+spin(i,iv(j+1))+spin(i,iv(j-1)),dp)
    DeltaH=2._dp*real(spin(i,j),dp)*neigh
  end function DeltaH
  
  function Hamiltfbc(spin)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp) :: Hamiltfbc,neigh
    integer(i4) :: Narr,i,j
    Narr=size(spin,dim=1)
    Hamiltfbc=0._dp
    do i=1,Narr
      do j=1,Narr
        if((i .ne. Narr) .and. (j==Narr)) then
          neigh=real(spin(i+1,j),dp)
        else if((i==Narr) .and. (j .ne. Narr)) then
          neigh=real(spin(i,j+1),dp)
        else if (i==Narr .and. j==Narr) then
          neigh=0._dp
        else
          neigh=real(spin(i+1,j)+spin(i,j+1),dp)
        end if
        Hamiltfbc=Hamiltfbc-real(spin(i,j),dp)*neigh
      end do
    end do
  end function Hamiltfbc

  function DeltaHfbc(spin,i,j)
    integer(i4), dimension(:,:), intent(in) :: spin
    integer(i4),intent(in) :: i,j
    real(dp) :: DeltaHfbc,neigh
    integer(i4) :: Narr
    Narr=size(spin,dim=1)    
    if( i==1 .and. j==1 ) then
      neigh=real(spin(i+1,j)+spin(i,j+1),dp)
    else if(i==1 .and. j==Narr) then
      neigh=real(spin(i+1,j)+spin(i,j-1),dp)
    else if(i==Narr .and. j==1) then
      neigh=real(spin(i-1,j)+spin(i,j+1),dp)
    else if(i==Narr .and. j==Narr) then
      neigh=real(spin(i-1,j)+spin(i,j-1),dp)
    else if(i==1 .and. (j .ne. 1) .and. (j .ne. Narr) ) then
      neigh=real(spin(i+1,j)+spin(i,j+1)+spin(i,j-1),dp)
    else if(i==Narr .and. (j .ne. 1) .and. (j .ne. Narr) ) then
      neigh=real(spin(i-1,j)+spin(i,j+1)+spin(i,j-1),dp)
    else if(j==1 .and. (i .ne. 1) .and. (i .ne. Narr) ) then
      neigh=real(spin(i+1,j)+spin(i-1,j)+spin(i,j+1),dp)
    else if(j==Narr .and. (i .ne. 1) .and. (i .ne. Narr) ) then
      neigh=real(spin(i+1,j)+spin(i-1,j)+spin(i,j-1),dp)
    else 
      neigh=real(spin(i+1,j)+spin(i-1,j)+spin(i,j+1)+spin(i,j-1),dp)
    end if
    DeltaHfbc=2._dp*real(spin(i,j),dp)*neigh
  end function DeltaHfbc

  function Magnet(spin)
    integer(i4), dimension(:,:), intent(in) :: spin
    integer(i4) :: i,j
    real(dp) :: Magnet
    Magnet=0._dp
    do i=1,size(spin,dim=1)
      do j=1,size(spin,dim=2)
        Magnet=Magnet+real(spin(i,j),dp)
      end do
    end do
  end function Magnet
  
  function qexp(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f
    real(dp) :: k1,k2
      k1=1._dp-q
      k2=1._dp+(1._dp -q)*x
      !if( k2 .le.  0._dp) then
      !  write(*,*) "Error"
      !  stop
      !end if
      f=k2**(1._dp/k1)
  end function qexp
  
  function p_metropolis(T,dH,E)
    real(dp), intent(in) :: T,dH,E
    real(dp) :: p_metropolis
    real(dp) :: T2,p1,p2
    T2=2.5_dp
    !p_metropolis=exp(-dH/T)
    !p_metropolis=qexp(-x)
    !p1=exp(-dH/T)/(1._dp+exp( (1._dp/T-1._dp/T2)*abs(E)))
    !p2=exp(-dH/T2)/(1._dp+exp( -(1._dp/T-1._dp/T2)*abs(E)))
    p1=exp(-dH/T)/(1._dp+exp( (1._dp/T-1._dp/T2)*E ))
    p2=exp(-dH/T2)/(1._dp+exp( -(1._dp/T-1._dp/T2)*E ))
    p_metropolis=min(1._dp,p1+p2)
  end function
  
  recursive function find(x,parent) result(out)
    integer(i4), intent(in) :: x
    integer(i4), intent(inout) :: parent(:)
    integer(i4) :: out
    if(parent(x) /= x) then
      out=find(parent(x),parent )
    else 
      out=x
    end if
  end function find

  subroutine union(x,y,parent)
    integer(i4),intent(in) :: x,y
    integer(i4),intent(inout) :: parent(:)
    integer :: root_x, root_y
    root_x=find(x,parent)
    root_y=find(y,parent)
    if (root_x /= root_y) then
      parent(root_y)=root_x
    end if
  end subroutine union


end module functions
