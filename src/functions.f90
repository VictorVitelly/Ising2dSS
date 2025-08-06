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
  
  function p_metropolis(x)
    real(dp), intent(in) :: x    
    real(dp) :: p_metropolis
    p_metropolis=exp(x)
    !p_metropolis=qexp(-x)
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
