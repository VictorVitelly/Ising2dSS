module arrays
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters, only : N
    implicit none

    integer(i4), allocatable :: spin(:,:)
    integer(i4),dimension(:),allocatable :: ip,im

contains

  subroutine init_vecs()
  integer(i4) :: i
  real(dp) :: x1,x2
  allocate(ip(N),im(N))
  do i=1,N-1
    ip(i)=i+1
  end do
  ip(N)=1
  do i=2,N
    im(i)=i-1
  end do
  im(1)=N
  end subroutine init_vecs

end module arrays
