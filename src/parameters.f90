module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=32,thermalization=2000,eachsweep=100,Nmsrs=200
    integer(i4),parameter :: Mbins=10,Nauto=15000,Nmsrs2=120,Mbin(4)=(/5,10,15,20/)
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real(dp), parameter :: q=2.0_dp
    real :: starting,ending

end module parameters
