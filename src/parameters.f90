module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=64,thermalization=1000,eachsweep=20,Nmsrs=100
    integer(i4),parameter :: Mbins=10,Nauto=15000,Nmsrs2=120,Mbin(4)=(/5,10,15,20/)
    real(dp), parameter :: q=1.2_dp
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real :: starting,ending

end module parameters
