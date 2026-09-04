module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    
    integer(i4), parameter :: N=8,thermalization=20000,eachsweep=500,Nmsrs=500
    integer(i4),parameter :: Mbins=10,Nauto=15000,Nmsrs2=120,Mbin(5)=(/4,5,10,15,20/)
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real(dp), parameter :: q=1.1_dp,vol=real(N*N,dp)
    real(dp), parameter :: PI=4._dp*Atan(1.0_dp)
    real :: starting,ending

end module parameters
