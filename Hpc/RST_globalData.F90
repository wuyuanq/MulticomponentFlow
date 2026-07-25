
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_globalData

    implicit none

    ! the physical parameters
    real(kind=8), parameter :: R = 8.314

    ! the program parameters
    integer, parameter :: MAX_BUF = 1.D5
    integer, parameter :: MAX_COMMBUF = 1.D6

    ! the model parameters
    integer :: pnrows
    integer :: pncols
    character :: phase ! 'l' means liquid, 'g' means gas
    integer :: Nc
    real(kind=8) :: Temp
    real(kind=8) :: Lx
    real(kind=8) :: Ly
    real(kind=8) :: timeEnd
    integer :: nx
    integer :: ny
    integer :: nt
    real(kind=8) :: gravX
    real(kind=8) :: gravY
    real(kind=8), dimension(:), pointer :: xs
    real(kind=8), dimension(:), pointer :: ys
    real(kind=8), dimension(:), pointer :: ts
    real(kind=8), dimension(:,:), pointer :: Kxx
    real(kind=8), dimension(:,:), pointer :: Kyy
    real(kind=8), dimension(:,:), pointer :: poro
    real(kind=8), dimension(:,:,:), pointer :: src
    integer, dimension(:,:), pointer :: isDiriX
    integer, dimension(:,:), pointer :: isDiriY
    real(kind=8), dimension(:,:), pointer :: pBdryX
    real(kind=8), dimension(:,:), pointer :: pBdryY
    real(kind=8), dimension(:,:), pointer :: pInit
    real(kind=8), dimension(:,:,:), pointer :: xBdryX
    real(kind=8), dimension(:,:,:), pointer :: xBdryY
    real(kind=8), dimension(:,:,:), pointer :: xInit
    real(kind=8), dimension(:,:), pointer :: uBdryX
    real(kind=8), dimension(:,:), pointer :: uBdryY
    real(kind=8), dimension(:), pointer :: ct
    real(kind=8), dimension(:), pointer :: cp
    real(kind=8), dimension(:), pointer :: af
    real(kind=8), dimension(:), pointer :: mw
    real(kind=8), dimension(:), pointer :: cv
    real(kind=8), dimension(:,:), pointer :: delta
    character(len = 10) :: soludoc

    ! the global variables
    real(kind=8), dimension(:,:), pointer :: p
    real(kind=8), dimension(:,:), pointer :: ux
    real(kind=8), dimension(:,:), pointer :: uy
    real(kind=8), dimension(:,:), pointer :: lambdax
    real(kind=8), dimension(:,:), pointer :: lambday
    real(kind=8), dimension(:,:), pointer :: Kxxbar
    real(kind=8), dimension(:,:), pointer :: Kyybar
    real(kind=8), dimension(:,:), pointer :: rho
    real(kind=8), dimension(:,:), pointer :: rhobarx
    real(kind=8), dimension(:,:), pointer :: rhobary
    real(kind=8), dimension(:,:,:), pointer :: x
    real(kind=8), dimension(:,:,:), pointer :: xbarx
    real(kind=8), dimension(:,:,:), pointer :: xbary
    real(kind=8), dimension(:,:), pointer :: xi
    real(kind=8), dimension(:,:), pointer :: xibarx
    real(kind=8), dimension(:,:), pointer :: xibary
    real(kind=8), dimension(:,:), pointer :: visc
    real(kind=8), dimension(:,:,:), pointer :: moleincell
    real(kind=8), dimension(:,:,:), pointer :: moleincell_old
    real(kind=8), dimension(:,:), pointer :: deri_xi_p
    real(kind=8), dimension(:,:,:), pointer :: deri_xi_n
    real(kind=8) :: totalmole
    integer :: t

    real(kind=8), dimension(:), pointer :: initial_x_guess
    real(kind=8) :: timestart
    integer :: num_procs, myid
    integer :: prow, pcol
    integer :: localnrows, localncols
    integer :: xlower, xupper, ylower, yupper

end module RST_globalData

