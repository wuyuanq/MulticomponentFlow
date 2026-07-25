
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_model

    implicit none

    type :: model
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
    end type model

end module RST_model
