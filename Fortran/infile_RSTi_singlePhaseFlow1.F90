
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

program infile_RSTi_singlePhaseFlow

    use RST_model
    use RST_singlePhaseFlowDriver
    implicit none

    type(model) :: modelCase
    integer :: i
    real(kind=8) :: uconst
    real(kind=8), dimension(:), pointer :: initPres

    modelCase%phase = 'g'
    modelCase%Nc = 2
    modelCase%Temp = 480.0
    modelCase%Lx = 4.0
    modelCase%Ly = 4.0
    modelCase%timeEnd = 0.1*365*24*3600.0
    modelCase%nx = 40
    modelCase%ny = 40
    modelCase%nt = 0.1*365*24
    allocate(modelCase%xs(modelCase%nx+1))
    do i = 1, modelCase%nx+1
        modelCase%xs(i) = (i-1)*modelCase%Lx/modelCase%nx
    end do
    allocate(modelCase%ys(modelCase%ny+1))
    do i = 1, modelCase%ny+1
        modelCase%ys(i) = (i-1)*modelCase%Ly/modelCase%ny
    end do
    allocate(modelCase%ts(modelCase%nt+1))
    do i = 1, modelCase%nt+1
        modelCase%ts(i) = (i-1)*modelCase%timeEnd/modelCase%nt
    end do
    modelCase%gravX = 0.0
    modelCase%gravY = -9.807
    allocate(modelCase%Kxx(modelCase%nx,modelCase%ny))
    modelCase%Kxx(:,:) = 9.8692327*1.D-15
    allocate(modelCase%Kyy(modelCase%nx,modelCase%ny))
    modelCase%Kyy = modelCase%Kxx
    allocate(modelCase%poro(modelCase%nx,modelCase%ny))
    modelCase%poro(:,:) = 0.2
    allocate(modelCase%src(modelCase%Nc,modelCase%nx,modelCase%ny))
    modelCase%src(:,:,:) = 0.0
    allocate(modelCase%isDiriX(2,modelCase%ny))
    modelCase%isDiriX(:,:) = 0
    modelCase%isDiriX(2,modelCase%ny) = 1
    allocate(modelCase%isDiriY(modelCase%nx,2))
    modelCase%isDiriY(:,:) = 0
    modelCase%isDiriY(modelCase%nx,2) = 1
    allocate(modelCase%pBdryX(2,modelCase%ny))
    allocate(initPres(modelCase%ny))
    open(70,file='../InitialReservoirData/P1.txt')
    do i = 1, modelCase%ny
        read(70,*) initPres(i)
    end do
    modelCase%pBdryX(:,:) = 0.0
    modelCase%pBdryX(2,modelCase%ny) = initPres(modelCase%ny)
    allocate(modelCase%pBdryY(modelCase%nx,2))
    modelCase%pBdryY(:,:) = 0.0
    modelCase%pBdryY(modelCase%nx,2) = initPres(modelCase%ny)
    allocate(modelCase%pInit(modelCase%nx,modelCase%ny))
    do i = 1, modelCase%nx
        modelCase%pInit(i, 1:modelCase%ny) = initPres(1:modelCase%ny)
    end do
    allocate(modelCase%xBdryX(modelCase%Nc,2,modelCase%ny))
    modelCase%xBdryX(:,:,:) = 0.0
    modelCase%xBdryX(1,1,1) = 1.0
    allocate(modelCase%xBdryY(modelCase%Nc,modelCase%nx,2))
    modelCase%xBdryY(:,:,:) = 0.0
    modelCase%xBdryY(1,1,1) = 1.0
    allocate(modelCase%xInit(modelCase%Nc,modelCase%nx,modelCase%ny))
    modelCase%xInit(1,:,:) = 0.0
    modelCase%xInit(2,:,:) = 1.0
    uconst = 2*1.D-6
    allocate(modelCase%uBdryX(2,modelCase%ny))
    modelCase%uBdryX(:,:) = 0.0
    modelCase%uBdryX(1,1) = -uconst
    allocate(modelCase%uBdryY(modelCase%nx,2))
    modelCase%uBdryY(:,:) = 0.0
    modelCase%uBdryY(1,1) = -uconst
    allocate(modelCase%ct(modelCase%Nc))
    modelCase%ct(1) = 190.0
    modelCase%ct(2) = 370.0
    allocate(modelCase%cp(modelCase%Nc))
    modelCase%cp(1) = 4.6*1.D6
    modelCase%cp(2) = 4.2*1.D6
    allocate(modelCase%af(modelCase%Nc))
    modelCase%af(1) = 0.01
    modelCase%af(2) = 0.15
    allocate(modelCase%mw(modelCase%Nc))
    modelCase%mw(1) = 0.016
    modelCase%mw(2) = 0.044
    allocate(modelCase%cv(modelCase%Nc))
    modelCase%cv(1) = 0.0062
    modelCase%cv(2) = 0.0045
    allocate(modelCase%delta(modelCase%Nc,modelCase%Nc))
    modelCase%delta(:,:) = 0.0
    modelCase%delta(1,2) = 0.036
    modelCase%delta(2,1) = modelCase%delta(1,2)
    modelCase%soludoc = 'case1'

    call driver(modelCase)

    deallocate(initPres)

end program infile_RSTi_singlePhaseFlow
