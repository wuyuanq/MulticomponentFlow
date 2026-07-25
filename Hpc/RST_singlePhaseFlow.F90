
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_singlePhaseFlow

    use RST_model
    use RST_globalData
    use RST_PREOS
    use RST_viscosity
    include 'mpif.h'

contains

    subroutine initialize(modelCase)

        implicit none

        type(model), intent(in out) :: modelCase

        logical :: alive
        character(len=50) :: fmhtxt, fmrtxt
        integer :: indexl, indexr, indexu, indexd
        integer :: i, j, k, m, ierr

        call MPI_INIT(ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()

        pnrows = modelCase%pnrows
        pncols = modelCase%pncols
        phase = modelCase%phase
        Nc = modelCase%Nc
        Temp = modelCase%Temp
        Lx = modelCase%Lx
        Ly = modelCase%Ly
        timeEnd = modelCase%timeEnd
        nx = modelCase%nx
        ny = modelCase%ny
        nt = modelCase%nt
        gravX = modelCase%gravX
        gravY = modelCase%gravY

        allocate(xs(nx+1))
        allocate(ys(ny+1))
        allocate(ts(nt+1))
        allocate(Kxx(nx,ny))
        allocate(Kyy(nx,ny))
        allocate(poro(nx,ny))
        allocate(src(Nc,nx,ny))
        allocate(isDiriX(2,ny))
        allocate(isDiriY(nx,2))
        allocate(pBdryX(2,ny))
        allocate(pBdryY(nx,2))
        allocate(pInit(nx,ny))
        allocate(xBdryX(Nc,2,ny))
        allocate(xBdryY(Nc,nx,2))
        allocate(xInit(Nc,nx,ny))
        allocate(uBdryX(2,ny))
        allocate(uBdryY(nx,2))
        allocate(ct(Nc))
        allocate(cp(Nc))
        allocate(af(Nc))
        allocate(mw(Nc))
        allocate(cv(Nc))
        allocate(delta(Nc,Nc))

        xs = modelCase%xs
        ys = modelCase%ys
        ts = modelCase%ts
        Kxx = modelCase%Kxx
        Kyy = modelCase%Kyy
        poro = modelCase%poro
        src = modelCase%src
        isDiriX = modelCase%isDiriX
        isDiriY = modelCase%isDiriY
        pBdryX = modelCase%pBdryX
        pBdryY = modelCase%pBdryY
        pInit = modelCase%pInit
        uBdryX = modelCase%uBdryX
        uBdryY = modelCase%uBdryY
        xBdryX = modelCase%xBdryX
        xBdryY = modelCase%xBdryY
        xInit = modelCase%xInit
        uBdryX = modelCase%uBdryX
        uBdryY = modelCase%uBdryY
        ct = modelCase%ct
        cp = modelCase%cp
        af = modelCase%af
        mw = modelCase%mw
        cv = modelCase%cv
        delta = modelCase%delta
        soludoc = modelCase%soludoc

        uBdryX(1, 1:ny) = -uBdryX(1, 1:ny)
        uBdryY(1:nx, 1) = -uBdryY(1:nx, 1)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

        prow = mod(myid,pnrows)+1
        pcol = myid/pnrows+1

        localnrows = ny/pnrows
        localncols = nx/pncols

        xlower = (pcol-1)*localncols+1
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows

        allocate(p(0:localnrows+1, 0:localncols+1))
        allocate(ux(localnrows, localncols+1))
        ux = 0
        allocate(uy(localnrows+1, localncols))
        uy = 0
        allocate(lambdax(localnrows, localncols+1))
        lambdax(:,:) = 0
        allocate(lambday(localnrows+1, localncols))
        lambday(:,:) = 0
        allocate(Kxxbar(localnrows, localncols+1))
        allocate(Kyybar(localnrows+1, localncols))
        allocate(rho(0:localnrows+1, 0:localncols+1))
        allocate(rhobarx(localnrows, localncols+1))
        rhobarx(:,:) = 0
        allocate(rhobary(localnrows+1, localncols))
        rhobary(:,:) = 0
        allocate(x(Nc, 0:localnrows+1, 0:localncols+1))
        allocate(xbarx(Nc, localnrows, localncols+1))
        xbarx(:,:,:) = 0
        allocate(xbary(Nc, localnrows+1, localncols))
        xbary(:,:,:) = 0
        allocate(xi(0:localnrows+1, 0:localncols+1))
        allocate(xibarx(localnrows, localncols+1))
        xibarx(:,:) = 0
        allocate(xibary(localnrows+1, localncols))
        xibary(:,:) = 0
        allocate(visc(0:localnrows+1, 0:localncols+1))
        allocate(moleincell(Nc, 0:localnrows+1, 0:localncols+1))
        allocate(moleincell_old(Nc, 0:localnrows+1, 0:localncols+1))
        allocate(deri_xi_p(0:localnrows+1, 0:localncols+1))
        allocate(deri_xi_n(Nc, 0:localnrows+1, 0:localncols+1))

        if(pcol == 1) then
            ux(1:localnrows, 1) = uBdryX(1, ylower:ylower+localnrows-1)
        end if
        if(pcol == pncols) then
            ux(1:localnrows, localncols+1) = uBdryX(2, ylower:ylower+localnrows-1)
        end if
        if(prow == 1) then
            uy(1, 1:localncols) = uBdryY(xlower:xlower+localncols-1, 1)
        end if
        if(prow == pnrows) then
            uy(localnrows+1, 1:localncols) = uBdryY(xlower:xlower+localncols-1, 2)
        end if

        indexd = 0
        indexu = localnrows+1
        indexl = 0
        indexr = localncols+1
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if

        do i = indexd, indexu
            do j = indexl, indexr
                p(i,j) = pInit(xlower+j-1,ylower+i-1)
            end do
        end do
        if(prow == 1) then
            do i = 1, localncols
                p(0, i) = pBdryY(xlower+i-1, 1)
            end do
        end if
        if(prow == pnrows) then
            do i = 1, localncols
                p(localnrows+1, i) = pBdryY(xlower+i-1, 2)
            end do
        end if
        if(pcol == 1) then
            do i = 1, localnrows
                p(i, 0) = pBdryX(1, ylower+i-1)
            end do
        end if
        if(pcol == pncols) then
            do i = 1, localnrows
                p(i, localncols+1) = pBdryX(2, ylower+i-1)
            end do
        end if

        if(pncols == 1) then ! only one process column
            do i = 1, localnrows
                Kxxbar(i,1) = Kxx(1,ylower+i-1)
            end do
            do i = 1, localnrows
                Kxxbar(i,localncols+1) = Kxx(nx,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 2, localncols
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(pcol == 1) then
            do i = 1, localnrows
                Kxxbar(i,1) = Kxx(1,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 2, localncols+1
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(pcol == pncols) then
            do i = 1, localnrows
                Kxxbar(i,localncols+1) = Kxx(nx,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 1, localncols
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        else
            do i = 1, localnrows
                do j = 1, localncols+1
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        end if

        if(pnrows == 1) then ! only one process row
            do i = 1, localncols
                Kyybar(1,i) = Kyy(xlower+i-1,1)
            end do
            do i = 1, localncols
                Kyybar(localnrows+1,i) = Kyy(xlower+i-1, ny)
            end do
            do i = 2, localnrows
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(prow == 1) then
            do i = 1, localncols
                Kyybar(1,i) = Kyy(xlower+i-1,1)
            end do
            do i = 2, localnrows+1
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(prow == pnrows) then
            do i = 1, localncols
                Kyybar(localnrows+1,i) = Kyy(xlower+i-1, ny)
            end do
            do i = 1, localnrows
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        else
            do i = 1, localnrows+1
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        end if

        do m = 1, Nc
            do i = indexd, indexu
                do j = indexl, indexr
                    x(m,i,j) = xInit(m,xlower+j-1,ylower+i-1)
                end do
            end do
        end do

        totalmole = 0.0
        t = 2

        allocate(initial_x_guess(nx*ny/num_procs))
        k = 1
        do i = 1, localnrows
            do j = 1, localncols
                initial_x_guess(k) = p(i,j)
                k = k + 1
            end do
        end do

        if(myid == 0) then
            inquire(file = soludoc, exist = alive)
            if(.not.alive) then
                call system('mkdir '//trim(adjustl(soludoc)))
            end if

            fmhtxt = trim(adjustl(soludoc))//'/soln_1PhFlw_moleHistory.txt'
            fmrtxt = trim(adjustl(soludoc))//'/soln_1PhFlw_moleRatioHistory.txt'

            open(unit=40, file=trim(adjustl(fmhtxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            end if
            open(unit=50, file=trim(adjustl(fmrtxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            end if
        end if

        deallocate(modelCase%xs)
        deallocate(modelCase%ys)
        deallocate(modelCase%ts)
        deallocate(modelCase%Kxx)
        deallocate(modelCase%Kyy)
        deallocate(modelCase%poro)
        deallocate(modelCase%src)
        deallocate(modelCase%isDiriX)
        deallocate(modelCase%isDiriY)
        deallocate(modelCase%pBdryX)
        deallocate(modelCase%pBdryY)
        deallocate(modelCase%pInit)
        deallocate(modelCase%xBdryX)
        deallocate(modelCase%xBdryY)
        deallocate(modelCase%xInit)
        deallocate(modelCase%uBdryX)
        deallocate(modelCase%uBdryY)
        deallocate(modelCase%ct)
        deallocate(modelCase%cp)
        deallocate(modelCase%af)
        deallocate(modelCase%mw)
        deallocate(modelCase%cv)
        deallocate(modelCase%delta)

    end subroutine initialize

    subroutine computeParameters()

        implicit none

        real(kind=8), dimension(:), pointer :: xtemp, deri_xi_ntemp, moleincelltemp
        real(kind=8) :: visctemp, deri_xi_ptemp
        real(kind=8), dimension(:), pointer :: leftmole
        real(kind=8) :: sum, leftmole1sum, totaldesiredleftmole, sentbuffer(2)
        integer :: indexd, indexu, indexl, indexr, prociteration
        logical :: isZero
        integer :: i, j, m, k, ierr

        integer :: status(MPI_STATUS_SIZE)
        integer :: request, position
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)

        allocate(xtemp(Nc))
        allocate(deri_xi_ntemp(Nc))
        allocate(moleincelltemp(Nc))

        indexd = 0
        indexu = localnrows+1
        indexl = 0
        indexr = localncols+1
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if

        do j = indexl, indexr
            do i = indexd, indexu
                if((.not.((i==0).and.(j==0))).and.(.not.((i==0).and.(j==localncols+1))).and.(.not. &!
                    ((i==localnrows+1).and.(j==0))).and.(.not.((i==localnrows+1).and.(j==localncols+1)))) then
                    if(t > 2) then
                        moleincell_old(1:Nc,i,j) = moleincell(1:Nc,i,j)
                    end if
                    xtemp(1:Nc) = x(1:Nc,i,j)
                    call PREOS( xtemp, p(i,j), xi(i,j), rho(i,j), deri_xi_p(i,j), deri_xi_ntemp, moleincelltemp )
                    deri_xi_n(1:Nc,i,j) = deri_xi_ntemp(1:Nc)
                    moleincell(1:Nc,i,j) = moleincelltemp(1:Nc)
                    if(t == 2) then
                        moleincell_old(1:Nc,i,j) = moleincell(1:Nc,i,j)
                    end if
                    visc(i,j) = viscosity( xtemp, xi(i,j), p(i,j) )

                end if
            end do
        end do

        allocate(leftmole(Nc))
        leftmole = 0
        do m = 1, Nc
            do j = 1, localncols
                do i = 1, localnrows
                    leftmole(m) = leftmole(m) + moleincell(m,i,j)
                end do
            end do
        end do

        totaldesiredleftmole = 0.0
        do m = 2, Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m)
        end do

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(myid /= 0) then

            position = 0
            call MPI_PACK(totaldesiredleftmole, 1, MPI_DOUBLE_PRECISION, sentbuffer, 2*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(leftmole(1), 1, MPI_DOUBLE_PRECISION, sentbuffer, 2*8, position, MPI_COMM_WORLD, ierr)
            call MPI_IBSEND(sentbuffer,2,MPI_DOUBLE_PRECISION,0,myid+num_procs,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)

        else

            sum = totaldesiredleftmole
            leftmole1sum = leftmole(1)

            do prociteration = 1, num_procs-1

                call MPI_RECV(sentbuffer, 2, MPI_DOUBLE_PRECISION, prociteration, prociteration+num_procs, &!
                    MPI_COMM_WORLD, status, ierr)

                sum = sum + sentbuffer(1)
                leftmole1sum = leftmole1sum + sentbuffer(2)

            end do

            if(t == 2) then
                totalmole = sum
            end if

            write(40, fmt='(es12.5)', iostat=ierr) (totalmole-sum)/totalmole
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            end if
            write(50, fmt='(es12.5)', iostat=ierr) sum/leftmole1sum
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            end if

        end if

        deallocate(leftmole)

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        if(pncols == 1) then
            do i = 1, localnrows
                isZero = .true.
                do m = 1, Nc
                    if(xBdryX(m,1,ylower+i-1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((ux(i,1)>0).and.(p(i,1)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryX(1:Nc,1,ylower+i-1)
                    call PREOS( xtemp, p(i,1), xibarx(i,1), rhobarx(i,1), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibarx(i,1), p(i,1) )
                    lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visctemp
                    xbarx(1:Nc,i,1) = xtemp(1:Nc)
                else
                    xibarx(i,1) = xi(i,1)
                    rhobarx(i,1) = rho(i,1)
                    lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visc(i,1)
                    xbarx(1:Nc,i,1) = x(1:Nc,i,1)
                end if
            end do
            do i = 1, localnrows


if(ux(i,localncols+1)<0) then
print *, 'a'
!pause
ux(i,localncols+1) = 0
end if


                isZero = .true.
                do m = 1, Nc
                    if(xBdryX(m,2,ylower+i-1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((ux(i,localncols+1)<0).and.(p(i,nx)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryX(1:Nc,2,ylower+i-1)
                    call PREOS( xtemp, p(i,nx), xibarx(i,localncols+1), rhobarx(i,localncols+1), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibarx(i,localncols+1), p(i,nx) )
                    lambdax(i,localncols+1) = Kxxbar(i,localncols+1)*xibarx(i,localncols+1)/visctemp
                    xbarx(1:Nc,i,localncols+1) = xtemp(1:Nc)
                else
                    xibarx(i,localncols+1) = xi(i,localncols)
                    rhobarx(i,localncols+1) = rho(i,localncols)
                    lambdax(i,localncols+1) = Kxxbar(i,localncols+1)*xibarx(i,localncols+1)/visc(i,localncols)
                    xbarx(1:Nc,i,localncols+1) = x(1:Nc,i,localncols)
                end if
            end do
            do j = 2, localncols
                do i = 1, localnrows
                    if(ux(i,j) > 0) then
                        xibarx(i,j) = xi(i,j-1)
                        rhobarx(i,j) = rho(i,j-1)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j-1)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j-1)
                    else
                        xibarx(i,j) = xi(i,j)
                        rhobarx(i,j) = rho(i,j)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(pcol == 1) then
            do i = 1, localnrows
                isZero = .true.
                do m = 1, Nc
                    if(xBdryX(m,1,ylower+i-1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((ux(i,1)>0).and.(p(i,1)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryX(1:Nc,1,ylower+i-1)
                    call PREOS( xtemp, p(i,1), xibarx(i,1), rhobarx(i,1), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibarx(i,1), p(i,1) )
                    lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visctemp
                    xbarx(1:Nc,i,1) = xtemp(1:Nc)
                else
                    xibarx(i,1) = xi(i,1)
                    rhobarx(i,1) = rho(i,1)
                    lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visc(i,1)
                    xbarx(1:Nc,i,1) = x(1:Nc,i,1)
                end if
            end do
            do i = 1, localnrows
                do j = 2, localncols+1
                    if(ux(i,j) > 0) then
                        xibarx(i,j) = xi(i,j-1)
                        rhobarx(i,j) = rho(i,j-1)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j-1)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j-1)
                    else
                        xibarx(i,j) = xi(i,j)
                        rhobarx(i,j) = rho(i,j)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(pcol==pncols) then

if(ux(i,localncols+1)<0) then
print *, 'b'
!pause
ux(i,localncols+1) = 0
end if

            do i = 1, localnrows
                isZero = .true.
                do m = 1, Nc
                    if(xBdryX(m,2,ylower+i-1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((ux(i,localncols+1)<0).and.(p(i,nx)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryX(1:Nc,2,ylower+i-1)
                    call PREOS( xtemp, p(i,nx), xibarx(i,localncols+1), rhobarx(i,localncols+1), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibarx(i,localncols+1), p(i,nx) )
                    lambdax(i,localncols+1) = Kxxbar(i,localncols+1)*xibarx(i,localncols+1)/visctemp
                    xbarx(1:Nc,i,localncols+1) = xtemp(1:Nc)
                else
                    xibarx(i,localncols+1) = xi(i,localncols)
                    rhobarx(i,localncols+1) = rho(i,localncols)
                    lambdax(i,localncols+1) = Kxxbar(i,localncols+1)*xibarx(i,localncols+1)/visc(i,localncols)
                    xbarx(1:Nc,i,localncols+1) = x(1:Nc,i,localncols)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    if(ux(i,j) > 0) then
                        xibarx(i,j) = xi(i,j-1)
                        rhobarx(i,j) = rho(i,j-1)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j-1)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j-1)
                    else
                        xibarx(i,j) = xi(i,j)
                        rhobarx(i,j) = rho(i,j)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        else
            do j = 1, localncols+1
                do i = 1, localnrows
                    if(ux(i,j) > 0) then
                        xibarx(i,j) = xi(i,j-1)
                        rhobarx(i,j) = rho(i,j-1)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j-1)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j-1)
                    else
                        xibarx(i,j) = xi(i,j)
                        rhobarx(i,j) = rho(i,j)
                        lambdax(i,j) = Kxxbar(i,j)*xibarx(i,j)/visc(i,j)
                        xbarx(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        end if

        if(pnrows == 1) then
            do i = 1, localncols

if(uy(1,i)>0) then
print *, 'c'
!pause
uy(1,i) = 0
end if

                isZero = .true.
                do m = 1, Nc
                    if(xBdryY(m,xlower+i-1,1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((uy(1,i)>0).and.(p(1,i)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryY(1:Nc,xlower+i-1,1)
                    call PREOS( xtemp, p(1,i), xibary(1,i), rhobary(1,i), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibary(1,i), p(1,i) )
                    lambday(1,i) = Kyybar(1,i)*xibary(1,i)/visctemp
                    xbary(1:Nc,1,i) = xtemp(1:Nc)
                else
                    xibary(1,i) = xi(1,i)
                    rhobary(1,i) = rho(1,i)
                    lambday(1,i) = Kyybar(1,i)*xibary(1,i)/visc(1,i)
                    xbary(1:Nc,1,i) = x(1:Nc,1,i)
                end if
            end do
            do i = 1, localncols
                isZero = .true.
                do m = 1, Nc
                    if(xBdryY(m,i,2) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((uy(localnrows+1,i)<0).and.(p(localnrows+1,i+1)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryY(1:Nc,i,2)
                    call PREOS( xtemp, p(localnrows+1,i+1), xibary(localnrows+1,i), rhobary(localnrows+1,i), &!
                        deri_xi_ptemp, deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibary(localnrows+1,i), p(localnrows+1,i+1) )
                    lambday(localnrows+1,i) = Kyybar(localnrows+1,i)*xibary(localnrows+1,i)/visctemp
                    xbary(1:Nc,localnrows+1,i) = xtemp(1:Nc)
                else
                    xibary(localnrows+1,i) = xi(localnrows,i)
                    rhobary(localnrows+1,i) = rho(localnrows,i)
                    lambday(localnrows+1,i) = Kyybar(localnrows+1,i)*xibary(localnrows+1,i)/visc(localnrows,i)
                    xbary(1:Nc,localnrows+1,i) = x(1:Nc,localnrows,i)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows
                    if(uy(i,j) > 0) then
                        xibary(i,j) = xi(i-1,j)
                        rhobary(i,j) = rho(i-1,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i-1,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i-1,j)
                    else
                        xibary(i,j) = xi(i,j)
                        rhobary(i,j) = rho(i,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(prow==1) then
            do i = 1, localncols

if(uy(1,i)>0) then
print *, 'd'
uy(1,i) = 0
!pause
end if

                isZero = .true.
                do m = 1, Nc
                    if(xBdryY(m,xlower+i-1,1) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((uy(1,i)>0).and.(p(1,i)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryY(1:Nc,xlower+i-1,1)
                    call PREOS( xtemp, p(1,i), xibary(1,i), rhobary(1,i), deri_xi_ptemp, &!
                        deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibary(1,i), p(1,i) )
                    lambday(1,i) = Kyybar(1,i)*xibary(1,i)/visctemp
                    xbary(1:Nc,1,i) = xtemp(1:Nc)
                else
                    xibary(1,i) = xi(1,i)
                    rhobary(1,i) = rho(1,i)
                    lambday(1,i) = Kyybar(1,i)*xibary(1,i)/visc(1,i)
                    xbary(1:Nc,1,i) = x(1:Nc,1,i)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows+1
                    if(uy(i,j) > 0) then
                        xibary(i,j) = xi(i-1,j)
                        rhobary(i,j) = rho(i-1,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i-1,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i-1,j)
                    else
                        xibary(i,j) = xi(i,j)
                        rhobary(i,j) = rho(i,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(prow==pnrows) then
            do i = 1, localncols
                isZero = .true.
                do m = 1, Nc
                    if(xBdryY(m,i,2) /= 0) then
                        isZero = .false.
                        exit
                    end if
                end do
                if((uy(localnrows+1,i)<0).and.(p(localnrows+1,i+1)/=0).and.(.not.isZero)) then
                    xtemp(1:Nc) = xBdryY(1:Nc,i,2)
                    call PREOS( xtemp, p(localnrows+1,i+1), xibary(localnrows+1,i), rhobary(localnrows+1,i), &!
                        deri_xi_ptemp, deri_xi_ntemp, moleincelltemp )
                    visctemp = viscosity( xtemp, xibary(localnrows+1,i), p(localnrows+1,i+1) )
                    lambday(localnrows+1,i) = Kyybar(localnrows+1,i)*xibary(localnrows+1,i)/visctemp
                    xbary(1:Nc,localnrows+1,i) = xtemp(1:Nc)
                else
                    xibary(localnrows+1,i) = xi(localnrows,i)
                    rhobary(localnrows+1,i) = rho(localnrows,i)
                    lambday(localnrows+1,i) = Kyybar(localnrows+1,i)*xibary(localnrows+1,i)/visc(localnrows,i)
                    xbary(1:Nc,localnrows+1,i) = x(1:Nc,localnrows,i)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    if(uy(i,j) > 0) then
                        xibary(i,j) = xi(i-1,j)
                        rhobary(i,j) = rho(i-1,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i-1,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i-1,j)
                    else
                        xibary(i,j) = xi(i,j)
                        rhobary(i,j) = rho(i,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        else
            do j = 1, localncols
                do i = 1, localnrows+1
                    if(uy(i,j) > 0) then
                        xibary(i,j) = xi(i-1,j)
                        rhobary(i,j) = rho(i-1,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i-1,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i-1,j)
                    else
                        xibary(i,j) = xi(i,j)
                        rhobary(i,j) = rho(i,j)
                        lambday(i,j) = Kyybar(i,j)*xibary(i,j)/visc(i,j)
                        xbary(1:Nc,i,j) = x(1:Nc,i,j)
                    end if
                end do
            end do
        end if

        deallocate(xtemp)
        deallocate(deri_xi_ntemp)
        deallocate(moleincelltemp)

    end subroutine computeParameters

    subroutine computePres()

        implicit none

        real(kind=8) :: xedge, yedge, ledge, redge, uedge, dedge
        real(kind=8) :: up, down, left, right
        real(kind=8), dimension(:), pointer :: psent
        real(kind=8), dimension(:), pointer :: recvbuffer
        integer :: recvbuffersize
        integer :: i, j, r, m, k

        integer :: ierr, num_iter, errorcode
        integer :: status(MPI_STATUS_SIZE)
        integer :: request, requestl, requestr, requestu, requestd
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)

        integer :: local_size
        integer(kind=8) :: grid
        integer(kind=8) :: stencil
        integer(kind=8) :: A
        integer(kind=8) :: b
        integer(kind=8) :: x
        integer(kind=8) :: solver
        integer :: ilower(2), iupper(2)
        integer :: myentry
        integer :: offsets(5,2) = reshape([0,-1,1,0,0,0,0,0,-1,1], [5,2])
        integer :: stencil_indices(5) = (/0, 1, 2, 3, 4/)
        real(kind=8), dimension(:), pointer :: values
        real(kind=8), dimension(:), pointer :: rhs_values, x_values

        local_size = nx*ny/num_procs

        allocate(rhs_values(local_size))
        allocate(x_values(local_size))
        allocate(values(local_size*5))
        values = 0

        ilower(1) = xlower
        iupper(1) = xupper
        ilower(2) = ylower
        iupper(2) = yupper

        call HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, grid, ierr)
        call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)
        call HYPRE_StructGridAssemble(grid, ierr)

        call HYPRE_StructStencilCreate(2, 5, stencil, ierr)
        do myentry = 1, 5
            call HYPRE_StructStencilSetElement(stencil, myentry-1, offsets(myentry,1:2), ierr)
        end do

        call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, A, ierr)
        call HYPRE_StructMatrixInitialize(A, ierr)

        k = 1
        r = 1
        do i = 1, localnrows
            do j = 1, localncols
            
                xedge = xs(xlower+j) - xs(xlower+j-1)
                yedge = ys(ylower+i) - ys(ylower+i-1)

                if((pcol==1).and.(j==1)) then
                    ledge = 0
                else
                    ledge = xs(xlower+j-1) - xs(xlower+j-2)
                end if

                if((pcol==pncols).and.(j==localncols)) then
                    redge = 0
                else
                    redge = xs(xlower+j+1) - xs(xlower+j)
                end if

                if((prow==1).and.(i==1)) then
                    dedge = 0
                else
                    dedge = ys(ylower+i-1) - ys(ylower+i-2)
                end if

                if((prow==pnrows).and.(i==localnrows)) then
                    uedge = 0
                else
                    uedge = ys(ylower+i+1) - ys(ylower+i)
                end if

                values(r) = poro(xlower+j-1,ylower+i-1)*deri_xi_p(i,j)/(timeEnd/nt)

                rhs_values(k) = poro(xlower+j-1,ylower+i-1)*deri_xi_p(i,j)/(timeEnd/nt)*p(i,j)
                do m = 1, Nc
                    rhs_values(k) = rhs_values(k) + src(m,xlower+j-1,ylower+i-1) - &!
                        poro(xlower+j-1,ylower+i-1)*deri_xi_n(m,i,j)* &!
                        (moleincell(m,i,j)-moleincell_old(m,i,j))/(timeEnd/nt)
                end do

                up = -2*lambday(i+1,j)/yedge/(yedge+uedge)
                down = -2*lambday(i,j)/yedge/(yedge+dedge)
                left = -2*lambdax(i,j)/xedge/(xedge+ledge)
                right = -2*lambdax(i,j+1)/xedge/(xedge+redge)

                if((prow == 1).and.(i == 1).and.(isDiriY(xlower+j-1,1) == 0)) then
                    rhs_values(k) = rhs_values(k) + uBdryY(xlower+j-1,1)*xibary(i,j)/yedge
                elseif((prow == 1).and.(i == 1).and.(isDiriY(xlower+j-1,1) == 1)) then
                    rhs_values(k) = rhs_values(k) - down*p(i-1,j) + lambday(i,j)*rhobary(i,j)*gravY/yedge
                    values(r) = values(r) - down
                else
                    values(r+3) = down
                    values(r) = values(r) - down
                    rhs_values(k) = rhs_values(k) + lambday(i,j)*rhobary(i,j)*gravY/yedge
                end if

                if((pcol == 1).and.(j == 1).and.(isDiriX(1,ylower+i-1) == 0)) then
                    rhs_values(k) = rhs_values(k) + uBdryX(1,ylower+i-1)*xibarx(i,j)/xedge
                else if((pcol == 1).and.(j == 1).and.(isDiriX(1,ylower+i-1) == 1)) then
                    rhs_values(k) = rhs_values(k) - left*p(i,j-1) + lambdax(i,j)*rhobarx(i,j)*gravX/xedge
                    values(r) = values(r) - left
                else
                    values(r+1) = left
                    values(r) = values(r) - left
                    rhs_values(k) = rhs_values(k) + lambdax(i,j)*rhobarx(i,j)*gravX/xedge
                end if

                if((pcol == pncols).and.(j == localncols).and.(isDiriX(2,ylower+i-1) == 0)) then
                    rhs_values(k) = rhs_values(k) - uBdryX(2,ylower+i-1)*xibarx(i,j+1)/xedge
                else if((pcol == pncols).and.(j == localncols).and.(isDiriX(2,ylower+i-1) == 1)) then
                    rhs_values(k) = rhs_values(k) - right*p(i,j+1) - lambdax(i,j+1)*rhobarx(i,j+1)*gravX/xedge
                    values(r) = values(r) - right
                else
                    values(r+2) = right
                    values(r) = values(r) - right
                    rhs_values(k) = rhs_values(k) - lambdax(i,j+1)*rhobarx(i,j+1)*gravX/xedge
                end if

                if((prow == pnrows).and.(i == localnrows).and.(isDiriY(xlower+j-1,2) == 0)) then
                    rhs_values(k) = rhs_values(k) - uBdryY(xlower+j-1,2)*xibary(i+1,j)/yedge
                else if((prow == pnrows).and.(i == localnrows).and.(isDiriY(xlower+j-1,2) == 1)) then
                    rhs_values(k) = rhs_values(k) - up*p(i+1,j) - lambday(i+1,j)*rhobary(i+1,j)*gravY/yedge
                    values(r) = values(r) - up
                else
                    values(r+4) = up
                    values(r) = values(r) - up
                    rhs_values(k) = rhs_values(k) - lambday(i+1,j)*rhobary(i+1,j)*gravY/yedge
                end if

                r = r + 5
                k = k + 1

            end do
        end do

        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 5, stencil_indices, values, ierr)
        call HYPRE_StructMatrixAssemble(A, ierr)

        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, b, ierr)
        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, x, ierr)
        call HYPRE_StructVectorInitialize(b, ierr)
        call HYPRE_StructVectorInitialize(x, ierr)
        call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs_values, ierr)
        x_values(:) = initial_x_guess(:)
        call HYPRE_StructVectorSetBoxValues(x, ilower, iupper, x_values, ierr)
        call HYPRE_StructVectorAssemble(b, ierr)
        call HYPRE_StructVectorAssemble(x, ierr)

        call HYPRE_StructSMGCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
        call HYPRE_StructSMGSetMaxIter(solver, 50, ierr)
        call HYPRE_StructSMGSetTol(solver, 1.0e-07, ierr)
        call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
        call HYPRE_StructSMGSetNumPreRelax(solver, 1, ierr)
        call HYPRE_StructSMGSetNumPostRelax(solver, 1, ierr)
        !call HYPRE_StructSMGSetLogging(solver, 1, ierr)
        !call HYPRE_StructSMGSetPrintLevel(solver, 3, ierr)

        call HYPRE_StructSMGSetup(solver, A, b, x, ierr)
        call HYPRE_StructSMGSolve(solver, A, b, x, ierr)
        call HYPRE_StructSMGGetNumIterations(solver, num_iter, ierr)
        if((ierr /= 0).or.(num_iter == 0)) then
            print *, 'The solver error.', ierr, num_iter
            call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
        end if

        call HYPRE_StructVectorGetBoxValues(x, ilower, iupper, x_values, ierr)

        k = 1
        do i = 1, localnrows
            do j = 1, localncols
                p(i,j) = x_values(k)
                k = k + 1
            end do
        end do

        initial_x_guess(:) = x_values(:)

        deallocate(rhs_values)
        deallocate(x_values)
        deallocate(values)

        call HYPRE_StructGridDestroy(grid, ierr)
        call HYPRE_StructStencilDestroy(stencil, ierr)
        call HYPRE_StructMatrixDestroy(A, ierr)
        call HYPRE_StructVectorDestroy(b, ierr)
        call HYPRE_StructVectorDestroy(x, ierr)
        call HYPRE_StructSMGDestroy(solver, ierr)

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(pcol /= 1) then
            allocate(psent(localnrows))
            psent = p(1:localnrows,1)
            call MPI_IBSEND(psent, localnrows, MPI_DOUBLE_PRECISION, myid-pnrows, myid+num_procs, MPI_COMM_WORLD, requestl, ierr)
            deallocate(psent)
        end if

        if(pcol /= pncols) then
            allocate(psent(localnrows))
            psent = p(1:localnrows,localncols)
            call MPI_IBSEND(psent, localnrows, MPI_DOUBLE_PRECISION, myid+pnrows, myid+num_procs, MPI_COMM_WORLD, requestr, ierr)
            deallocate(psent)
        end if

        if(prow /= 1) then
            allocate(psent(localncols))
            psent = p(1,1:localncols)
            call MPI_IBSEND(psent, localncols, MPI_DOUBLE_PRECISION, myid-1, myid+num_procs, MPI_COMM_WORLD, requestd, ierr)
            deallocate(psent)
        end if

        if(prow /= pnrows) then
            allocate(psent(localncols))
            psent = p(localnrows,1:localncols)
            call MPI_IBSEND(psent, localncols, MPI_DOUBLE_PRECISION, myid+1, myid+num_procs, MPI_COMM_WORLD, requestu, ierr)
            deallocate(psent)
        end if

        if(pcol /= 1) then
            recvbuffersize = localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-pnrows, myid-pnrows+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            p(1:localnrows,0) = recvbuffer
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            recvbuffersize = localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+pnrows, myid+pnrows+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            p(1:localnrows,localncols+1) = recvbuffer
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            recvbuffersize = localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-1, myid-1+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            p(0,1:localncols) = recvbuffer
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            recvbuffersize = localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+1, myid+1+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            p(localnrows+1,1:localncols) = recvbuffer
            deallocate(recvbuffer)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

    end subroutine computePres

    subroutine computeMoleFrac()

        implicit none

        real(kind=8), dimension(:), pointer :: c
        real(kind=8) :: div, ctotal
        real(kind=8), dimension(:), pointer :: xsent, recvbuffer
        integer :: recvbuffersize
        integer :: requestl, requestr, requestu, requestd
        integer :: i, j, m, k

        integer :: ierr
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: position
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)

        allocate(c(Nc))

        do j = 1, localncols
            do i = 1, localnrows
                do m = 1, Nc
                    div = (xbarx(m,i,j+1)*ux(i,j+1)*xibarx(i,j+1) - xbarx(m,i,j)*ux(i,j)*xibarx(i,j)) &!
                        /(xs(j+1)-xs(j)) + (xbary(m,i+1,j)*uy(i+1,j)*xibary(i+1,j) - &!
                        xbary(m,i,j)*uy(i,j)*xibary(i,j))/(ys(i+1)-ys(i))
                    c(m) = (src(m,xlower+j-1,ylower+i-1)-div)*(timeEnd/nt)/poro(xlower+j-1,ylower+i-1) + x(m,i,j)*xi(i,j)
                    if(c(m)<0) then
                        print *, 'Please tune the time step.', c(m),m,i,j,myid
                        call MPI_Abort(MPI_COMM_WORLD,ierr)
                    end if
                end do
                ctotal = 0
                do m = 1, Nc
                    ctotal = ctotal + c(m)
                end do
                do m = 1, Nc
                    x(m,i,j) = c(m)/ctotal
                end do
                do m =1, Nc
                    if(x(m,i,j) < 1.D-99) then
                        x(m,i,j) = 0.0
                    end if
                end do
            end do
        end do

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)
        
        if(pcol /= 1) then
            allocate(xsent(Nc*localnrows))
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    xsent(k) = x(m,i,1)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(xsent, Nc*localnrows, MPI_DOUBLE_PRECISION, myid-pnrows, myid+num_procs, MPI_COMM_WORLD, requestl, ierr)
            deallocate(xsent)
        end if

        if(pcol /= pncols) then
            allocate(xsent(Nc*localnrows))
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    xsent(k) = x(m,i,localncols)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(xsent, Nc*localnrows, MPI_DOUBLE_PRECISION, myid+pnrows, myid+num_procs, MPI_COMM_WORLD, requestr, ierr)
            deallocate(xsent)
        end if

        if(prow /= 1) then
            allocate(xsent(Nc*localncols))
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    xsent(k) = x(m,1,i)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(xsent, Nc*localncols, MPI_DOUBLE_PRECISION, myid-1, myid+num_procs, MPI_COMM_WORLD, requestd, ierr)
            deallocate(xsent)
        end if

        if(prow /= pnrows) then
            allocate(xsent(Nc*localncols))
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    xsent(k) = x(m,localnrows,i)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(xsent, Nc*localncols, MPI_DOUBLE_PRECISION, myid+1, myid+num_procs, MPI_COMM_WORLD, requestu, ierr)
            deallocate(xsent)
        end if

        if(pcol /= 1) then
            recvbuffersize = Nc*localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-pnrows, myid-pnrows+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    x(m,i,0) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            recvbuffersize = Nc*localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+pnrows, myid+pnrows+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    x(m,i,localncols+1) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            recvbuffersize = Nc*localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-1, myid-1+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    x(m,0,i) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            recvbuffersize = Nc*localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+1, myid+1+num_procs, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    x(m,localnrows+1,i) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        deallocate(c)

    end subroutine computeMoleFrac

    subroutine computeVel()

        implicit none

        integer :: i, j, ierr

        ! compute the total new velocities in x direction
        if(pncols == 1) then
            do i = 1, localnrows
                if(isDiriX(1,ylower+i-1) == 1) then
                    ux(i,1) = -lambdax(i,1)/xibarx(i,1)*((p(i,1)-p(i,0))*2/(xs(2)-xs(1)) - rhobarx(i,1)*gravX)
                else
                    ux(i,1) = uBdryX(1,ylower+i-1)
                end if
            end do
            do i = 1, localnrows
                if(isDiriX(2,ylower+i-1) == 1) then
                    ux(i,localncols+1) = -lambdax(i,localncols+1)/xibarx(i,localncols+1)*((p(i,localncols+1) &!
                        -p(i,localncols))*2/(xs(nx+1)-xs(nx)) - rhobarx(i,localncols+1)*gravX)
                else
                    ux(i,localncols+1) = uBdryX(2,ylower+i-1)
                end if
            end do
            do j = 2, localncols
                do i = 1, localnrows
                    ux(i,j) = -lambdax(i,j)/xibarx(i,j)*((p(i,j)-p(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - rhobarx(i,j)*gravX)
                end do
            end do
        elseif(pcol==1) then
            do i = 1, localnrows
                if(isDiriX(1,ylower+i-1) == 1) then
                    ux(i,1) = -lambdax(i,1)/xibarx(i,1)*((p(i,1)-p(i,0))*2/(xs(2)-xs(1)) - rhobarx(i,1)*gravX)
                else
                    ux(i,1) = uBdryX(1,ylower+i-1)
                end if
            end do
            do i = 1, localnrows
                do j = 2, localncols+1
                    ux(i,j) = -lambdax(i,j)/xibarx(i,j)*((p(i,j)-p(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - rhobarx(i,j)*gravX)
                end do
            end do
        elseif(pcol==pncols) then
            do i = 1, localnrows
                if(isDiriX(2,ylower+i-1) == 1) then
                    ux(i,localncols+1) = -lambdax(i,localncols+1)/xibarx(i,localncols+1)*((p(i,localncols+1)-p(i,localncols))*2/ &!
                        (xs(nx+1)-xs(nx)) - rhobarx(i,localncols+1)*gravX)
                else
                    ux(i,localncols+1) = uBdryX(2,ylower+i-1)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    ux(i,j) = -lambdax(i,j)/xibarx(i,j)*((p(i,j)-p(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - rhobarx(i,j)*gravX)
                end do
            end do
        else
            do j = 1, localncols+1
                do i = 1, localnrows
                    ux(i,j) = -lambdax(i,j)/xibarx(i,j)*((p(i,j)-p(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - rhobarx(i,j)*gravX)
                end do
            end do
        end if

        ! compute the total new velocities in y direction
        if(pnrows == 1) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,1) == 1) then
                    uy(1,j) = -lambday(1,j)/xibary(1,j)*((p(1,j)-p(0,j))*2/(ys(2)-ys(1)) - rhobary(1,j)*gravY)
                else
                    uy(1,j) = uBdryY(xlower+j-1,1)
                end if
            end do
            do j = 1, localncols
                if(isDiriY(xlower+j-1,2) == 1) then
                    uy(localnrows+1,j) = -lambday(localnrows+1,j)/xibary(localnrows+1,j)*((p(localnrows+1,j) &!
                        -p(localnrows,j))*2/(ys(ny+1)-ys(ny)) - rhobary(localnrows+1,j)*gravY)
                else
                    uy(localnrows+1,j) = uBdryY(xlower+j-1,2)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows
                    uy(i,j) = -lambday(i,j)/xibary(i,j)*((p(i,j)-p(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - rhobary(i,j)*gravY)
                end do
            end do
        elseif(prow==1) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,1) == 1) then
                    uy(1,j) = -lambday(1,j)/xibary(1,j)*((p(1,j)-p(0,j))*2/(ys(2)-ys(1)) - rhobary(1,j)*gravY)
                else
                    uy(1,j) = uBdryY(xlower+j-1,1)
                end if
            end do
            do i = 2, localnrows+1
                do j = 1, localncols
                    uy(i,j) = -lambday(i,j)/xibary(i,j)*((p(i,j)-p(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - rhobary(i,j)*gravY)
                end do
            end do
        elseif(prow==pnrows) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,2) == 1) then
                    uy(localnrows+1,j) = -lambday(localnrows+1,j)/xibary(localnrows+1,j)*((p(localnrows+1,j)-p(localnrows,j))*2/ &!
                        (ys(ny+1)-ys(ny)) - rhobary(localnrows+1,j)*gravY)
                else
                    uy(localnrows+1,j) = uBdryY(xlower+j-1,2)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    uy(i,j) = -lambday(i,j)/xibary(i,j)*((p(i,j)-p(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - rhobary(i,j)*gravY)
                end do
            end do
        else
            do j = 1, localncols
                do i = 1, localnrows+1
                    uy(i,j) = -lambday(i,j)/xibary(i,j)*((p(i,j)-p(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - rhobary(i,j)*gravY)
                end do
            end do
        end if

    end subroutine computeVel

    subroutine finalize()

        implicit none

        integer :: ierr
        real(kind=8) :: timefinish

        deallocate(xs)
        deallocate(ys)
        deallocate(ts)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(poro)
        deallocate(src)
        deallocate(isDiriX)
        deallocate(isDiriY)
        deallocate(pBdryX)
        deallocate(pBdryY)
        deallocate(pInit)
        deallocate(xBdryX)
        deallocate(xBdryY)
        deallocate(xInit)
        deallocate(uBdryX)
        deallocate(uBdryY)
        deallocate(ct)
        deallocate(cp)
        deallocate(af)
        deallocate(mw)
        deallocate(cv)
        deallocate(delta)

        deallocate(p)
        deallocate(ux)
        deallocate(uy)
        deallocate(lambdax)
        deallocate(lambday)
        deallocate(Kxxbar)
        deallocate(Kyybar)
        deallocate(rho)
        deallocate(rhobarx)
        deallocate(rhobary)
        deallocate(x)
        deallocate(xbarx)
        deallocate(xbary)
        deallocate(xi)
        deallocate(xibarx)
        deallocate(xibary)
        deallocate(visc)
        deallocate(moleincell)
        deallocate(moleincell_old)
        deallocate(deri_xi_p)
        deallocate(deri_xi_n)
        deallocate(initial_x_guess)

        if(myid == 0) then
            close(40)
            close(50)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timefinish = MPI_Wtime()
        print *, 'Elapsed time = ', timefinish-timestart, ' seconds.'

        call MPI_Finalize(ierr)

    end subroutine finalize

end module RST_singlePhaseFlow
