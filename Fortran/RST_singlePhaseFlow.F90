
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

    implicit none

contains

    subroutine initialize(modelCase)

        type(model), intent(in out) :: modelCase
       
        logical :: alive
        character(len=50) :: fmhtxt, fmrtxt
        integer :: i, j, m, ierr

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

        allocate(p(ny+2, nx+2))
        allocate(ux(ny, nx+1))
        ux = 0
        allocate(uy(ny+1, nx))
        uy = 0
        allocate(lambdax(ny, nx+1))
        lambdax(:,:) = 0
        allocate(lambday(ny+1, nx))
        lambday(:,:) = 0
        allocate(Kxxbar(ny, nx+1))
        allocate(Kyybar(ny+1, nx))
        allocate(rho(ny, nx))
        allocate(rhobarx(ny, nx+1))
        rhobarx(:,:) = 0
        allocate(rhobary(ny+1, nx))
        rhobary(:,:) = 0
        allocate(x(Nc, ny, nx))
        allocate(xbarx(Nc, ny, nx+1))
        xbarx(:,:,:) = 0
        allocate(xbary(Nc, ny+1, nx))
        xbary(:,:,:) = 0
        allocate(xi(ny, nx))
        allocate(xibarx(ny, nx+1))
        xibarx(:,:) = 0
        allocate(xibary(ny+1, nx))
        xibary(:,:) = 0
        allocate(visc(ny, nx))
        allocate(moleincell(Nc, ny, nx))
        allocate(moleincell_old(Nc, ny, nx))
        allocate(deri_xi_p(ny, nx))
        allocate(deri_xi_n(Nc,ny,nx))

        p(1, 2:nx+1) = pBdryY(1:nx, 1)
        p(ny+2, 2:nx+1) = pBdryY(1:nx, 2)
        p(2:ny+1, 1) = pBdryX(1, 1:ny)
        p(2:ny+1, nx+2) = pBdryX(2, 1:ny)
        p(2:ny+1,2:nx+1) = transpose(pInit(1:nx, 1:ny))

        uBdryX(1, 1:ny) = -uBdryX(1, 1:ny)
        uBdryY(1:nx, 1) = -uBdryY(1:nx, 1)
        ux(1:ny, 1) = uBdryX(1, 1:ny)
        ux(1:ny, nx+1) = uBdryX(2, 1:ny)
        uy(1, 1:nx) = uBdryY(1:nx, 1)
        uy(ny+1, 1:nx) = uBdryY(1:nx, 2)

        Kxxbar(1:ny,1) = Kxx(1,1:ny)
        Kxxbar(1:ny,nx+1) = Kxx(nx, 1:ny)
        do i = 1, ny
            do j = 2, nx
                Kxxbar(i,j) = (xs(j+1)-xs(j-1)) / ((xs(j)-xs(j-1))/Kxx(j-1,i)+(xs(j+1)-xs(j))/Kxx(j,i))
            end do
        end do

        Kyybar(1,1:nx) = Kyy(1:nx,1)
        Kyybar(ny+1,1:nx) = Kyy(1:nx, ny)
        do i = 2, ny
            do j = 1, nx
                Kyybar(i,j) = (ys(i+1)-ys(i-1)) / ((ys(i)-ys(i-1))/Kyy(j,i-1)+(ys(i+1)-ys(i))/Kyy(j,i))
            end do
        end do

        do m = 1, Nc
            do i = 1, ny
                do j = 1, nx
                    x(m,i,j) = xInit(m,j,i)
                end do
            end do
        end do

        totalmole = 0.0
        t = 2

        inquire(file = soludoc, exist = alive)
        if(.not.alive) then
            call system('mkdir '//trim(adjustl(soludoc)))
        end if

        fmhtxt = trim(adjustl(soludoc))//'/soln_1PhFlw_moleHistory.txt'
        fmrtxt = trim(adjustl(soludoc))//'/soln_1PhFlw_moleRatioHistory.txt'

        open(unit=40, file=trim(adjustl(fmhtxt)), status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if
        open(unit=50, file=trim(adjustl(fmrtxt)), status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
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

        real(kind=8), dimension(:), pointer :: xtemp, deri_xi_ntemp, moleincelltemp
        real(kind=8), dimension(:), pointer :: leftmole
        real(kind=8) :: totaldesiredleftmole
        real(kind=8) :: visctemp, ignore1
        real(kind=8), dimension(:), pointer  :: ignore2, ignore3
        logical :: isZero
        integer :: i, j, m

        allocate(xtemp(Nc))
        allocate(deri_xi_ntemp(Nc))
        allocate(moleincelltemp(Nc))
        do j = 1, nx
            do i = 1, ny
                if(t > 2) then
                    moleincell_old(1:Nc,i,j) = moleincell(1:Nc,i,j)
                end if
                xtemp(1:Nc) = x(1:Nc,i,j)
                call PREOS( xtemp, p(i+1,j+1), xi(i,j), rho(i,j), deri_xi_p(i,j), deri_xi_ntemp, moleincelltemp )
                deri_xi_n(1:Nc,i,j) = deri_xi_ntemp(1:Nc)
                moleincell(1:Nc,i,j) = moleincelltemp(1:Nc)
                if(t == 2) then
                    moleincell_old(1:Nc,i,j) = moleincell(1:Nc,i,j)
                end if
                visc(i,j) = viscosity( xtemp, xi(i,j), p(i+1,j+1) )
            end do
        end do

        allocate(leftmole(Nc))
        leftmole = 0
        do j = 1, nx
            do i = 1, ny
                do m = 1, Nc
                    leftmole(m) = leftmole(m) + moleincell(m,i,j)
                end do
            end do
        end do

        totaldesiredleftmole = 0.0
        do m = 2, Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m)
        end do

        if(t == 2) then
            totalmole = totaldesiredleftmole
        end if

        write(40, fmt="(es12.5)") (totalmole-totaldesiredleftmole)/totalmole
        write(50, fmt="(es12.5)") totaldesiredleftmole/leftmole(1)

        deallocate(leftmole)

        allocate(ignore2(Nc))
        allocate(ignore3(Nc))
        do i = 1, ny
            isZero = .true.
            do m = 1, Nc
                if(xBdryX(m,1,i) /= 0) then
                    isZero = .false.
                    exit
                end if
            end do
            if((ux(i,1)>0).and.(p(i+1,2)/=0).and.(.not.isZero)) then
                xtemp(1:Nc) = xBdryX(1:Nc,1,i)
                call PREOS( xtemp, p(i+1,2), xibarx(i,1), rhobarx(i,1), ignore1, ignore2, ignore3 )
                visctemp = viscosity( xtemp, xibarx(i,1), p(i+1,2) )
                lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visctemp
                xbarx(1:Nc,i,1) = xtemp(1:Nc)
            else
                xibarx(i,1) = xi(i,1)
                rhobarx(i,1) = rho(i,1)
                lambdax(i,1) = Kxxbar(i,1)*xibarx(i,1)/visc(i,1)
                xbarx(1:Nc,i,1) = x(1:Nc,i,1)
            end if
        end do

        do i = 1, ny
            isZero = .true.
            do m = 1, Nc
                if(xBdryX(m,2,i) /= 0) then
                    isZero = .false.
                    exit
                end if
            end do
            if((ux(i,nx+1)<0).and.(p(i+1,nx+1)/=0).and.(.not.isZero)) then
                xtemp(1:Nc) = xBdryX(1:Nc,2,i)
                call PREOS( xtemp, p(i+1,nx+1), xibarx(i,nx+1), rhobarx(i,nx+1), ignore1, ignore2, ignore3 )
                visctemp = viscosity( xtemp, xibarx(i,nx+1), p(i+1,nx+1) )
                lambdax(i,nx+1) = Kxxbar(i,nx+1)*xibarx(i,nx+1)/visctemp
                xbarx(1:Nc,i,nx+1) = xtemp(1:Nc)
            else
                xibarx(i,nx+1) = xi(i,nx)
                rhobarx(i,nx+1) = rho(i,nx)
                lambdax(i,nx+1) = Kxxbar(i,nx+1)*xibarx(i,nx+1)/visc(i,nx)
                xbarx(1:Nc,i,nx+1) = x(1:Nc,i,nx)
            end if
        end do

        do j = 2, nx
            do i = 1, ny
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

        do j = 1, nx
            isZero = .true.
            do m = 1, Nc
                if(xBdryY(m,j,1) /= 0) then
                    isZero = .false.
                    exit
                end if
            end do
            if((uy(1,j)>0).and.(p(2,j+1)/=0).and.(.not.isZero)) then
                xtemp(1:Nc) = xBdryY(1:Nc,j,1)
                call PREOS( xtemp, p(2,j+1), xibary(1,j), rhobary(1,j), ignore1, ignore2, ignore3 )
                visctemp = viscosity( xtemp, xibary(1,j), p(2,j+1) )
                lambday(1,j) = Kyybar(1,j)*xibary(1,j)/visctemp
                xbary(1:Nc,1,j) = xtemp(1:Nc)
            else
                xibary(1,j) = xi(1,j)
                rhobary(1,j) = rho(1,j)
                lambday(1,j) = Kyybar(1,j)*xibary(1,j)/visc(1,j)
                xbary(1:Nc,1,j) = x(1:Nc,1,j)
            end if
        end do

        do j = 1, nx
            isZero = .true.
            do m = 1, Nc
                if(xBdryY(m,j,2) /= 0) then
                    isZero = .false.
                    exit
                end if
            end do
            if((uy(ny+1,j)<0).and.(p(ny+1,j+1)/=0).and.(.not.isZero)) then
                xtemp(1:Nc) = xBdryY(1:Nc,j,2)
                call PREOS( xtemp, p(ny+1,j+1), xibary(ny+1,j), rhobary(ny+1,j), ignore1, ignore2, ignore3 )
                visctemp = viscosity( xtemp, xibary(ny+1,j), p(ny+1,j+1) )
                lambday(ny+1,j) = Kyybar(ny+1,j)*xibary(ny+1,j)/visctemp
                xbary(1:Nc,ny+1,j) = xtemp(1:Nc)
            else
                xibary(ny+1,j) = xi(ny,j)
                rhobary(ny+1,j) = rho(ny,j)
                lambday(ny+1,j) = Kyybar(ny+1,j)*xibary(ny+1,j)/visc(ny,j)
                xbary(1:Nc,ny+1,j) = x(1:Nc,ny,j)
            end if
        end do

        do j = 1, nx
            do i = 2, ny
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

        deallocate(ignore2)
        deallocate(ignore3)
        deallocate(xtemp)
        deallocate(deri_xi_ntemp)
        deallocate(moleincelltemp)

    end subroutine computeParameters

    subroutine computePres()

        real(kind=8), dimension(:,:), pointer :: A
        real(kind=8), dimension(:), pointer :: b
        real(kind=8) :: xedge, yedge, ledge, redge, uedge, dedge
        real(kind=8) :: up, down, left, right
        integer :: i, j, r, m
        integer :: INFO
        integer, dimension(:), pointer :: IPIV

        allocate(A(nx*ny,nx*ny))
        allocate(b(nx*ny))
        allocate(IPIV(nx*ny))
        A(:,:) = 0.0
        b(:) = 0.0
        IPIV(:) = 0.0

        do i = 2, ny+1
            do j = 2, nx+1
            
                r = (i-2)*nx + (j-1)

                xedge = xs(j) - xs(j-1)
                yedge = ys(i) - ys(i-1)

                if(j /= 2) then
                    ledge = xs(j-1) - xs(j-2)
                else
                    ledge = 0
                end if

                if(j /= nx+1) then
                    redge = xs(j+1) - xs(j)
                else
                    redge = 0
                end if

                if(i /= ny+1) then
                    uedge = ys(i+1) - ys(i)
                else
                    uedge = 0
                end if

                if(i /= 2) then
                    dedge = ys(i-1) - ys(i-2)
                else
                    dedge = 0
                end if

                A(r,r) = poro(j-1,i-1)*deri_xi_p(i-1,j-1)/(timeEnd/nt)
                b(r) = poro(j-1,i-1)*deri_xi_p(i-1,j-1)/(timeEnd/nt)*p(i,j)
                do m = 1, Nc
                    b(r) = b(r) + src(m,j-1,i-1) - poro(j-1,i-1)*deri_xi_n(m,j-1,i-1)* &!
                        (moleincell(m,j-1,i-1)-moleincell_old(m,j-1,i-1))/(timeEnd/nt)
                end do

                up = -2*lambday(i,j-1)/yedge/(yedge+uedge)
                down = -2*lambday(i-1,j-1)/yedge/(yedge+dedge)
                left = -2*lambdax(i-1,j-1)/xedge/(xedge+ledge)
                right = -2*lambdax(i-1,j)/xedge/(xedge+redge)

                if((i == ny+1).and.(isDiriY(j-1,2) == 0)) then
                    b(r) = b(r) - uBdryY(j-1,2)*xibary(i,j-1)/yedge
                elseif((i == ny+1).and.(isDiriY(j-1,2) == 1)) then
                    b(r) = b(r) - up*p(i+1,j) - lambday(i,j-1)*rhobary(i,j-1)*gravY/yedge
                    A(r,r) = A(r,r) - up
                else
                    A(r,r+nx) = up
                    A(r,r) = A(r,r) - up
                    b(r) = b(r) - lambday(i,j-1)*rhobary(i,j-1)*gravY/yedge
                end if

                if((i == 2).and.(isDiriY(j-1,1) == 0)) then
                    b(r) = b(r) + uBdryY(j-1,1)*xibary(i-1,j-1)/yedge
                elseif((i == 2).and.(isDiriY(j-1,1) == 1)) then
                    b(r) = b(r) - down*p(i-1,j) + lambday(i-1,j-1)*rhobary(i-1,j-1)*gravY/yedge
                    A(r,r) = A(r,r) - down
                else
                    A(r,r-nx) = down
                    A(r,r) = A(r,r) - down
                    b(r) = b(r) + lambday(i-1,j-1)*rhobary(i-1,j-1)*gravY/yedge
                end if

                if((j == 2).and.(isDiriX(1,i-1) == 0)) then
                    b(r) = b(r) + uBdryX(1,i-1)*xibarx(i-1,j-1)/xedge
                elseif((j == 2).and.(isDiriX(1,i-1) == 1)) then
                    b(r) = b(r) - left*p(i,j-1) + lambdax(i-1,j-1)*rhobarx(i-1,j-1)*gravX/xedge
                    A(r,r) = A(r,r) - left
                else
                    A(r,r-1) = left
                    A(r,r) = A(r,r) - left
                    b(r) = b(r) + lambdax(i-1,j-1)*rhobarx(i-1,j-1)*gravX/xedge
                end if

                if((j == nx+1).and.(isDiriX(2,i-1) == 0)) then
                    b(r) = b(r) - uBdryX(2,i-1)*xibarx(i-1,j)/xedge
                elseif((j == nx+1).and.(isDiriX(2,i-1) == 1)) then
                    b(r) = b(r) - right*p(i,j+1) - lambdax(i-1,j)*rhobarx(i-1,j)*gravX/xedge
                    A(r,r) = A(r,r) - right
                else
                    A(r,r+1) = right
                    A(r,r) = A(r,r) - right
                    b(r) = b(r) - lambdax(i-1,j)*rhobarx(i-1,j)*gravX/xedge
                end if

            end do
        end do

        call dgesv(nx*ny, 1, A, nx*ny, IPIV, b, nx*ny, INFO)

        do i = 2, ny+1
            do j = 2, nx+1
                p(i,j) = b((i-2)*nx+(j-1))
            end do
        end do

        deallocate(A)
        deallocate(b)
        deallocate(IPIV)

    end subroutine computePres

    subroutine computeMoleFrac()

        real(kind=8), dimension(:), pointer :: c
        real(kind=8) :: div, ctotal
        integer :: i, j, m
        
        allocate(c(Nc))

        do j = 1, nx
            do i = 1, ny
                do m = 1, Nc
                    div = (xbarx(m,i,j+1)*ux(i,j+1)*xibarx(i,j+1) - xbarx(m,i,j)*ux(i,j)*xibarx(i,j)) &!
                        /(xs(j+1)-xs(j)) + (xbary(m,i+1,j)*uy(i+1,j)*xibary(i+1,j) - &!
                        xbary(m,i,j)*uy(i,j)*xibary(i,j))/(ys(i+1)-ys(i))
                    c(m) = (src(m,j,i)-div)*(timeEnd/nt)/poro(j,i) + x(m,i,j)*xi(i,j)
                    if(c(m)<0) then
                        print *, 'Please tune the time step.', c(m)
                        stop
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

        deallocate(c)

    end subroutine computeMoleFrac

    subroutine computeVel()

        integer :: i, j

        do i = 1, ny
            if(isDiriX(1,i) == 1) then
                ux(i,1) = -lambdax(i,1)/xibarx(i,1) * ((p(i+1,2)-p(i+1,1))*2/(xs(2)-xs(1)) &!
                    - rhobarx(i,1)*gravX)
            else
                ux(i,1) = uBdryX(1,i)
            end if
        end do

        do i = 1, ny
            if(isDiriX(2,i) == 1) then
                ux(i,nx+1) = -lambdax(i,nx+1)/xibarx(i,nx+1) * ((p(i+1,nx+2)-p(i+1,nx+1))*2/(xs(nx+1)-xs(nx))&!
                    - rhobarx(i,nx+1)*gravX)
            else
                ux(i,nx+1) = uBdryX(2,i)
            end if
        end do

        do j = 2, nx
            do i = 1, ny
                ux(i,j) = -lambdax(i,j)/xibarx(i,j) * ((p(i+1,j+1)-p(i+1,j))*2/(xs(j+1)-xs(j-1)) - &!
                    rhobarx(i,j)*gravX)
            end do
        end do

        do j = 1, nx
            if(isDiriY(j,1) == 1) then
                uy(1,j) = -lambday(1,j)/xibary(1,j) * ((p(2,j+1)-p(1,j+1))*2/(ys(2)-ys(1)) - &!
                    rhobary(1,j)*gravY)
            else
                uy(1,j) = uBdryY(j,1)
            end if
        end do

        do j = 1, nx
            if(isDiriY(j,2) == 1) then
                uy(ny+1,j) = -lambday(ny+1,j)/xibary(ny+1,j) * ((p(ny+2,j+1)-p(ny+1,j+1))*2/(ys(ny+1)-ys(ny))&!
                    - rhobary(ny+1,j)*gravY)
            else
                uy(ny+1,j) = uBdryY(j,2)
            end if
        end do

        do j = 1, nx
            do i = 2, ny
                uy(i,j) = -lambday(i,j)/xibary(i,j) * ((p(i+1,j+1)-p(i,j+1))*2/(ys(i+1)-ys(i-1)) - &!
                    rhobary(i,j)*gravY)
            end do
        end do

    end subroutine computeVel

    subroutine finalize()

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

        close(40)
        close(50)

    end subroutine finalize

end module RST_singlePhaseFlow
