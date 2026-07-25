
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_writeFile

    use RST_model
    use RST_globalData
    include 'mpif.h'

contains

    subroutine writeFile()

        implicit none

        character(len=50) :: fptxt, fuxtxt, fuytxt
        character(len=50), dimension(:), pointer :: fmftxt
        character :: charm
        real(kind=8), dimension(:,:), pointer :: allp
        real(kind=8), dimension(:,:), pointer :: allux
        real(kind=8), dimension(:,:), pointer :: alluy
        real(kind=8), dimension(:,:,:), pointer :: allx
        integer :: p_prow, p_pcol
        integer :: p_xlower, p_xupper, p_ylower, p_yupper
        integer :: local_size, buffersize
        real(kind=8), dimension(:), pointer :: buffer
        real(kind=8), dimension(:), pointer :: psent, uxsent, uysent, xsent
        integer :: position
        integer :: i, j, m, n, prociteration

        integer :: ierr
        integer :: request
        integer :: status(MPI_STATUS_SIZE)        
        integer :: commbuffer_size = MAX_COMMBUF
        real(kind=8) :: commbuffer(MAX_COMMBUF)

        local_size = nx*ny/num_procs
        buffersize = (Nc+3)*local_size

        allocate(buffer(buffersize))
        call MPI_BUFFER_ATTACH(commbuffer,commbuffer_size,ierr)

        if(myid /= 0) then

            allocate(psent(local_size))
            allocate(uxsent(local_size))
            allocate(uysent(local_size))
            allocate(xsent(local_size*Nc))

            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    psent(n) = p(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    uxsent(n) = ux(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    uysent(n) = uy(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    do m = 1, Nc
                        xsent(n) = x(m,i,j)
                        n = n + 1
                    end do
                end do
            end do

            position = 0
            call MPI_PACK(psent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(uxsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(uysent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(xsent, local_size*Nc, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)

            call MPI_IBSEND(buffer,buffersize,MPI_DOUBLE_PRECISION,0,myid+num_procs,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)

            deallocate(psent)
            deallocate(uxsent)
            deallocate(uysent)
            deallocate(xsent)

        else

            allocate(allP(ny,nx))
            allocate(allux(ny,nx+1))
            allocate(alluy(ny+1,nx))
            allocate(allx(Nc,ny,nx))

            allux(1:ny,nx+1) = uBdryX(2, 1:ny)
            alluy(ny+1,1:nx) = uBdryY(1:nx, 2)

            do j = 1, localncols
                do i = 1, localnrows
                    allP(i,j) = p(i,j)
                    allux(i,j) = ux(i,j)
                    alluy(i,j) = uy(i,j)
                end do
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    do m = 1, Nc
                        allx(m,i,j) = x(m,i,j)
                    end do
                end do
            end do

            do prociteration = 1, num_procs-1

                p_prow = mod(prociteration,pnrows)+1
                p_pcol = prociteration/pnrows+1

                p_xlower = (p_pcol-1)*localncols+1
                p_xupper = p_pcol*localncols
                p_ylower = (p_prow-1)*localnrows+1
                p_yupper = p_prow*localnrows

                call MPI_RECV(buffer, buffersize, MPI_DOUBLE_PRECISION, prociteration, &!
                    prociteration+num_procs, MPI_COMM_WORLD, status, ierr)

                n = 1
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allp(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allux(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        alluy(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        do m = 1, Nc
                            allx(m,i,j) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do
            end do

            fptxt = trim(adjustl(soludoc))//'/soln_1PhFlw_P_raw.txt'
            fuxtxt = trim(adjustl(soludoc))//'/soln_1PhFlw_Ux_raw.txt'
            fuytxt = trim(adjustl(soludoc))//'/soln_1PhFlw_Uy_raw.txt'
            allocate(fmftxt(Nc))
            do m = 1, Nc
                write(charm,'(i1)', iostat=ierr) m
                if(ierr /= 0) then
                    print *, 'write file error. ', ierr
                    stop
                end if
                fmftxt(m) = trim(adjustl(soludoc))//'/soln_1PhFlw_X'//charm//'_raw.txt'
            end do

            open(unit=10, file=trim(adjustl(fptxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                stop
            end if
            do j = 1, nx
                do i = 1, ny
                    write(10, fmt='(es15.8)', iostat=ierr) allp(i,j)
                    if(ierr /= 0) then
                        print *, 'write file error. ', ierr
                        stop
                    end if
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuxtxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                stop
            end if
            do j = 1, nx+1
                do i = 1, ny
                    write(10, fmt='(es12.5)', iostat=ierr) allux(i,j)
                    if(ierr /= 0) then
                        print *, 'write file error. ', ierr
                        stop
                    end if
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuytxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                stop
            end if
            do j = 1, nx
                do i = 1, ny+1
                    write(10, fmt='(es12.5)', iostat=ierr) alluy(i,j)
                    if(ierr /= 0) then
                        print *, 'write file error. ', ierr
                        stop
                    end if
                end do
            end do
            close(10)

            do m = 1, Nc
                open(unit=10, file=trim(adjustl(fmftxt(m))), status='replace', iostat=ierr)
                if(ierr /= 0) then
                    print *, 'open file error. ', ierr
                    stop
                end if
                do j = 1, nx
                    do i = 1, ny
                        write(10, fmt='(es12.5)', iostat=ierr) allx(m,i,j)
                        if(ierr /= 0) then
                            print *, 'write file error. ', ierr
                            stop
                        end if
                    end do
                end do
                close(10)
            end do

            deallocate(fmftxt)

            deallocate(allp)
            deallocate(allux)
            deallocate(alluy)
            deallocate(allx)

        end if

        deallocate(buffer)
        call MPI_BUFFER_DETACH(commbuffer,commbuffer_size,ierr)

    end subroutine writeFile

end module RST_writeFile
