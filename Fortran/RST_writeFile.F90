
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

    implicit none

contains

    subroutine writeFile()

        integer :: i, j, m, ierr
        character(len=50) :: fptxt, fuxtxt, fuytxt
        character(len=50), dimension(:), pointer :: fmftxt
        character :: charm

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

        open(unit=10, file=fptxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if
        do j = 2, nx+1
            do i = 2, ny+1
                write(10, fmt='(es15.8)', iostat=ierr) p(i,j)
                if(ierr /= 0) then
                    print *, 'write file error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fuxtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if
        do j = 1, nx+1
            do i = 1, ny
                write(10, fmt='(es12.5)', iostat=ierr) ux(i,j)
                if(ierr /= 0) then
                    print *, 'write file error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fuytxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if
        do j = 1, nx
            do i = 1, ny+1
                write(10, fmt='(es12.5)', iostat=ierr) uy(i,j)
                if(ierr /= 0) then
                    print *, 'write file error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        do m = 1, Nc
            open(unit=10, file=fmftxt(m), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                stop
            end if
            do j = 1, nx
                do i = 1, ny
                    write(10, fmt='(es12.5)', iostat=ierr) x(m,i,j)
                    if(ierr /= 0) then
                        print *, 'write file error. ', ierr
                        stop
                    end if
                end do
            end do
            close(10)
        end do

        deallocate(fmftxt)

    end subroutine writeFile

end module RST_writeFile
