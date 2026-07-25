
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_genPlotfileMatlab

    use RST_model
    use RST_globalData

    implicit none

contains

    subroutine genPlotfileMatlab()

        character(len=30) :: fpmatlabm
        character(len=2) :: charNc
        character(len=8) :: charLx, charLy
        character(len=10) :: chartimeEnd, charnt
        character(len=5) :: charnx, charny
        logical :: alive

        fpmatlabm = trim(adjustl(soludoc))//"/matlabplot.m"

        open(unit=10, file=trim(adjustl(fpmatlabm)), status='replace')

        write(10, fmt="(a)") "path('/Users/wuy/Dropbox/Research/MulticomponentFlow_matlab/2D', path);"
        write(charNc,'(i2)') Nc
        write(10, fmt="(a)") "model.Nc = "//charNc//";"
        write(charLx,'(f8.2)') Lx
        write(10, fmt="(a)") "Lx = "//charLx//";"
        write(charLy,'(f8.2)') Ly
        write(10, fmt="(a)") "Ly = "//charLy//";"
        write(chartimeEnd,'(f10.1)') timeEnd
        write(10, fmt="(a)") "timeEnd = "//chartimeEnd//";"
        write(10, fmt="(a)") "m4sat = 2;"
        write(charnx,'(i5)') nx
        write(10, fmt="(a)") "nx = "//charnx//";"
        write(charny,'(i5)') ny
        write(10, fmt="(a)") "ny = "//charny//";"
        write(charnt,'(i10)') nt
        write(10, fmt="(a)") "nt = "//charnt//";"
        write(10, fmt="(a)") "model.nx = nx;"
        write(10, fmt="(a)") "model.ny = ny;"
        write(10, fmt="(a)") "model.nt = nt;"
        write(10, fmt="(a)") "model.xs = (0:nx)*Lx/nx;"
        write(10, fmt="(a)") "model.ys = (0:ny)*Ly/ny;"
        write(10, fmt="(a)") "model.ts = (0:nt)*timeEnd/nt;"
        write(10, fmt="(a)") "fptxt = 'soln_1PhFlw_P_raw.txt';"
        write(10, fmt="(a)") "fuxtxt = 'soln_1PhFlw_Ux_raw.txt';"
        write(10, fmt="(a)") "fuytxt = 'soln_1PhFlw_Uy_raw.txt';"
        write(10, fmt="(a)") "fmftxt = [];"
        write(10, fmt="(a)") "for m = 1 : model.Nc"
        write(10, fmt="(a)") "    fk = ['soln_1PhFlw_X', num2str(m), '_raw.txt'];"
        write(10, fmt="(a)") "    fmftxt = [fmftxt; fk];"
        write(10, fmt="(a)") "end"
        write(10, fmt="(a)") "fmhtxt = 'soln_1PhFlw_moleHistory.txt';"
        write(10, fmt="(a)") "fmrtxt = 'soln_1PhFlw_moleRatioHistory.txt';"
        write(10, fmt="(a)") "model.soludoc = 'matlabplots';"
        inquire(file = trim(adjustl(soludoc))//'/matlabplots', exist = alive)
        if(.not.alive) then
            call system("mkdir "//trim(adjustl(soludoc))//"/matlabplots")
        end if
        write(10, fmt="(a)") "RST_plot( model, fptxt, fuxtxt, fuytxt, fmftxt, fmhtxt, fmrtxt);"
        write(10, fmt="(a)") "rmpath('/Users/wuy/Dropbox/Research/MulticomponentFlow_matlab/2D');"

        close(10)

    end subroutine genPlotfileMatlab

end module RST_genPlotfileMatlab