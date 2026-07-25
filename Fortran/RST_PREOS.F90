
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_PREOS

    use RST_model
    use RST_globalData

    implicit none

contains

    subroutine PREOS(local_x, local_p, local_xi, local_rho, local_deri_xi_p, local_deri_xi_n, local_moleincell)

        real(kind=8), dimension(:), pointer, intent(in) :: local_x
        real(kind=8), intent(in) :: local_p
        real(kind=8), intent(out) :: local_xi
        real(kind=8), intent(out) :: local_rho
        real(kind=8), intent(out) :: local_deri_xi_p
        real(kind=8), dimension(:), pointer, intent(out) :: local_deri_xi_n
        real(kind=8), dimension(:), pointer, intent(out) :: local_moleincell
        real(kind=8) :: ZG, ZL, Ztemp
        real(kind=8), dimension(:), pointer :: am, bm
        real(kind=8) :: a, b, bigA, bigB
        real(kind=8) :: mm, C2, C3, ctotal, xi_bs
        real(kind=8), dimension(:), pointer :: lambda, alpha, c, C1toC2, ZRP, ZRM
        integer :: i, j, m, number, size
        real(kind=8) :: molesum
        real(kind=8) :: deri_A_p, deri_B_p, deri_Z_p, deri_xi_p_bs
        real(kind=8), dimension(:), pointer :: deri_a_n, deri_b_n, deri_bigA_n, deri_bigB_n, deri_Z_n
        real(kind=8) :: PP(4), Z(3)

        allocate(lambda(Nc))
        do i = 1, Nc
            lambda(i) = 0.37464 + 1.5423*af(i) - 0.26992*af(i)**2
        end do

        allocate(alpha(Nc))
        do i = 1, Nc
            alpha(i) = (1+lambda(i)*(1-dsqrt(1.D0*Temp/ct(i))))**2
        end do

        allocate(am(Nc))
        allocate(bm(Nc))
        do i = 1, Nc
            am(i) = 0.45724*alpha(i)*R**2*ct(i)**2/cp(i)
            bm(i) = 0.077796*R*ct(i)/cp(i)
        end do

        a = 0
        do i = 1, Nc
            do j = 1, Nc
                a = a + local_x(i)*local_x(j)*(1-delta(i,j))*dsqrt(1.D0*am(i)*am(j))
            end do
        end do

        b = 0
        do i = 1, Nc
            b = b + local_x(i)*bm(i)
        end do

        deallocate(lambda)
        deallocate(alpha)

        bigA = a*local_p/(R*Temp)**2
        bigB = b*local_p/(R*Temp)

        PP(1) = 1.0
        PP(2) = -(1-bigB)
        PP(3) = bigA-3*bigB**2-2*bigB
        PP(4) = -(bigA*bigB-bigB**2-bigB**3)
        call getCubicRoot(PP, Z, size)

        if(phase == 'l') then

            allocate(ZRP(size))
            number = 0
            do i = 1, size
                if(Z(i)>0) then
                    ZRP(number+1) = Z(i)
                    number = number + 1
                end if
            end do
            allocate(ZRM(number))
            do i = 1, number
                ZRM(i) = ZRP(i)
            end do
            ZL = ZRM(1)
            do i = 2, number
                if(ZRM(i) < ZL) then
                    ZL = ZRM(i)
                end if
            end do
            deallocate(ZRM)
            deallocate(ZRP)

            xi_bs = local_p/(R*Temp*ZL)

            allocate(c(Nc))
            allocate(C1toC2(Nc))
            do i = 1, Nc
                C1toC2(i) = 110.07*af(i)**4 - 83.807*af(i)**3 + 18.926* &!
                    af(i)**2 - 1.6348*af(i) - 0.0066
            end do
            C2 = 2.013645*1.D-3
            C3 = 0.89
            do i = 1, Nc
                c(i) = C1toC2(i)*C2 + C2*(Temp/ct(i)-C3)**2
            end do
            ctotal = 0
            do i = 1, Nc
                ctotal = ctotal + local_x(i)*c(i)*mw(i)
            end do
            local_xi = 1/(1/xi_bs + ctotal)

            local_rho = 0
            do m = 1, Nc
                local_rho = local_rho + local_x(m)*mw(m)
            end do
            local_rho = local_rho*local_xi

            deri_A_p = a/(R*Temp)**2
            deri_B_p = b/(R*Temp)
            deri_Z_p = -(deri_B_p*ZL**2+(deri_A_p-2*(1+3*bigB)*deri_B_p)*ZL - (deri_A_p*bigB+(bigA-2*bigB-3*bigB**2)*deri_B_p))/ &!
                (3*ZL**2-2*(1-bigB)*ZL+(bigA-2*bigB-3*bigB**2))
            deri_xi_p_bs = 1/(R*Temp*ZL) - local_p/(R*Temp*ZL**2)*deri_Z_p
            local_deri_xi_p = (1/(1+ctotal*xi_bs)**2)*deri_xi_p_bs

            deallocate(c)
            deallocate(C1toC2)

            Ztemp = ZL

        else

            ZG = Z(1)
            do i = 2, size
                if(Z(i) > ZG) then
                    ZG = Z(i)
                end if
            end do
            local_xi = local_p/(R*Temp*ZG)
            local_rho = 0
            do m = 1, Nc
                local_rho = local_rho + local_x(m)*mw(m)
            end do
            local_rho = local_rho*local_xi
            deri_A_p = a/(R*Temp)**2
            deri_B_p = b/(R*Temp)
            deri_Z_p = -(deri_B_p*ZG**2+(deri_A_p-2*(1+3*bigB)*deri_B_p)*ZG - (deri_A_p*bigB+(bigA-2*bigB-3*bigB**2)*deri_B_p))/ &!
                (3*ZG**2-2*(1-bigB)*ZG+(bigA-2*bigB-3*bigB**2))
            local_deri_xi_p = 1/(R*Temp*ZG) - local_p/(R*Temp*ZG**2)*deri_Z_p

            Ztemp = ZG

        end if

        molesum = 0
        do m = 1, Nc
            local_moleincell(m) = local_xi*local_x(m)*(xs(2)-xs(1))*(ys(2)-ys(1))*poro(2,2)!!!!!!!!!!!!!!!
            molesum = molesum + local_moleincell(m)
        end do

        allocate(deri_a_n(Nc))
        allocate(deri_b_n(Nc))

        do i = 1, Nc
            deri_a_n(i) = 0
            do j = 1, Nc
                deri_a_n(i) = deri_a_n(i) + 2*local_x(j)/molesum**2*(1-delta(i,j))*dsqrt(1.D0*am(i)*am(j))
            end do
        end do

        do i = 1, Nc
            deri_b_n(i) = bm(i)/molesum
        end do

        allocate(deri_bigA_n(Nc))
        allocate(deri_bigB_n(Nc))
        allocate(deri_Z_n(Nc))

        do i = 1, Nc

            deri_bigA_n(i) = local_p*deri_a_n(i)/R**2/Temp**2

            deri_bigB_n(i) = local_p*deri_b_n(i)/R/Temp

            deri_Z_n(i) = (-Ztemp**2*deri_bigB_n(i)-Ztemp*(deri_bigA_n(i)-2*deri_bigB_n(i)-6*bigB*deri_bigB_n(i))+ &!
                bigA*deri_bigB_n(i)+bigB*deri_bigA_n(i)-2*bigB*deri_bigB_n(i)-3*bigB**2*deri_bigB_n(i))/ &!
                (3*Ztemp**2-2*Ztemp*(1-bigB)+bigA-2*bigB-3*bigB**2)

            local_deri_xi_n(i) = (-local_p*deri_Z_n(i))/R/Temp/Ztemp**2

        end do

        deallocate(am)
        deallocate(bm)
        deallocate(deri_a_n)
        deallocate(deri_b_n)
        deallocate(deri_bigA_n)
        deallocate(deri_bigB_n)
        deallocate(deri_Z_n)

    end subroutine PREOS

    subroutine getCubicRoot(local_p, local_X, size)

        real(kind=8), parameter :: TwoPi = 8.D0*atan(1.D0)
        real(kind=8), intent(in) :: local_p(4)
        real(kind=8), intent(in out) :: local_X(3)
        integer, intent(out) :: size
        real(kind=8) :: a, b, c, d, Alph, Beta, Delt, R1, R2, tht
        
        local_X = 0.D0
        a = local_p(1)
        b = local_p(2)/(3.D0*a)
        c = local_p(3)/(6.D0*a)
        d = local_p(4)/(2.D0*a)

        Alph = -b*b*b + 3.D0*b*c - d
        Beta =  b*b - 2.D0*c
        Delt = Alph*Alph-Beta*Beta*Beta

        if(Delt > 0.D0) then
            tht = Alph+dsqrt(1.D0*Delt); R1 = sign(abs(tht)**(1.D0/3.D0), 1.D0*tht)
            tht = Alph-dsqrt(1.D0*Delt); R2 = sign(abs(tht)**(1.D0/3.D0), 1.D0*tht)
            local_X(1) = -b+R1+R2
            size = 1
        else if(Delt == 0.D0) then
            R1 = sign(abs(Alph)**(1.D0/3.D0), 1.D0*Alph)
            if(R1 == 0.D0) then
                local_X(1) = -b
                size = 1
            else
                local_X(1) = -b+2.D0*R1
                local_X(2) = -b-R1
                size = 2
            end if
        else if(Delt < 0.D0) then
            tht = acos(Alph/(dsqrt(1.D0*Beta)*Beta))
            local_X(1)  = -b+2.D0*dsqrt(1.D0*Beta)*cos(tht/3.D0)
            local_X(2)  = -b+2.D0*dsqrt(1.D0*Beta)*cos((tht+TwoPi)/3.D0)
            local_X(3)  = -b+2.D0*dsqrt(1.D0*Beta)*cos((tht-TwoPi)/3.D0)
            size = 3
        end if

    end subroutine getCubicRoot

end module RST_PREOS
