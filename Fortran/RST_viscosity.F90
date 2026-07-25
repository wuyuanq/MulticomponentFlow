
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_viscosity

    use RST_model
    use RST_globalData

    implicit none

contains

    function viscosity( local_x, local_xi, local_p ) result(mu)

        real(kind=8), dimension(:), pointer, intent(in) :: local_x
        real(kind=8), intent(in) :: local_xi
        real(kind=8), intent(in) :: local_p
        real(kind=8) :: mu

        real(kind=8), dimension(:), pointer :: local_cp, local_mw
        integer :: i, j, m
        real(kind=8), dimension(:), pointer :: CI, XMIU
        real(kind=8) :: RO, TR, S, S1, VS, SV, Z7PLUS, MW7PLUS, GM7PLUS, V7PLUS, ROC, ROR, ST, SM, SPI, F, COR
        logical :: isHeavy
        real(kind=8) :: A(16), B(4)
        real(kind=8) :: PPC, TPC, GML, VS1, PPR, TPR, CE, VISR

        allocate(local_cp(Nc))
        allocate(local_mw(Nc))

        local_cp(1:Nc) = 9.8692*1.D-6*cp(1:Nc)
        local_mw(1:Nc) = mw(1:Nc)*1000.0

        if(phase == 'l') then

            allocate(CI(Nc))
            do i = 1, Nc
                CI(i) = 1/(ct(i)**(1.D0/6.D0)/sqrt(local_mw(i))/local_cp(i)**0.666666666666D0)
            end do

            RO = local_xi*1.d-3

            allocate(XMIU(Nc))
            do i = 1, Nc
                TR = Temp/ct(i)
                if (TR <= 1.5D0) then
                    XMIU(i)=34.D-5*CI(i)*TR**.94D0
                else
                    XMIU(i)=17.78D-5*CI(i)*(4.58D0*TR-1.67D0)**0.625D0
                end if
            end do

            S = 0.D0
            S1 = 0.D0
            GML = 0.D0

            do i = 1, Nc
                S = S + local_x(i)*XMIU(i)*sqrt(local_mw(i))
                S1 = S1 + local_x(i)*sqrt(local_mw(i))
                GML = GML + local_x(i)*local_mw(i)
            end do
            RO = RO*GML
            VS = S/S1
            SV = 0.D0

            isHeavy = .false.
            if(.not.isHeavy) then
                do i = 1, Nc
                    SV = SV + local_x(i)*cv(i)*local_mw(i)
                end do
            else
                do i = 1, 6
                    SV = SV + local_x(i)*cv(i)
                end do

                Z7PLUS = 0.D0
                do i = 7, Nc
                    Z7PLUS = Z7PLUS + local_x(i)
                end do
                MW7PLUS = 430.7
                GM7PLUS = 0.988 ! AS GIVEN IN THE HEAVY OIL CASE
                V7PLUS = 21.573 + 0.015122*MW7PLUS - 27.656*GM7PLUS + 0.070615*MW7PLUS*GM7PLUS
                V7PLUS = V7PLUS*0.3048**3/0.453592 ! cub-ft/lb-mol -> cub-m/kg-mol

                SV = SV + Z7PLUS*V7PLUS
            end if

            ROC = 1.D0/SV
            if(.not.isHeavy) THEN
                ROR = RO/ROC/GML
            else
                ROR = RO/ROC
            end if

            ST = 0.D0
            SM = 0.D0
            SPI = 0.D0

            do i = 1, Nc
                SM = SM + local_x(i)*local_mw(i)
                ST = ST + local_x(i)*ct(i)
                SPI = SPI + local_x(i)*local_cp(i)
            end do

            CE = ST**(1.D0/6.D0)/sqrt(SM)/SPI**(2.D0/3.D0)

            ! AUGUST 16, 2007: If ROR is greater than 10, F is suspiciously large ...
            F = .1023D0+.023364D0*ROR+.058533D0*ROR**2- 0.040758D0*ROR**3  +.0093324D0*ROR**4

            COR = (F**4-1.D-4)/CE
            mu = VS + COR

            deallocate(CI)
            deallocate(XMIU)

        else if(phase == 'g') then

            A(1) = -2.4621182
            A(2) = 2.97054714
            A(3) = -0.286264054
            A(4) = 8.05420522*1.D-3
            A(5) = 2.80860949
            A(6) = -3.49803305
            A(7) = 0.36037302
            A(8) = -0.0104432413
            A(9) = -0.793385684
            A(10) = 1.39643306
            A(11) = -0.149144925
            A(12) = 4.41015512*1.D-3
            A(13) = 0.0839387178
            A(14) = -0.186408848
            A(15) = 0.0203367881
            A(16) = -6.09579263*1.D-4

            B(1) = 4.0
            B(2) = 4.0
            B(3) = 4.0
            B(4) = 4.0

            PPC = 0.0
            TPC = 0.0
            GML = 0.0

            do m = 1, Nc
                PPC = PPC + local_x(m)*local_cp(m)
                TPC = TPC + local_x(m)*ct(m)
                GML = GML + local_x(m)*local_mw(m)
            end do

            VS1 = (7.43+0.0133*GML)*(1.8*Temp)**1.5 &!
                /(1.8*Temp+75.4+13.9*GML)*1.0*1.D-4
            PPR = local_p*1.D-5/PPC
            TPR = Temp/TPC
            do i = 1, 4
                j = 4*i
                B(i) = A(j-3)+(A(j-2)+(A(j-1)+A(j)*PPR)*PPR)*PPR
            end do

            CE = B(1)+(B(2)+(B(3)+B(4)*TPR)*TPR)*TPR
            VISR = dexp(CE*1.D0)/TPR
            mu = VISR*VS1

        end if

        mu = mu * 0.001

        deallocate(local_cp)
        deallocate(local_mw)

    end function viscosity

end module RST_viscosity
