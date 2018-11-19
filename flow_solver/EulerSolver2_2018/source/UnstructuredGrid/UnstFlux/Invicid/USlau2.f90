!***********************************/
!	Name:SLAU2で非粘性流束を計算するプログラム
!	Alias:USlau2
!	Description:隣接する2つのセルから基礎変数を受け取って流束を計算する
!	Type:
!	Input:
!	Output:
!	Note:基礎変数MUSCLとの併用を前提としているため，入力は基礎変数
!	Author:Akitaka Toyota
!	Date:2018.11.19
!	Update:
!	Other:
!***********************************/
subroutine USlau2(UConf, UG, UCC, UCE) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(in) :: UCC
    type(CellEdge), intent(inout) :: UCE

    double precision, allocatable :: PhaiLR(:, :), Normal(:)    ! (1:密度，x速度，y速度，z速度, エンタルピー, 2:L(minus)=1, R(plus)=2
    double precision :: PresL, PresR, VelNormL, VelNormR
    double precision :: SoundV_Ave !,SoundV_L, SoundV_R
    double precision :: NormalVelocity_L, NormalVelocity_R
    double precision :: AbsNormalVelocity_ave, AbsNormalVelocity_P, AbsNormalVelocity_M   ! PMはプラスマイナスの意味
    double precision :: Mach_L, Mach_R, Mach_hat
    double precision :: param_Chi, param_g   ! Chi =　カイ = xの筆記体は文献によってLambdaと表記されることがある
    double precision :: pres_P_a, pres_M_a, gamma_hr
    double precision :: mass_flux, pressure_flux
    iDim = UG%GM%Dimension

    allocate(PhaiLR(5, 2), Normal(5))
    Normal = 0.0d0

    do iEdge=1, UG%GI%Edges !すべての界面について
        ! set quantity of cell between surface iEdge
        !if(UConf%UseMUSCL == 1) then
        PresL = UCE%RebuildQunatity(5,1,1,2,iEdge)
        VelNormL = sum(UCE%RebuildQunatity(2:4,1,1,2,iEdge)**2)

        PhaiLR(1,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)
        PhaiLR(2,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(2,1,1,2,iEdge)
        PhaiLR(3,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(3,1,1,2,iEdge)
        PhaiLR(4,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(4,1,1,2,iEdge)
        PhaiLR(5,1) = InverseGmin1 * Gamma * PresL / UCE%RebuildQunatity(1,1,1,2,iEdge) + 0.5d0 * VelNormL

        PresR = UCE%RebuildQunatity(5,1,1,1,iEdge)
        VelNormR = sum(UCE%RebuildQunatity(2:4,1,1,1,iEdge)**2)

        PhaiLR(1,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)
        PhaiLR(2,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(2,1,1,1,iEdge)
        PhaiLR(3,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(3,1,1,1,iEdge)
        PhaiLR(4,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(4,1,1,1,iEdge)
        PhaiLR(5,2) = InverseGmin1 * Gamma * PresR / UCE%RebuildQunatity(1,1,1,1,iEdge) + 0.5d0 * VelNormR

        !else   ! UCC%PrimitiveVariablesからPhaiを計算する場合はこっちを利用(空間1次精度計算用)
        !end if
        Normal(2:4) = UG%GM%Normal(iEdge, :)

        ! get sound veocity average
        SoundV_Ave = 0.5d0 * (sqrt(Gamma * PresL / PhaiLR(1, 1)) + sqrt(Gamma * PresR / PhaiLR(1, 2)))

        ! get absolute normal velocity of each cell
        NormalVelocity_L = dot_product(PhaiLR(2:4, 1), UG%GM%Normal(iEdge, :))
        NormalVelocity_R = dot_product(PhaiLR(2:4, 2), UG%GM%Normal(iEdge, :))

        AbsNormalVelocity_ave = (PhaiLR(1, 1) * abs(NormalVelocity_L) + PhaiLR(1, 2) * abs(NormalVelocity_R)) &
                                & / (PhaiLR(1, 1) + PhaiLR(1, 2))

        ! calc local mach number
        Mach_L = NormalVelocity_L / SoundV_Ave
        Mach_R = NormalVelocity_R / SoundV_Ave
        Mach_hat = min(1.0d0, sqrt(0.5d0 * (VelNormL + VelNormR)) / SoundV_Ave)

        ! calc parameter param_Chi and param_g
        param_Chi = (1.0d0 - Mach_hat) ** 2
        param_g = -max(min(Mach_L, 0.0d0), -1.0d0) * min(max(Mach_R, 0.0d0), 1.0d0)


        ! calc absolute normal velocity Plus and Minus
        AbsNormalVelocity_P = (1.0d0 - param_g) * AbsNormalVelocity_ave + param_g * abs(NormalVelocity_L)
        AbsNormalVelocity_M = (1.0d0 - param_g) * AbsNormalVelocity_ave + param_g * abs(NormalVelocity_R)

        ! get mass_flux
        mass_flux = 0.5d0 * (PhaiLR(1, 1) * (NormalVelocity_L + AbsNormalVelocity_P) &
                         & + PhaiLR(1, 2) * (NormalVelocity_R - AbsNormalVelocity_M) &
                         & - param_Chi / SoundV_Ave * (PhaiLR(5, 2) - PhaiLR(5, 1)))    ! 最後の圧力の差分式は、プラスマイナス逆である可能性あり

        ! get pressure flux
        if(abs(Mach_L) > 1) then
            pres_P_a = 0.5d0 * (1.0d0 + sign(1.0d0, Mach_L))
        else
            pres_P_a = 0.25d0 * ((Mach_L + 1.0d0) ** 2) * (2.0d0 - Mach_L)
        end if

        if(abs(Mach_R) > 1) then
            pres_M_a = 0.5d0 * (1.0d0 - sign(1.0d0, Mach_R))
        else
            pres_M_a = 0.25d0 * ((Mach_R - 1.0d0) ** 2) * (2.0d0 + Mach_R)
        end if

        !gamma_hr = ! high resolution化は後回し

        pressure_flux = 0.5d0 * (PresL + PresR) + 0.5d0 * (pres_P_a - pres_M_a) * (PresL - PresR) &
                        & + sqrt(0.5d0 * (VelNormL + VelNormR) * (pres_P_a + pres_M_a) * (0.5d0*(PhaiLR(1,1) + PhaiLR(1,2))) * SoundV_Ave)


        UCE%NormalFluxDiff(:,1,1,1,iEdge) = &
        &   0.5d0 * (mass_flux + abs(mass_flux)) * PhaiLR(:, 1) &
        & + 0.5d0 * (mass_flux - abs(mass_flux)) * PhaiLR(:, 2) &
        & + pressure_flux * Normal

    end do

    return
end subroutine USlau2
