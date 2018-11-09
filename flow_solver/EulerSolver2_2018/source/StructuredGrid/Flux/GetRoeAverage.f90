!***********************************/
!	Name:Roe平均を計算するプログラム
!	Alias:GetRoeAverage
!	Description:隣接する2つのセルから保存変数を受け取ってRoe平均を計算する
!	Type:RoeAverage
!	Input:Geometory,RoeAverage
!	Output:RoeAverage
!	Note:内部関数としてRoeの荷重平均を求める関数が入っている
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元に対応
!***********************************/
subroutine GetRoeAverage(Geom,RA,SP4O)
    use StructVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none

    type(Geometry), intent(in) :: Geom
    type(RoeAverage), intent(inout) :: RA
    type(PrivateVar4OMP) :: SP4O


    SP4O%GRA%SqrtDensity(:) = dsqrt(RA%SideQuantity(1,:))
    SP4O%GRA%InverseDensitySum = 1.0d0 / (SP4O%GRA%SqrtDensity(1) + SP4O%GRA%SqrtDensity(2))

        !1:音速,2~:速度,最後:エンタルピ
        if(Geom%Dimension == 1) then
            SP4O%GRA%Velocity(1,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:)

            SP4O%GRA%Enthalpy(:) = Gamma*RA%SideQuantity(3,:)/ RA%SideQuantity(1,:) - Gmin1*0.5d0*SP4O%GRA%Velocity(1,:)**2

            RA%RoeAverage(2) = RoeLoadAverage(SP4O%GRA%Velocity(1,1),SP4O%GRA%Velocity(1,2))
            RA%RoeAverage(3) = RoeLoadAverage(SP4O%GRA%Enthalpy(1),SP4O%GRA%Enthalpy(2))
            RA%RoeAverage(1) = dsqrt(Gmin1*RA%RoeAverage(3)-0.5d0*RA%RoeAverage(2)**2)

        else if(Geom%Dimension == 2) then
            SP4O%GRA%Velocity(1,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:)
            SP4O%GRA%Velocity(2,:) = RA%SideQuantity(3,:)/RA%SideQuantity(1,:)
            SP4O%GRA%Enthalpy(:) = Gamma*RA%SideQuantity(4,:)/ RA%SideQuantity(1,:) &
                &   - Gmin1*0.5d0*(SP4O%GRA%Velocity(1,:)**2 + SP4O%GRA%Velocity(2,:)**2)

            RA%RoeAverage(2) = RoeLoadAverage(SP4O%GRA%Velocity(1,1),SP4O%GRA%Velocity(1,2))
            RA%RoeAverage(3) = RoeLoadAverage(SP4O%GRA%Velocity(2,1),SP4O%GRA%Velocity(2,2))
            RA%RoeAverage(4) = RoeLoadAverage(SP4O%GRA%Enthalpy(1),SP4O%GRA%Enthalpy(2))
            RA%RoeAverage(1) = dsqrt(Gmin1*RA%RoeAverage(4)-0.5d0*(RA%RoeAverage(2)**2 + RA%RoeAverage(3)**2))
        end if

    return

contains
    function RoeLoadAverage(QuantityL,QuantityR) result(Average)
    implicit none
    double precision,intent(in) :: QuantityL,QuantityR
    double precision :: Average

    Average =  (SP4O%GRA%SqrtDensity(1)*QuantityL + SP4O%GRA%SqrtDensity(2)*QuantityR)*SP4O%GRA%InverseDensitySum

    return
    end function RoeLoadAverage

end subroutine GetRoeAverage
