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
!	Update:2017.12.27
!	Other:1次元と2次元に対応 !12/26 UGetRoeAverageの引数からUGを削除
!***********************************/
subroutine UGetRoeAverage(RA,P4O)
    use StructVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none

    type(RoeAverage), intent(inout) :: RA
    type(PrivateVar4OMP), intent(inout) :: P4O

        if(RA%SideQuantity(1,1) < 0.0d0)  then
            write(6,*) P4O%Cell1,P4O%Cell2
            write(6,*) RA%SideQuantity(1,1),RA%SideQuantity(1,2)
        end if

    P4O%GRA%SqrtDensity(:) = dsqrt(RA%SideQuantity(1,:))
    P4O%GRA%InverseDensitySum = 1.0d0 / (P4O%GRA%SqrtDensity(1) + P4O%GRA%SqrtDensity(2))

    !保存量からのRoe平均生成
        !1:音速,2~:速度,最後:エンタルピ

            P4O%GRA%Velocity(1,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:)
            P4O%GRA%Velocity(2,:) = RA%SideQuantity(3,:)/RA%SideQuantity(1,:)
            P4O%GRA%Velocity(3,:) = RA%SideQuantity(4,:)/RA%SideQuantity(1,:)
            P4O%GRA%Enthalpy(:) = Gamma*RA%SideQuantity(5,:)/ RA%SideQuantity(1,:) &
                &   - Gmin1*0.5d0*(P4O%GRA%Velocity(1,:)**2 + P4O%GRA%Velocity(2,:)**2 + P4O%GRA%Velocity(3,:)**2)

            RA%RoeAverage(2) = RoeLoadAverage(P4O%GRA%Velocity(1,1),P4O%GRA%Velocity(1,2))
            RA%RoeAverage(3) = RoeLoadAverage(P4O%GRA%Velocity(2,1),P4O%GRA%Velocity(2,2))
            RA%RoeAverage(4) = RoeLoadAverage(P4O%GRA%Velocity(3,1),P4O%GRA%Velocity(3,2))
            RA%RoeAverage(5) = RoeLoadAverage(P4O%GRA%Enthalpy(1),P4O%GRA%Enthalpy(2))

            if(Gmin1*RA%RoeAverage(5)-0.5d0*(sum(RA%RoeAverage(2:4)**2)) < 0) then
                write(6,*) "UGetRoeAverage Line42"
                write(6,*) "RightSide"
                write(6,*) "Vel",P4O%GRA%Velocity(1,1), P4O%GRA%Velocity(2,1),P4O%GRA%Velocity(3,1)
                write(6,*) "Ent",P4O%GRA%Enthalpy(1)
                write(6,*) "SqD",P4O%GRA%SqrtDensity(1)
                write(6,*) "LeftSide"
                write(6,*) "Vel",P4O%GRA%Velocity(1,2), P4O%GRA%Velocity(2,2),P4O%GRA%Velocity(3,2)
                write(6,*) "Ent",P4O%GRA%Enthalpy(2)
                write(6,*) "SqD",P4O%GRA%SqrtDensity(2)
                write(6,*) "RoeAverage"
                write(6,*) "DRL",RA%SideQuantity(1,1),RA%SideQuantity(1,2)
                write(6,*) "a^2",Gmin1*RA%RoeAverage(5)-0.5d0*(sum(RA%RoeAverage(2:4)**2))
                write(6,*) "Rh",RA%RoeAverage(5)
                write(6,*) "RVx",RA%RoeAverage(2)
                write(6,*) "RVy",RA%RoeAverage(3)
                write(6,*) "RVz",RA%RoeAverage(4)
                RetryFlag = 1
             end if


                RA%RoeAverage(1) = dsqrt(Gmin1*RA%RoeAverage(5)-0.5d0*(sum(RA%RoeAverage(2:4)**2)))

    return
contains
    function RoeLoadAverage(QuantityL,QuantityR) result(Average)
    implicit none
    double precision,intent(in) :: QuantityL,QuantityR
    double precision :: Average

    Average =  (P4O%GRA%SqrtDensity(1)*QuantityL + P4O%GRA%SqrtDensity(2)*QuantityR)*P4O%GRA%InverseDensitySum

    return
    end function RoeLoadAverage

end subroutine UGetRoeAverage
