!***********************************/
!	Name:Roeスキームによって中心差分の流束を計算するプログラム
!	Alias:GetGeneralFlux
!	Description:計算結果はRoe流束の第1項に相当する
!	Type:RoeAverage
!	Input:Geometory,RoeAverage
!	Output:RoeAverage
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元に対応
!***********************************/
subroutine GetGeneralFlux(Geom,RA)
    use StructVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Geometry), intent(in) :: Geom
    type(RoeAverage), intent(inout) :: RA


        if(Geom%Dimension == 1) then !次元1のとき
            RA%GeneralFlux(1,1,:) = RA%SideQuantity(2,:)

            RA%GeneralFlux(2,1,:) = RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:) &
                &   + Gmin1*(RA%SideQuantity(3,:)-0.5d0*RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:))

            RA%GeneralFlux(3,1,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:) &
                &   * (Gamma*RA%SideQuantity(3,:)-0.5d0*Gmin1*RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:))

        else if(Geom%Dimension == 2) then !次元2のとき
            if(RA%Direction == 1) then !x方向流束
                RA%GeneralFlux(1,1,:) = RA%SideQuantity(2,:)

                RA%GeneralFlux(2,1,:) = RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:) &
                    &   + Gmin1*(RA%SideQuantity(4,:) &
                    &   - 0.5d0*(RA%SideQuantity(2,:)**2 + RA%SideQuantity(3,:)**2)/RA%SideQuantity(1,:))

                RA%GeneralFlux(3,1,:) = RA%SideQuantity(2,:)*RA%SideQuantity(3,:)/RA%SideQuantity(1,:)

                RA%GeneralFlux(4,1,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:) &
                    &   * (Gamma*RA%SideQuantity(4,:) &
                    &   - 0.5d0*Gmin1*(RA%SideQuantity(2,:)**2 + RA%SideQuantity(3,:)**2)/RA%SideQuantity(1,:))

            else if(RA%Direction == 2) then !y方向流束
                RA%GeneralFlux(1,1,:) = RA%SideQuantity(3,:)

                RA%GeneralFlux(2,1,:) = RA%SideQuantity(2,:)*RA%SideQuantity(3,:)/RA%SideQuantity(1,:)

                RA%GeneralFlux(3,1,:) = RA%SideQuantity(3,:)**2/RA%SideQuantity(1,:) &
                    &   + Gmin1*(RA%SideQuantity(4,:) &
                    &   - 0.5d0*(RA%SideQuantity(2,:)**2 + RA%SideQuantity(3,:)**2)/RA%SideQuantity(1,:))

                RA%GeneralFlux(4,1,:) = RA%SideQuantity(3,:)/RA%SideQuantity(1,:) &
                    &   * (Gamma*RA%SideQuantity(4,:) &
                    &   - 0.5d0*Gmin1*(RA%SideQuantity(2,:)**2 + RA%SideQuantity(3,:)**2)/RA%SideQuantity(1,:))
            end if
        !else if(Geom%Dimension == 3) then !次元3のとき
        end if

    return
end subroutine GetGeneralFlux
