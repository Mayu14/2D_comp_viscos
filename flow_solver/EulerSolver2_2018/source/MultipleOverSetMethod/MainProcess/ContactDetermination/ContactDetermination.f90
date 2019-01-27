!***********************************/
!	Name:ローカル格子同士の接触判定プログラム
!	Alias:OContactDetermination
!	Description:ローカル格子の接触を2段階で判定する
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.20
!	Update:2017.12.20
!	Other:
!***********************************/
subroutine OContactDetermination(UG1,UG2,MG1,MG2,iTimeCombination)
    use StructVar_Mod
!    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG1,UG2 !Target,Donor
    type(MoveGrid), intent(in) :: MG1,MG2 !Target,Donor
    integer, intent(in) :: iTimeCombination !N+1段階とN+1段階のとき0，N段階が入っていると1
    integer :: iToF, BoundarySphereMethod

!    write(6,*) "Determination ""CONTACT"" between TargetGrid and DonorGrid"
!	1.Lp格子と，n+1段階のLx格子(x<p)またはn段階のLy格子(p=<y)格子について接触判定を行う
!   L1格子は2~N，L2格子は3~Nのように接触判定を行う
    iToF = BoundarySphereMethod(UG1%RoughRadius,UG2%RoughRadius,MG1%RC,MG2%RC,iTimeCombination)

    !if(iToF == 1) then
        !N+1段階とN+1段階の接触判定ならiKind=0
        !N+1段階とN段階の接触判定ならiKind=1
        !call GJKAlgorithm(iTimeCombination,MG1,MG2,UG1,UG2,iToF)
        !どうも間違いがあるようなのであとで修正する(自格子同士のときは問題なく接触判定できてるので，あからさまに接触していないケースが怪しいのだが...)

    !end if

    if(iToF == 1) then !格子が重なっているとき
        if(iTimeCombination == 0) then !格子番号的にUG2 < UG1
            UG1%Overlap(iTimeCombination+1,UG2%Number) = 1 !n+1段階の上位格子と被っているとき1
            UG2%Overlap(iTimeCombination+1,UG1%Number) = 2 !一括判定を行う内部境界を持ち，下側から参照する側として

        else !格子番号的にUG1 <= UG2 N+1とNの接触判定
            UG1%Overlap(iTimeCombination+1,UG2%Number) = 3 !n段階の下位格子と被っているとき2
        end if
    else
        UG1%Overlap(iTimeCombination+1,UG2%Number) = 0 !格子が重なっていない(or考慮しなくてよい)とき0
        if(iTimeCombination == 0) then !格子番号的にUG2 < UG1
            UG2%Overlap(iTimeCombination+1,UG1%Number) = 0 !一括判定を行う内部境界を持ち，下側から参照する側として
        end if
    end if


return
end subroutine OContactDetermination
