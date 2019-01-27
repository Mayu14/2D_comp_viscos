!***********************************/
!	Name:CountUpVertex型変数への代入処理
!	Alias:OInputCUV
!	Description:CountUpVertex型変数への代入処理
!	Type:
!	Input:
!	Output:
!	Note:第1引数でtarget側かdonor側かを決め，第2引数でN+1ステップかNステップかを決める(Donorのみ有効)
!	Author:Akitaka Toyota
!	Date:2017.12.23
!	Update:
!	Other:
!***********************************/
    subroutine InputCUV(iTargetOrDonor,iNp1OrN,iCellNum,UG,NSR,CUV,D_MG)
    use StructVar_Mod
    implicit none
    integer, intent(in) :: iTargetOrDonor, iNp1OrN !target要素かDonor要素か，N段階かN+1段階か
    integer, intent(in) :: iCellNum
    type(UnstructuredGrid), intent(in) :: UG
    type(NextStepRelation), intent(in)  :: NSR
    type(CountUpVertex), intent(inout) :: CUV

    type(MoveGrid), intent(in) :: D_MG

    if(iTargetOrDonor == 1) then !ターゲット格子であるとき(常にN+1段階のみ)
        CUV%iTargetNumber(0) = iCellNum
        CUV%iTargetNumber(1:3) = UG%Tri%Point(iCellNum,1:3) !
        CUV%TargetCoords(0,1:3) = NSR%C%NextCoordsInG(iCellNum,1:3) !更新後のグローバル系での値
        CUV%TargetCoords(1,1:3) = NSR%P%NextCoordsInG(CUV%iTargetNumber(1),1:3)
        CUV%TargetCoords(2,1:3) = NSR%P%NextCoordsInG(CUV%iTargetNumber(2),1:3)
        CUV%TargetCoords(3,1:3) = NSR%P%NextCoordsInG(CUV%iTargetNumber(3),1:3)
        CUV%TargetVolume = UG%GM%Volume(iCellNum)

        CUV%iTotalDonor = 0 !initialize

    else !ドナー格子であるとき 2
        CUV%iDonorNumber(0) = iCellNum !保有要素の要素番号を確保し
        CUV%iDonorNumber(1:3) = UG%Tri%Point(iCellNum,1:3)
        CUV%DonorVolume = UG%GM%Volume(iCellNum)

        !if(iNp1OrN == 0) then !N+1時間段階であるとき 0
        !    CUV%DonorCoords(0,1:3) = NSR%C%NextCoordsInG(CUV%iDonorNumber(0),1:3)
        !    CUV%DonorCoords(1,1:3) = NSR%P%NextCoordsInG(CUV%iDonorNumber(1),1:3)
        !    CUV%DonorCoords(2,1:3) = NSR%P%NextCoordsInG(CUV%iDonorNumber(2),1:3)
        !    CUV%DonorCoords(3,1:3) = NSR%P%NextCoordsInG(CUV%iDonorNumber(3),1:3)

        !else !N段階であるとき 1
        !    CUV%DonorCoords(0,1:3) = NSR%C%PresentCoordsInG(CUV%iDonorNumber(0),1:3)
        !    CUV%DonorCoords(1,1:3) = NSR%P%PresentCoordsInG(CUV%iDonorNumber(1),1:3)
        !    CUV%DonorCoords(2,1:3) = NSR%P%PresentCoordsInG(CUV%iDonorNumber(2),1:3)
        !    CUV%DonorCoords(3,1:3) = NSR%P%PresentCoordsInG(CUV%iDonorNumber(3),1:3)
        !end if
            CUV%DonorCoords(0,1:3) = D_MG%RC%Cell(CUV%iDonorNumber(0),1:3)
            CUV%DonorCoords(1,1:3) = D_MG%RC%Point(CUV%iDonorNumber(1),1:3)
            CUV%DonorCoords(2,1:3) = D_MG%RC%Point(CUV%iDonorNumber(2),1:3)
            CUV%DonorCoords(3,1:3) = D_MG%RC%Point(CUV%iDonorNumber(3),1:3)

    end if

    return
    end subroutine InputCUV
