!***********************************/
!	Name:
!	Alias:ReRegistBoundaryCell
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:To delete relationship of interpolation on the Outer-Boundary-Cell
!	Author:Akitaka Toyota
!	Date:2018.02.06
!	Update:
!	Other:
!***********************************/
    subroutine ReRegistBoundaryCell(UG,MG,Geom) !(iGrid,GG,iGlobalCellNumber,Geom,CUV)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(UnstructuredGrid), intent(inout) :: UG
    type(MoveGrid), intent(inout) :: MG
    type(Geometry), intent(inout) :: Geom

    integer, allocatable :: LocalGridNum(:) !iCenterX, iCenterY

    integer ::  iCell!, iRealCell
    type(DonorElement), pointer :: DE
    type(DonorElement), pointer :: TempDE, NextDE
    integer :: iCount

    allocate(LocalGridNum(2))

    !do iCell = UG%GI%RealCells+1, UG%GI%AllCells
     !   if(UG%GM%CellType(iCell,1,1) == 4) then !Overset boundary only
    !        iRealCell = UG%VC%Cell(iCell,1)
        do iCell = 1, UG%GI%RealCells

            if(UG%GM%CellType(iCell,1,1) /= 2) then !not internal wall boudanry
                call GetLocalNumber2D(MG%NSR%C%NextCoordsInG(iCell,1:2),2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),LocalGridNum(1:2))

                if(Geom%Interpolated(LocalGridNum(1),LocalGridNum(2),1) == 0) then
                    !Delete
                    DE => MG%ADE(iCell)%DE
                    do while(associated(DE))
                        TempDE => DE%Next
                            nullify(DE%Next) !?
                        deallocate(DE)
                        DE => TempDE
                    end do
                    nullify(MG%ADE(iCell)%DE) !?

                    !Addition
                    allocate(DE)
                    DE%Cell = Local2Global(LocalGridNum(1),LocalGridNum(2),Geom%CellNumber(1))
                    DE%Grid = 0
                    nullify(DE%Next)
                    DE%SharedArea = 1.0d0
                    MG%ADE(iCell)%DE => DE

                    UG%GM%Interpolated(iCell,1,1) = 1
                end if

            else
                DE => MG%ADE(iCell)%DE
                iCount = 0
                NextDE => DE%Next
                if(associated(NextDE)) then
                    if(DE%Grid == 0) then
                        deallocate(DE)
                        MG%ADE(iCell)%DE => NextDE
                    end if
                end if

                do while(associated(NextDE))
                    if(NextDE%Grid == 0) then
                        DE%Next => NextDE%Next
                        deallocate(NextDE)
                    end if
                    DE => DE%Next
                    NextDE => NextDE%Next
                end do

            end if
        end do


    return
    end subroutine ReRegistBoundaryCell
