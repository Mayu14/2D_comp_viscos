!***********************************/
!	Name:影響域分離法のための配列確保
!	Alias:AllocCC4MB
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.31
!	Update:
!	Other:
!***********************************/
    subroutine AllocCC4MB(MConf,CC,CC4MB)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: MConf
    type(CellCenter), intent(in) :: CC !reference
    type(CellCenter), intent(inout) :: CC4MB

    allocate(CC4MB%ConservedQuantity(5,lbound(CC%ConservedQuantity,2):ubound(CC%ConservedQuantity,2),1,1)) !for calculation
    allocate(CC4MB%PrimitiveVariable(5,lbound(CC%PrimitiveVariable,2):ubound(CC%PrimitiveVariable,2),1,1))

    allocate(CC4MB%InterpolatedQuantity(5,lbound(CC%InterpolatedQuantity,2):ubound(CC%InterpolatedQuantity,2),1,1))

    !CC4MB%PrimitiveVariable !for output/muscl
    if(MConf%UseMUSCL == 1) then
        allocate(CC4MB%GradientOfVariable(5,lbound(CC%GradientOfVariable,2):ubound(CC%GradientOfVariable,2),1,1,3))
        allocate(CC4MB%LimiterFunction(5,lbound(CC%LimiterFunction,2):ubound(CC%LimiterFunction,2),1,1))
        allocate(CC4MB%VariationOfVariable(2,5,lbound(CC%VariationOfVariable,3):ubound(CC%VariationOfVariable,3),1,1))
    end if

    !if(MConf%UseOverSet == 1) then !for overset??? たぶんいらない
    !    allocate(CC4MB%InterpolatedQuantity(5,lbound(CC%InterpolatedQuantity,2):ubound(CC%InterpolatedQuantity,2),1,1)
    !end if

    if(MConf%UseRRK2 == 1) then
        allocate(CC4MB%PreviousQuantity(5,lbound(CC%PreviousQuantity,2):ubound(CC%PreviousQuantity,2),1,1))
        allocate(CC4MB%RungeKuttaG1(5,lbound(CC%RungeKuttaG1,2):ubound(CC%RungeKuttaG1,2),1,1))
        allocate(CC4MB%RungeKuttaG3(5,lbound(CC%RungeKuttaG3,2):ubound(CC%RungeKuttaG3,2),1,1))
    end if

    allocate(CC4MB%TimeWidth(lbound(CC%TimeWidth,1):ubound(CC%TimeWidth,1),1,1))
    allocate(CC4MB%TmpTimeWidth(lbound(CC%TmpTimeWidth,1):ubound(CC%TmpTimeWidth,1),1,1,1)) !for variable time step

    return
    end subroutine AllocCC4MB
