!***********************************/
!	Name:流束計算プログラムを制御するためのプログラム
!	Alias:CalcFlux
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%NormalFluxDiff
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元,1次精度と2次精度のみに対応
!***********************************/
subroutine CalcFlux(Conf,Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    interface
        subroutine UpwindFlux_Dim1(Geom,CE,CC) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
            use StructVar_Mod
            type(Geometry), intent(in) :: Geom
            type(CellEdge), intent(inout) :: CE
            type(CellCenter), intent(in), optional :: CC
        end subroutine UpwindFlux_Dim1
    end interface

    interface
        subroutine UpwindFlux_Dim2(Geom,CE,CC) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
            use StructVar_Mod
            type(Geometry), intent(in) :: Geom
            type(CellEdge), intent(inout) :: CE
            type(CellCenter), intent(in), optional :: CC
        end subroutine UpwindFlux_Dim2
    end interface

    if(Conf%UseMUSCL == 1) then
        call MUSCL(Conf,Geom,CC,CE)
    end if

    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)
    iVar => iVariable
    iFX => iFaceX
    iFY => iFaceY
    iFZ => iFaceZ
    iCX => iFaceX
    iCY => iFaceY
    iCZ => iFaceZ

    if(Conf%UseMUSCL == 0) then
        if(Geom%Dimension == 1) call UpwindFlux_Dim1(Geom,CE,CC)
        if(Geom%Dimension == 2) call UpwindFlux_Dim2(Geom,CE,CC)
        !if(Geom%Dimension == 3) call UpwindFlux_Dim3(Geom,CE,CC)
    else if(Conf%UseMUSCL == 1) then
        if(Geom%Dimension == 1) call UpwindFlux_Dim1(Geom,CE)
        if(Geom%Dimension == 2) call UpwindFlux_Dim2(Geom,CE)
    end if
    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)

return
end subroutine CalcFlux
