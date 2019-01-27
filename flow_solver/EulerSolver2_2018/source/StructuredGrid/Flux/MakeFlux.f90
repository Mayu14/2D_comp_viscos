subroutine MakeFlux(Geom,CC,CE)!没
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
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


    call Conserve2Primitive(Geom,CC)

    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)
    iVar => iVariable
    iFX => iFaceX
    iFY => iFaceY
    iFZ => iFaceZ
    iCX => iFaceX
    iCY => iFaceY
    iCZ => iFaceZ

    if(Geom%Dimension == 1) call UpwindFlux_Dim1(Geom,CE,CC)
    !if(Geom%Dimension == 2) call FluxDim2
    !if(Geom%Dimension == 3) call FluxDim3
    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)
    return

end subroutine MakeFlux
