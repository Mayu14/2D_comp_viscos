module NumericalEx
implicit none
    interface
        subroutine NumericalExperimentation(iStep,MConf,NE,CC,Geom,UG) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
            use StructVar_Mod
            integer, intent(in) :: iStep
            type(Configulation), intent(in) :: MConf
            type(NumericalExperiment), intent(inout) :: NE
            type(CellCenter), intent(inout) :: CC
            type(Geometry), intent(in) :: Geom
            type(UnstructuredGrid), intent(in), optional :: UG
        end subroutine NumericalExperimentation
    end interface

end module NumericalEx
