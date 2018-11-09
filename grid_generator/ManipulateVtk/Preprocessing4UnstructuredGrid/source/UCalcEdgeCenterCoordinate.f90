subroutine UCalcEdgeCenterCoordinate(UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    iELM => iElement
    call GravityCenterOfTriangle
    call MidPointOfLine

return
contains
!_________________________________________________
    subroutine GravityCenterOfTriangle
    implicit none
        do iElement=1,UG%Tri%Total
            iVertex1 = UG%Tri%Point(iELM,1) !三角形の頂点の点番号
            iVertex2 = UG%Tri%Point(iELM,2)
            iVertex3 = UG%Tri%Point(iELM,3)
            !重心
            UG%CD%Cell(iELM,:) = (UG%CD%Point(iVertex1,:) + UG%CD%Point(iVertex2,:) &
                                        &  + UG%CD%Point(iVertex3,:))/3.0d0
        end do
    return
    end subroutine GravityCenterOfTriangle

!_________________________________________________

    subroutine MidPointOfLine
    implicit none
        do iElement=1, UG%Line%Total
            iEndPoint1 = UG%Line%Point(iELM,1)
            iEndPoint2 = UG%Line%Point(iELM,2)
            UG%CD%Edge(iELM,:) = 0.5d0*(UG%CD%Point(iEndPoint1,:)+UG%CD%Point(iEndPoint2,:))
        end do
    return
    end subroutine MidPointOfLine

end subroutine UCalcEdgeCenterCoordinate
