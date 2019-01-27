!***********************************/
!	Name:動的割当が正しく動作しているかどうかのテストプログラム
!	Alias:AllocTest
!	Description:とりあえず前配列の下限と上限を表示してみる
!	Type:Geom,CC,CE
!	Input:Geometory,CellEdge,CellCenter
!	Output:
!	Note:新しい変数を作ることがあればここに追記すること
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine AllocTest(CC,CE,Geom)
use StructVar_Mod
implicit none
    type(CellCenter),intent(in) :: CC
    type(CellEdge),intent(in) :: CE
    type(Geometry),intent(in) :: Geom

    write(6,*) "lbound(CC%ConservedQuantity)"
    write(6,*) lbound(CC%ConservedQuantity)
    write(6,*) ubound(CC%ConservedQuantity)
    write(6,*) "lbound(CC%GradientOfVariable)"
    write(6,*) lbound(CC%GradientOfVariable)
    write(6,*) ubound(CC%GradientOfVariable)
    write(6,*) "lbound(CC%PreviousQuantity)"
    write(6,*) lbound(CC%PreviousQuantity)
    write(6,*) ubound(CC%PreviousQuantity)
    write(6,*) "lbound(CC%PrimitiveVariable)"
    write(6,*) lbound(CC%PrimitiveVariable)
    write(6,*) ubound(CC%PrimitiveVariable)
    write(6,*) "lbound(CC%RungeKuttaG1)"
    write(6,*) lbound(CC%RungeKuttaG1)
    write(6,*) ubound(CC%RungeKuttaG1)
    write(6,*) "lbound(CC%RungeKuttaG3)"
    write(6,*) lbound(CC%RungeKuttaG3)
    write(6,*) ubound(CC%RungeKuttaG3)
    write(6,*) "lbound(CC%TimeWidth)"
    write(6,*) lbound(CC%TimeWidth)
    write(6,*) ubound(CC%TimeWidth)
    write(6,*) "lbound(CC%TmpTimeWidth)"
    write(6,*) lbound(CC%TmpTimeWidth)
    write(6,*) ubound(CC%TmpTimeWidth)
    write(6,*) "lbound(CC%VariationOfVariable)"
    write(6,*) lbound(CC%VariationOfVariable)
    write(6,*) ubound(CC%VariationOfVariable)
    write(6,*) "lbound(CE%LimiterFunction)"
    write(6,*) lbound(CC%LimiterFunction)
    write(6,*) ubound(CC%LimiterFunction)


    write(6,*) "lbound(CE%NormalFluxDiff)"
    write(6,*) lbound(CE%NormalFluxDiff)
    write(6,*) ubound(CE%NormalFluxDiff)
    write(6,*) "lbound(CE%NormalGradient)"
    write(6,*) lbound(CE%NormalGradient)
    write(6,*) ubound(CE%NormalGradient)
    write(6,*) "lbound(CE%RebuildQunatity)"
    write(6,*) lbound(CE%RebuildQunatity)
    write(6,*) ubound(CE%RebuildQunatity)
    write(6,*) "lbound(CE%TmpLimiterFunction)"
    write(6,*) lbound(CE%TmpLimiterFunction)
    write(6,*) ubound(CE%TmpLimiterFunction)

    write(6,*) "lbound(Geom%Area)"
    write(6,*) lbound(Geom%Area)
    write(6,*) ubound(Geom%Area)
    write(6,*) "lbound(Geom%Bound)"
    write(6,*) lbound(Geom%Bound)
    write(6,*) ubound(Geom%Bound)
    write(6,*) "lbound(Geom%CellNumber)"
    write(6,*) lbound(Geom%CellNumber)
    write(6,*) ubound(Geom%CellNumber)
    !write(6,*) "lbound(Geom%Vector)"
    !write(6,*) lbound(Geom%Vector)
    !write(6,*) ubound(Geom%Vector)
    write(6,*) "lbound(Geom%Volume)"
    write(6,*) lbound(Geom%Volume)
    write(6,*) ubound(Geom%Volume)
    write(6,*) "lbound(Geom%Width)"
    write(6,*) lbound(Geom%Width)
    write(6,*) ubound(Geom%Width)
    return
end subroutine AllocTest
