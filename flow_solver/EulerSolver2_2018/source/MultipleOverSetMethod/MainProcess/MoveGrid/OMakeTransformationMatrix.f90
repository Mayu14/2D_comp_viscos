!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OOverSetS2U
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.24
!	Update:2017.11.24
!	Other:
!***********************************/
subroutine OMakeTransformationMatrix(MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Geometry), intent(in) :: Geom
    type(MoveGrid), intent(inout) :: MG
    double precision, allocatable :: TmpTransformation(:,:)
    !double precision, allocatable :: RotationalMatrix(:,:,:)

!    write(6,*) "Make Transformation Matrix of Moving Grid"
    allocate(TmpTransformation(4,4))
    if(Geom%Dimension == 2) then

        TmpTransformation = 0.0d0

        TmpTransformation(1,1) = cos(MG%GCD%Rotation(3))
        TmpTransformation(1,2) = sin(MG%GCD%Rotation(3))
        TmpTransformation(2,1) = -sin(MG%GCD%Rotation(3))
        TmpTransformation(2,2) = cos(MG%GCD%Rotation(3))

        TmpTransformation(3,3) = 1.0d0

        TmpTransformation(1:3,4) = MG%GCD%Translation(:)
        TmpTransformation(4,4) = 1.0d0

        MG%TM%Next2Current = TmpTransformation

        MG%TM%Current2Top = MG%TM%Next2Top

        MG%TM%Next2Top = matmul(MG%TM%Current2Top,MG%TM%Next2Current)

        call MakeInverseMatrix
    !else
    !    allocate(RotationalMatrix(3,3,3))
    !
    end if

return
contains

    subroutine MakeInverseMatrix
    implicit none

    MG%TM%Current2Next(1:3,1:3) = transpose(MG%TM%Next2Current(1:3,1:3)) !回転行列を反転させる
    MG%TM%Current2Next(1:3,4) = - matmul(transpose(MG%TM%Next2Current(1:3,1:3)),MG%TM%Next2Current(1:3,4)) !平行移動ベクトルの座標を変換する
    MG%TM%Current2Next(4,1:4) = MG%TM%Next2Current(4,1:4) !最後の行は同じ

    MG%TM%Top2Current = MG%TM%Top2Next

    MG%TM%Top2Next = matmul(MG%TM%Current2Next,MG%TM%Top2Current)

    return
    end subroutine MakeInverseMatrix

end subroutine OMakeTransformationMatrix
