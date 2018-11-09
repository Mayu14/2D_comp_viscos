subroutine CalcVector(Geom)
    use StructVar_Mod
    implicit none
    type(Geometry), intent(inout) :: Geom

    if(Geom%Dimension == 1) then
        Geom%Vector(1,1) = 1.0d0
    end if

    if(Geom%Dimension == 2) then
        Geom%Vector(1,1) = 1.0d0 !x方向面の法線ベクトル
        Geom%Vector(1,2) = 0.0d0
        Geom%Vector(2,1) = 0.0d0 !y方向面の法線ベクトル
        Geom%Vector(2,2) = 1.0d0
    end if

    if(Geom%Dimension == 3) then !Def:(接線ベクトル1).cross.(接線ベクトル2)=(法線ベクトル) (※crossは外積を示す)
        Geom%Vector(1,1) = 1.0d0 !x方向面の法線ベクトル,y方向面の接線ベクトル2,z方向面の接線ベクトル1
        Geom%Vector(2,1) = 0.0d0
        Geom%Vector(3,1) = 0.0d0
        Geom%Vector(1,2) = 0.0d0 !y方向面の法線ベクトル,x方向面の接線ベクトル1,z方向面の接線ベクトル2
        Geom%Vector(2,2) = 1.0d0
        Geom%Vector(3,2) = 0.0d0
        Geom%Vector(1,3) = 0.0d0 !z方向面の法線ベクトル,y方向面の接線ベクトル1,x方向面の接線ベクトル2
        Geom%Vector(2,3) = 0.0d0
        Geom%Vector(3,3) = 1.0d0
    end if

return
end subroutine CalcVector
