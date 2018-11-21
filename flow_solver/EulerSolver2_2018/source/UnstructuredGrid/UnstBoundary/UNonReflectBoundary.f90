!***********************************/
!	Name:無反射境界条件
!	Alias:UNonReflectBoundary
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.11.22
!	Update:
!	Other:
!***********************************/
subroutine UNonReflectBoundary(UG,UCC,iWorkCell)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iWorkCell

    integer :: iAdjacentLoaclEdge
    double precision :: length

    double precision, allocatable :: Velocity(:), AdjVelocity(:)    ! 仮想セル内速度, 隣接実セル内速度
    double precision, allocatable :: dUdn(:), dUdx(:), dUdy(:)   ! 境界を挟んだ基礎変数勾配の;n:外向き法線方向成分, x:x方向, y:y方向
    double precision :: alpha1=0.0d0, alpha2=0.0d0, kappa, r0, soundV   ! 変数名は元の論文(Taniguchi et.al, 2005)に準ずる (soundV = c)
    double precision, allocatable :: r(:), drdt(:, :) ! (r1, r2, r3), ((x or y), (1,2,3))
    double precision, allocatable :: A(:,:)
    double precision :: u1, u2, d1, d2, d3, m1, m2, m3, m4, l1, l2, e_tilde
    double precision :: rho

        allocate(Velocity(3), AdjVelocity(3), dUdn(5), dUdx(5), dUdy(5), r(3),drdt(2,3), A(4,4))

        call GradientBetweenBoundary

        call set_alpha1_alpha2

        soundV = Gamma * Gmin1 * (UCC%ConservedQuantity(5,iWorkCell, 1, 1)/UCC%ConservedQuantity(1,iWorkCell,1,1) - 0.5d0 * sum(Velocity**2))

        kappa = sqrt(alpha1 ** 2 + alpha2 ** 2)
        r0 = alpha1 * Velocity(1) + alpha2 * Velocity(2)
        r(1) = r0 - kappa * soundV
        r(2) = r0 + kappa * soundV
        r(3) = r0

        !(dxdt)123, (dydt)123
        drdt(1, 1) = Velocity(1) - sign(1.0d0, alpha1) / sqrt(alpha1**2 + alpha2**2) * soundV
        drdt(2, 1) = Velocity(2) - sign(1.0d0, alpha2) / sqrt(alpha1**(-2) + alpha2**(-2)) * soundV

        drdt(1, 2) = Velocity(1) + sign(1.0d0, alpha1) / sqrt(alpha1**2 + alpha2**2) * soundV
        drdt(2, 2) = Velocity(2) + sign(1.0d0, alpha2) / sqrt(alpha1**(-2) + alpha2**(-2)) * soundV

        drdt(:, 3) = Velocity(:)

        do iLoop = 1, 3
            if(dot_product(drdt(1:2, iLoop), UG%GM%Normal(iEdge, 1:2)) <= 0.0d0) then
                r(iLoop) = 0.0d0
            end if
        end do

        !make matrix A_tilde
        u1 = Velocity(1)
        u2 = Velocity(2)

        d1 = 0.5d0 * (r(2) + r(1))
        d2 = 0.5d0 * (r(2) - r(1))
        d3 = r(3)

        m1 = rho * (d1 - d3) / (kappa**2)
        m2 = (d1 - d3) / (soundV**2)
        m3 = d2 / (soundV * kappa)
        m4 = rho * d3

        e_tilde = 0.5d0 * (sum(Velocity**2) + InverseGmin1 * soundV**2)

        l1 = rho * d2 / (soundV * kappa)
        l2 = r0 * m1 + rho * d2 * e_tilde / (soundV * kappa)

        call set_matrix_A_tilde

        UCC%ConservedQuantity(:, iWorkCell, 1, 1) = UCC%ConservedQuantity(:, iWorkCell, 1, 1) - FixedTimeStep &
                & * matmul(A, sign(1.0d0, alpha1) * dUdx + sign(1.0d0, alpha2) * dUdy) / (abs(alpha1) + abs(alpha2))


    return

contains

    subroutine GradientBetweenBoundary
        implicit none

        iAdjacentCell = UG%VC%Cell(iCell,1) !隣接する実セル
        iAdjacentEdge = UG%VC%Edge(iCell) !共有界面
        iAdjacentLoaclEdge = UG%Line%Cell(iEdge, 1, 2)

        length = 2.0d0 * sqrt(dot_product(UG%GM%Width(iAdjacentCell, iAdjacentLoaclEdge, :), UG%GM%Width(iAdjacentCell, iAdjacentLoaclEdge, :)))

        rho = UCC%ConservedQuantity(1, iWorkCell, 1, 1)
        Velocity = UCC%ConservedQuantity(2:4, iWorkCell, 1, 1) / UCC%ConservedQuantity(1, iWorkCell, 1, 1)
        AdjVelocity = UCC%ConservedQuantity(2:4, iAdjacentCell, 1, 1) / UCC%ConservedQuantity(1, iAdjacentCell, 1, 1)

        dUdn(1) = (UCC%ConservedQuantity(1, iWorkCell, 1, 1) - UCC%ConservedQuantity(1, iAdjacentCell, 1, 1)) / length
        dUdn(2:4) = (Velocity - AdjVelocity) / length

        dUdn(5) = Gmin1 * &
            &   ((UCC%ConservedQuantity(5,iWorkCell,1,1) - 0.5d0 * UCC%ConservedQuantity(1,iWorkCell,1,1) &
            &    * dot_product(Velocity, Velocity)) &
            &  - (UCC%ConservedQuantity(5,iAdjacentCell,1,1) - 0.5d0 * UCC%ConservedQuantity(1,iAdjacentCell,1,1) &
            &    * dot_product(AdjVelocity,AdjVelocity))) / length

        dUdx = dUdn * UG%GM%Normal(iEdge, 1)
        dUdy = -dUdn * UG%GM%Normal(iEdge, 2)

        return
    end subroutine GradientBetweenBoundary


    subroutine set_alpha1_alpha2
        implicit none

        do iLoop = 1, 5
            if(dUdx(iLoop) == 0.0d0) then
                alpha2 = alpha2 + 1000000.0d0
            else
                alpha2 = alpha2 + dUdy(iLoop) / dUdx(iLoop)
            end if

            if(dUdy(iLoop) == 0.0d0) then
                alpha1 = alpha1 + 1000000.0d0
            else
                alpha1 = alpha1 + dUdx(iLoop) / dUdy(iLoop)
            end if
        end do

        if(abs(alpha2) > abs(alpha1)) then
            alpha1 = 0.25d0 * alpha1
            alpha2 = 1.0d0
        else
            alpha1 = 1.0d0
            alpha2 = 0.25d0 * alpha2
        end if

        return
    end subroutine set_alpha1_alpha2

    subroutine set_matrix_A_tilde
        implicit none

        A(1,1) = d3
        A(1,2) = alpha1 * l1
        A(1,3) = alpha2 * l1
        A(1,4) = m2
        A(2,1) = u1 * d3
        A(2,2) = alpha1 * (u1 * l1 + alpha1 * m1) + m4
        A(2,3) = alpha2 * (u1 * l1 + alpha1 * m1)
        A(2,4) = u1 * m2 + alpha1 * m3
        A(3,1) = u2 * d3
        A(3,2) = alpha1 * (u2 * l1 + alpha2 * m1)
        A(3,3) = alpha2 * (u2 * l1 + alpha2 * m1) + m4
        A(3,4) = u2 * m2 + alpha2 * m3
        A(4,1) = d3 * (e_tilde - soundV**2 * InverseGmin1)
        A(4,2) = alpha1 * l2 + u1 * rho * d3
        A(4,3) = alpha2 * l2 + u2 * rho * d3
        A(4,4) = e_tilde * m2 + r0 * m3 + d3 * InverseGmin1

        return
    end subroutine set_matrix_A_tilde

end subroutine UNonReflectBoundary
