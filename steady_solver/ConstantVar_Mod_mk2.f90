!***********************************/
!	Name:定数の宣言省略用モジュール
!	Alias:ConstantVar_Mod
!	Description:物理定数やその逆数，Runge-Kutta法のパラメータなどを収録
!	Type:double precision, integer
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2017.11.08
!	Other:
!***********************************/
module ConstantVar_Mod_mk2
implicit none

!Physical Constant
    double precision, parameter :: dPi = dacos(-1.0d0)
    double precision, parameter :: SpecificOfHeatRatio = 1.4d0

!Calculate Default
    double precision :: FixedTimeStep
    double precision :: DefaultTimeStep
    integer :: OutputInterval !結果の出力間隔
    integer :: IterationNumber !計算の反復回数
    integer :: CoreNumberOfCPU = 1
    integer :: GridNumber = 1!重合格子法にて使用する格子の数

!Constant value for Calculate
    double precision, parameter :: Gmin1 = SpecificOfHeatRatio - 1.0d0
    double precision, parameter :: InverseGmin1 = 1.0d0/Gmin1

    double precision, parameter :: EntropyCorrection = 0.2d0 !FDS用エントロピー補正量

    double precision, parameter :: RungeKuttaParameterB1 = 2.0d0
    double precision, parameter :: RungeKuttaParameterB2 = -1.0d0
    double precision, parameter :: RungeKuttaParameterC2 = 0.5d0

    double precision, parameter :: VenkatakrishnanParameterK = 0.3d0
    integer, target :: iZero=0, iOne=1, iTwo=2

    double precision, parameter :: CourantFriedrichsLewyCondition = 1.0d0

end module ConstantVar_Mod
