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
module ConstantVar_Mod
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
    integer :: RetryFlag = 0 !音速が負になる等の問題が発生するとこのフラグが点灯し，再計算を実施する．フラグが点灯した状態で再度フラグ点灯しようとすると強制終了させる
!Constant value for Calculate
    double precision, parameter :: Gmin1 = SpecificOfHeatRatio - 1.0d0
    double precision, parameter :: InverseGmin1 = 1.0d0/Gmin1

    double precision :: EntropyCorrection = 0.1d0 !FDS用エントロピー補正量

    double precision, parameter :: RungeKuttaParameterB1 = 2.0d0
    double precision, parameter :: RungeKuttaParameterB2 = -1.0d0
    double precision, parameter :: RungeKuttaParameterC2 = 0.5d0

    double precision :: VenkatakrishnanParameterK = 0.3d0
    integer, target :: iZero=0, iOne=1, iTwo=2

    double precision :: CourantFriedrichsLewyCondition = 1.0d0

    double precision :: SigmoidGain
    double precision :: SigmoidRange = 34.538776394910684d0 !=log((1.0d0 - 10.0d0**(-15))/10.0d0**(-15))

end module ConstantVar_Mod
