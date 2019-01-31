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
    integer :: CheckNaNInterval = 500 !NaNになってないか確認する間隔
    integer :: CheckNaNInterval_default = 500 !NaNになってないか確認する間隔
    integer :: IterationNumber !計算の反復回数
    integer :: CoreNumberOfCPU = 1
    integer :: GridNumber = 1!重合格子法にて使用する格子の数
    integer :: RetryFlag = 0 !音速が負になる等の問題が発生するとこのフラグが点灯し，再計算を実施する．フラグが点灯した状態で再度フラグ点灯しようとすると強制終了させる
    integer :: DetailedReport = 0   ! 0で画面出力なし，1で開始・終了報告のみ，2でリトライの通知を追加，3で反復回数の通知を追加，4で残差通知を追加

    double precision :: dAttackAngle = 0.0d0

    !Constant value for Calculate
    logical :: invicid = .true.
    double precision, parameter :: Gmin1 = SpecificOfHeatRatio - 1.0d0
    double precision, parameter :: InverseGmin1 = 1.0d0/Gmin1
    complex(kind(0d0)), parameter :: ImaginaryNumber = cmplx(0.0d0, 1.0d0, kind(0d0))
    double precision :: Converge_tolerance = 10.0d0**(-5.0d0)

    ! CFL Number
    double precision :: CourantFriedrichsLewyCondition = 1.0d0
    double precision :: CFL_default = 0.5d0

    ! Physical Property Value
    double precision :: ReynoldsNumber = 0.0d0
    double precision :: MachNumber = 1.0d0
    double precision, parameter :: LaminarPrandtlNumber = 0.72d0
    double precision, parameter :: TurbulentPrandtlNumber = 0.9d0

    ! Dimensional Parameter (高度10km前後の遷音速巡行を想定)
    double precision :: InfinityDensity = 0.4135d0    ! [kg/m3]
    double precision :: InfinityPressure = 2.650d0 * (10.0d0 ** 4)  ! [Pa]
    double precision :: InfinityTemperature = 223.25d0 ! [K]
    double precision :: InfinityStaticViscosity = 2.6566d0 * (10.0d0 ** (-5))
    double precision :: InfinityKineticViscosity = 6.42469d0 * (10.0d0 ** (-5))

    ! Roe's FDS parameter
    double precision :: EntropyCorrection = 0.3d0 !FDS用エントロピー補正量

    ! Rational Runge-Kutta 2 parameter
    double precision, parameter :: RungeKuttaParameterB1 = 2.0d0
    double precision, parameter :: RungeKuttaParameterB2 = -1.0d0
    double precision, parameter :: RungeKuttaParameterC2 = 0.5d0

    ! Venkatakrishnan Limiter parameter
    double precision :: VenkatakrishnanParameterK = 0.3d0

    ! Sutherland's Law parameter (for gas)
    double precision, parameter :: SutherlandTemperature_S = 110.4d0
    double precision, parameter :: ReferenceTemperature_Tref = 273.15d0
    double precision, parameter :: ReferenceViscosity_Mu0 = 0.00001716d0
    double precision, parameter :: SutherlandCoefficient1 = 0.000001458d0 ! Mu0 / (Tref)**(3/2) * (Tref + S)

    ! Baldwin-Lomax parameter
    double precision, parameter :: KarmanConstant = 0.4d0   ! eq.(6)
    double precision, parameter :: ClauserConstant = 0.0168d0   ! eq.(7)
    double precision, parameter :: C_mutm = 14.0d0  ! eq.(11)

    integer, target :: iZero=0, iOne=1, iTwo=2
    ! Influence Area Splitting Method (Toyota, 2018)
    double precision :: SigmoidGain
    double precision :: SigmoidRange = 34.538776394910684d0 !=log((1.0d0 - 10.0d0**(-15))/10.0d0**(-15))

end module ConstantVar_Mod
