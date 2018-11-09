module ExactSol
implicit none

!///////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////

! shock tube problem - exact solution

!//////////////////////////////////////////////////////////
!//////////////////////////////////////////////////////////
!任意の時間，位置，速度，圧力，密度を格納するための配列を入力すれば使える
!空間刻みは密度用の配列の上限から取るため，配列は1:imaxで確保しておくこと
contains
    subroutine ExactSolver(t,x,u,p,rho)
    implicit none

    double precision, intent(in) :: t
    double precision, intent(inout) :: x(:),u(:),p(:),rho(:)

    double precision :: gamma
    double precision :: rho_l, rho_r, u_l, u_r, p_l, p_r, x_l, x_r, x_0, c_l, c_r
    double precision :: xi

    integer :: imx, i
    double precision :: dltx

    double precision :: gamm1, gamp1, gam2, gamg2, gamg2i
    double precision :: ratc, ratpre
    double precision :: coe1, coe2, coe3

    double precision :: pp, pp1, pp2, res, res2
    double precision :: beta, u_2, v_sh, v_cd
    double precision :: bcr, bc12, bc23, bc34, bc45, bcl, c_4

!----------------------------------------------------------
! initial condition 1
!----------------------------------------------------------

      gamma = 1.4d0
      rho_l = 1.0d0
      rho_r = 0.1d0
      u_l   = 0.0d0
      u_r   = 0.0d0
      p_l   = rho_l/gamma
      p_r   = rho_r/gamma

      !imx   = 501                  !格子点数
      imx = ubound(rho,1) - lbound(rho,1) + 1
!      t     = 0.50d0               !計算する時間
      !x_l   = -1.0d0                !左側境界
      x_l   = 0.0d0                !左側境界
      x_r   = 1.0d0                !右側境界
      x_0   = (x_r+x_l)*0.5d0      !膜の位置


!----------------------------------------------------------
! initial condition 2
!----------------------------------------------------------

      c_l=dsqrt(gamma*p_l/rho_l)   !音速
      c_r=dsqrt(gamma*p_r/rho_r)   !音速

!----------------------------------------------------------
! make grid
!----------------------------------------------------------

      dltx=(x_r-x_l)/dble(imx-1)   !格子幅
      do i=1,imx
       x(i)=x_l+dltx*dble(i-1)
      end do
!
!----------------------------------------------------------
! output initial condition
!----------------------------------------------------------
!
      do i=1,imx
       xi=x(i)
!
! --- region 1 ---(低圧部)
       if ((xi<=x_r).and.(xi>=x_0)) then
        u(i)=u_r
        p(i)=p_r
        rho(i)=rho_r
! --- region 5 ---(高圧部)
       else if ((xi<x_0).and.(xi>=x_l)) then
        u(i)=u_l
        p(i)=p_l
        rho(i)=rho_l
! --- otherwise ---
       else
        u(i)=0.0d0
        p(i)=0.0d0
        rho(i)=0.0d0
       end if
!
!d       write(10,'(1P,2E25.15)') x(i),u(i)
!d	 write(11,'(1P,2E25.15)') x(i),p(i)
!d	 write(12,'(1P,2E25.15)') x(i),rho(i)
      end do
!
! for gnuplot recognizing
!
!d      write(10,*)
!d	write(11,*)
!d	write(12,*)
!
!----------------------------------------------------------
! bisection method (二分法)
!----------------------------------------------------------
!
      gamm1 = gamma-1.0d0
      gamp1 = gamma+1.0d0
      gam2  = 2.0d0/gamm1
      gamg2 = gam2*gamma
      gamg2i= 1.0d0/gamg2
      ratc= c_l/c_r            !音速比
      ratpre= p_r/p_l          !圧力比
      coe1  = gam2*ratc
      coe2  = (u_l-u_r)/c_r
      coe3  = dsqrt(2.0d0/gamma)
!
      pp1=0.1d0*p_l/p_r
      pp2=p_l/p_r
!
      do i=1,10000
       pp=(pp1+pp2)*0.5d0
       res = func(pp)
       res2 = func(pp2)
    if(res*res2.le.0.0d0) then
        pp1=pp
       else
        pp2=pp
       end if
      end do
!
!----------------------------------------------------------
! beta, speed of shockwave and contact discontinuity
!----------------------------------------------------------
!
      beta  = dsqrt(gamm1+(gamp1)*pp)
      u_2   = u_r+c_r*coe3*(pp-1.0d0)/beta              !領域2での速度
      v_sh  = u_r+(pp-1.0d0)*(c_r**2)/(gamma*(u_2-u_r)) !衝撃波速度
      v_cd  = u_2                                       !接触速度
!
!----------------------------------------------------------
! output state
!----------------------------------------------------------
!
!	do j=5,25,5
!	 t=j/100.0d0
!	t=0.25d0
!
      bcr   = x_r
      bc12  = x_0+v_sh*t
      bc23  = x_0+v_cd*t
      bc34  = x_0+((gamp1*v_cd*0.5d0)-c_l-(gamm1*u_l*0.5d0))*t
      bc45  = x_0+(u_l-c_l)*t
      bcl   = x_l
!
      do i=1,imx
       xi=x(i)
!
! ---over left---
!
       if (xi>bcr) then
        u(i)=0.0d0
        p(i)=0.0d0
        rho(i)=0.0d0
!
!--- region 1 ---(u_1,p_1,rho_1)(低圧部)
!
       else if (xi>=bc12) then
        u(i)=u_r
        p(i)=p_r
        rho(i)=rho_r
!
!--- region 2 ---(u_2,p_2,rho_2) (衝撃波－接触不連続)
!
       else if (xi>=bc23) then
        u(i)=u_2
        p(i)=p_r*pp
        rho(i)=rho_r*(gamm1+(gamp1)*pp)/(gamp1+(gamm1)*pp)
!	  c0=c_r+gamm1*(u(i)-u_r)*0.5d0
!	  c1=dsqrt(gamma*p(i)/rho(i))
!	  write(6,*) u(i),p(i),rho(i)
!
!--- region 3 ---(u_3,p_3,rho_3) (接触不連続－膨張波)
!
       else if (xi>=bc34) then
        u(i)=u_2
        p(i)=p_r*pp
        rho(i)=rho_l*((p(i)/p_l)**(1.0d0/gamma))
!	  write(6,*) u(i),p(i),rho(i)
!
!--- region 4 ---(u_4,p_4,rho_4) (膨張波)
!
       else if (xi>=bc45) then
        u(i)=(2.0d0/gamp1) &
     &   *(((x(i)-x_0)/t)+c_l+(gamm1/2.0d0)*u_l)
        c_4=c_l-((gamm1)/2.0d0)*(u(i)-u_l)
        p(i)=p_l*((c_4/c_l)**gamg2)
        rho(i)=rho_l*((p(i)/p_l)**(1.0d0/gamma))
!	  write(6,*) u(i),p(i),rho(i)
!
!--- region 5 ---(u_5,p_5,rho_5) (高圧部)
!
       else if (xi>=bcl) then
        u(i)=u_l
        p(i)=p_l
        rho(i)=rho_l
!
! ---over left---
!
       else
        u(i)=0.0d0
        p(i)=0.0d0
        rho(i)=0.0d0
       end if
!
    end do
!
    return
    contains
!
!------------------------------------------------
! function for bisection algorithm
!------------------------------------------------
!
      function func(pp) result(res)
        double precision, intent(in) :: pp
        double precision :: res
!
      res=coe1*(1.0d0-(pp*ratpre)**(gamg2i)) &
     &   +coe2 &
     &   -coe3*(pp-1.0d0)/dsqrt(gamm1+(gamp1)*pp)
!
        return
      end function func

    end subroutine ExactSolver

end module ExactSol
