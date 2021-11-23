!-----------------------------------------------------------------------
! Copyright 2018 Sagar Masuti 
!
! This file is part of RHEOLOGY 
!
! RHEOLOGY is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! RHEOLOGY is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with RHEOLOGY.  If not, see <http://www.gnu.org/licenses/>.
!
! \author   Sylvain Barbot (2017)
! \author   Sagar Masuti
!
! 01-10-2019: Added predictor-corrector method 
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE ode

  USE types

  IMPLICIT NONE

  PUBLIC
  REAL(8), ALLOCATABLE :: AK2(:),AK3(:),AK4(:),AK5(:),AK6(:)
  REAL(8), ALLOCATABLE :: yrkck(:)

  ! numerical accuracy
  REAL(8) :: epsilon = 1.d-6

  ! maximum time step
  REAL*8 :: maximumTimeStep = 1.7976931348623158d308

CONTAINS

  !------------------------------------------------------------------------------
  !       FIFTH-ORDER RUNGE-KUTTA  : STEPPER ROUTINE
  !       SEE NUMERICAL RECIPES 2ND. ED. P.712
  !------------------------------------------------------------------------------
  SUBROUTINE RKQSm(n,t,y,dydt,yscal,ycand,yerr,ytmp,htry,hdid,hnext,odefun,config,ind,dims,m)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(IN)    :: htry,yscal(n)
    REAL*8, INTENT(INOUT) :: ytmp(n)
    REAL*8, INTENT(INOUT) :: t
    REAL*8, INTENT(INOUT) :: dydt(n),y(n)
    REAL*8, INTENT(INOUT) :: ycand(n),yerr(n)
    REAL*8, INTENT(INOUT) :: hdid,hnext
    EXTERNAL              :: odefun
    TYPE(CONFIG_STRUCT), INTENT(IN) :: config
    INTEGER, INTENT(IN)   :: ind
    INTEGER, INTENT(IN) :: dims
    REAL*8, DIMENSION(dims) :: m

    REAL*8, PARAMETER    :: SAFETY=0.7d0,PGROW=-0.2d0,PSHRNK=-0.15d0,ERRCON=1.89d-4
    INTEGER              :: ierr
    REAL*8               :: h,errmax,errmax_thrd,tnew

    h=htry

    DO WHILE (.TRUE.)
      
       CALL RKCKm(n,t,y,dydt,h,ycand,yerr,ytmp,odefun,config,ind,dims,m)
       errmax =MAXVAL(ABS(yerr(1:n)/yscal(1:n)))
       errmax = errmax/epsilon

       IF (errmax .GT. 1.d0) THEN
          h=h*SAFETY*(errmax**PSHRNK)
          IF (h .LT. 0.5d0*h) THEN
             h=0.5d0*h
          ENDIF
          tnew=t+h
          IF (tnew .EQ. t) THEN
             PRINT *, m
             WRITE(STDERR,'("# stepsize underflow in RKQS")')
             STOP 1
          END IF
       ELSE
          IF (errmax .GT. ERRCON) THEN
             hnext=SAFETY*h*(errmax**PGROW)
          ELSE
             hnext=1.d1*h
          END IF
          IF (hnext .GT. maximumTimeStep) THEN
             hnext=maximumTimeStep
          END IF

          hdid=h
          t=t+h
          y(1:n)=ycand(1:n)
          EXIT 
       ENDIF
    END DO

  END SUBROUTINE RKQSm

  !------------------------------------------------------------------------------
  !       FIFTH-ORDER RUNGE-KUTTA  : ALGORITHM ROUTINE
  !       SEE NUMERICAL RECIPES 2ND. ED. P.713
  !
  ! uses the BLAS library
  !------------------------------------------------------------------------------
  SUBROUTINE RKCKm(n,t,yin,dydt,h,yout,yerr,ytmp,odefun,config,ind,dims,m)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(IN)    :: t,h
    REAL*8, INTENT(INOUT) :: yin(n),dydt(n), &
                             yout(n),yerr(n),ytmp(n)
    EXTERNAL :: odefun
    TYPE(CONFIG_STRUCT), INTENT(IN) :: config
    INTEGER, INTENT(IN)   :: ind
    INTEGER, INTENT(IN) :: dims
    REAL*8, DIMENSION(dims) :: m

    REAL*8, PARAMETER :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.d0, A6=0.875d0
    REAL*8, PARAMETER :: B21=0.2d0, B31=3.d0/40.d0, B32=9.d0/40.d0
    REAL*8, PARAMETER :: B41=0.3d0, B42=-0.9d0, B43=1.2d0
    REAL*8, PARAMETER :: B51=-11.d0/54.d0, B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0  
    REAL*8, PARAMETER :: B61=1631.d0/55296.d0,   B62=175.d0/512.d0, B63=575.d0/13824.d0, &
                         B64=44275.d0/110592.d0, B65=253.d0/4096.d0
    REAL*8, PARAMETER :: C1=37.d0/378.d0, C3=250.d0/621.d0, C4=125.d0/594.d0, C6=512.d0/1771.d0
    REAL*8, PARAMETER :: DC1=C1-2825.d0/27648.d0,  DC3=C3-18575.d0/48384.d0, &
                         DC4=C4-13525.d0/55296.d0, DC5=-277.d0/14336.d0, DC6=C6-0.25d0

    ytmp = yin + h*B21*dydt
    CALL odefun(n,t+A2*h,ytmp,AK2,config,ind,dims,m)

    ytmp = yin + h*(B31*dydt + B32*AK2)
    CALL odefun(n,t+A3*h,ytmp,AK3,config,ind,dims,m)

    ytmp = yin + h*(B41*dydt + B42*AK2 + B43*AK3)
    CALL odefun(n,t+A4*h,ytmp,AK4,config,ind,dims,m)

    ytmp = yin + h*(B51*dydt + B52*AK2 + B53*AK3 + B54*AK4)
    CALL odefun(n,t+A5*h,ytmp,AK5,config,ind,dims,m)

    ytmp = yin + h*(B61*dydt + B62*AK2 + B63*AK3 + B64*AK4 + B65*AK5)
    CALL odefun(n,t+A6*h,ytmp,AK6,config,ind,dims,m)

    yout = yin + h*(C1*dydt + C3*AK3 + C4*AK4 + C6*AK6)
    yerr = h*(DC1*dydt + DC3*AK3 + DC4*AK4 + DC5*AK5 + DC6*AK6)

  END SUBROUTINE RKCKm
  
  !------------------------------------------------------------------------------
  ! Predictor-corrector method
  !------------------------------------------------------------------------------
  SUBROUTINE ode_pec(n,t,y,ytmp,htry,odefun,config,ind,dims,m)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(INOUT) :: ytmp(n)
    REAL*8, INTENT(INOUT) :: t
    REAL*8, INTENT(INOUT) :: y(n)
    EXTERNAL              :: odefun
    TYPE(CONFIG_STRUCT), INTENT(IN) :: config
    REAL*8, INTENT(IN)    :: htry
    INTEGER, INTENT(IN)   :: ind
    INTEGER, INTENT(IN)   :: dims
    REAL*8, DIMENSION(dims) :: m
    REAL*8               :: h

    h=htry
     
    CALL odefun(n,t,y,AK2,config,ind,dims,m)
    ytmp = y + h*AK2

    t = t+h

    CALL odefun(n,t,ytmp,AK3,config,ind,dims,m)
    y = y + 0.5*h*(AK2 + AK3)

  END SUBROUTINE ode_pec 

END MODULE ode
