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
!-----------------------------------------------------------------------

MODULE forward

#include "include.f90"
#include "macros.f90"

   USE types 
   USE ode
   USE, INTRINSIC :: IEEE_ARITHMETIC

   IMPLICIT NONE 

CONTAINS

   !-----------------------------------------------------------------------
   !> subroutine Forward_constant_stress 
   !  Forward function called from the sampling routine. Performs the 0d forward
   !  calculation and returns the misfit of forward model against the data. 
   !
   ! INPUT & OUTPUT:
   ! @param config - pointer to the config file.
   ! @param mis    - misfit of model vs data. (OUTPUT) 
   ! @param dims   - number of dimensions/parameters to optimize. 
   ! @param m      - the current sample for which forward calculations are needed. 
   !
   ! \author : Sagar Masuti
   !----------------------------------------------------------------------

   SUBROUTINE forward_constant_strain_rate(config,mis,dims,m)

      IMPLICIT NONE 

      TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config
      REAL*8, INTENT(OUT) :: mis
      INTEGER, INTENT(IN) :: dims
      REAL*8, DIMENSION(dims), INTENT(IN) :: m
      
      INTEGER :: i,j,ierr,iostatus
      CHARACTER(256) :: file1
      REAL*8 :: epsilonik,epsilonim
      REAL*8, DIMENSION(STATE_VECTOR_DGF) :: y0

      TYPE(DATA_STRUCT) :: dpre
      TYPE(DATA_STRUCT) :: pre

      mis=0.0_8
      DO i=1,config%n

         IF (config%debug) THEN
            file1="output" // "/" // config%exps(i)%filename
            OPEN (UNIT=15,FILE=file1,IOSTAT=iostatus,FORM="FORMATTED")
            IF (iostatus>0) STOP "could not open point file for writing"
            WRITE (15,'("#strain         stress")') 
         END IF 

         IF (config%exps(i)%instrain .EQ. 0)THEN
            epsilonik=0._8
            epsilonim=0._8
         END IF
         y0(STATE_VECTOR_MAXWELL)=epsilonim
         y0(STATE_VECTOR_STRESS)=config%exps(i)%sigma
         y0(STATE_VECTOR_KELVIN)=epsilonik
         y0(STATE_VECTOR_TSTRAIN)=config%exps(i)%instrain

         CALL ode_solver(config,i,dpre,dims,m,epsilonik,epsilonim,y0)

         !PRINT *, 'n: ',config%d(i)%n
         ALLOCATE(pre%stress(config%d(i)%n),STAT=ierr)
         IF (ierr>0) STOP "could not allocate pre"

         ! interpolate 
         CALL interp1(dpre%n,dpre%strain,dpre%stress,config%d(i)%n, &
                      config%d(i)%strain,pre%stress) 

         IF (.TRUE. .EQV. config%debug) THEN
            PRINT *, config%exps(i)%filename
         END IF

         ! compute misfit
         DO j=1,config%d(i)%n
            IF (.TRUE. .EQV. config%debug) THEN
               WRITE (15,'(2ES10.2E2)') config%d(i)%strain(j), &
                                        pre%stress(j)
               PRINT '(2ES10.2E2)',config%d(i)%strain(j), &
                                   pre%stress(j)
            END IF
            IF (ISNAN(pre%stress(j)) .NEQV. .TRUE.) THEN
               mis=mis+ &
                   REAL(((pre%stress(j)-config%d(i)%stress(j))**2)/(config%d(i)%error(j)**2))
            ELSE
               mis=mis+1e6
            END IF 
         END DO

         DEALLOCATE(pre%stress,dpre%strain,dpre%stress)
         IF (.TRUE. .EQV. config%debug) THEN
            CLOSE(15)
         END IF
      END DO
   END SUBROUTINE 

   !-----------------------------------------------------------------------
   !> subroutine ode_solver 
   !  Calls the ode Runge-Kutta/predictor-corrector solver 
   !  
   ! INPUT & OUTPUT:
   ! @param config    - pointer to the config file.
   ! @param ind       - index of the study for which forward calcution is needed. 
   ! @param dpre      - strain from the solver. (OUTPUT) 
   ! @param dims      - number of dimensions/parameters to optimize. 
   ! @param m         - the current sample for which forward calculations are needed. 
   ! @param epsilonik - previous epsilonik for initialization
   ! @param y0        - initial state vector
   !
   ! \author : Sagar Masuti
   !----------------------------------------------------------------------
   SUBROUTINE ode_solver(config,ind,dpre,dims,m,epsilonik,epsilonim,y0)

      IMPLICIT NONE 

      TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config
      INTEGER, INTENT(IN) :: ind
      TYPE(DATA_STRUCT), INTENT(INOUT) :: dpre 
      INTEGER, INTENT(IN) :: dims
      REAL*8, DIMENSION(dims) :: m
      REAL*8, INTENT(INOUT) :: epsilonik,epsilonim
      REAL*8, DIMENSION(STATE_VECTOR_DGF) :: y0

      ! error flag
      INTEGER :: ierr

      ! state vector
      REAL*8, DIMENSION(STATE_VECTOR_DGF) :: y

      ! rate of change of state vector
      REAL*8, DIMENSION(STATE_VECTOR_DGF) :: dydt,yscal

      ! temporary variables
      REAL*8, DIMENSION(STATE_VECTOR_DGF) :: ytmp,ytmp1,ytmp2,ytmp3

      ! time
      REAL*8 :: time,t0,dt_try,dt_next,dt_done
      
      INTEGER :: i
      
      REAL*8, DIMENSION(:), ALLOCATABLE :: simepsilonik

      ! maximum number of time steps (default)
      INTEGER :: maximumIterations
   

      ! initialize the y vector
      CALL initStateVector(STATE_VECTOR_DGF,y,y0)

      ! start time
      time=0.d0
      ! initial tentative time step
      dt_next=0.5

      maximumIterations=NINT((config%d(ind)%strain(config%d(ind)%n)-config%d(ind)%strain(1))/config%exps(ind)%strain_rate)
      maximumIterations=maximumIterations*2

      dpre%n=maximumIterations
      ALLOCATE(dpre%strain(maximumIterations),  &
               dpre%stress(maximumIterations),  &
               simepsilonik(maximumIterations),STAT=ierr)
      IF (ierr>0) STOP "could not allocate dpre and simepsilonik"

      DO i=1,maximumIterations
         t0=0.d0
         dt_try=dt_next
#ifdef PEC
         dt_done=dt_try
         CALL ode_pec(STATE_VECTOR_DGF,time,y,ytmp1,dt_try,odefun, &
                     config,ind,dims,m)
#else
         CALL odefun(STATE_VECTOR_DGF,time,y,dydt,config,ind,dims,m)

         yscal(:)=abs(y(:))+abs(dt_try*dydt(:))+TINY

         CALL RKQSm(STATE_VECTOR_DGF,t0,y,dydt, &
                    yscal,ytmp1,ytmp2,ytmp3,dt_try,&
                    dt_done,dt_next,odefun,config,ind,dims,m)
         IF (config%interval .LE. time) THEN
            EXIT
         END IF
#endif
         time=time+dt_done
         dpre%strain(i)=y(STATE_VECTOR_TSTRAIN)
         dpre%stress(i)=y(STATE_VECTOR_STRESS)
         simepsilonik(i)=y(STATE_VECTOR_KELVIN)
     END DO

     IF (ind+1 .LE. config%n) THEN
        IF(config%exps(ind+1)%instrain .NE. 0._8) THEN
           IF (config%exps(ind+1)%sigma .LT. (config%d(ind)%stress(config%d(ind)%n)-10)) THEN 
               CALL interp1(i-1,dpre%stress,simepsilonik,1, &
                            config%exps(ind+1)%sigma,epsilonik)     
           ELSE
               epsilonik=simepsilonik(i-1)
           END IF
        ELSE 
           epsilonik=0._8
        END IF
     END IF

     DEALLOCATE(simepsilonik)

   END SUBROUTINE 

   !-----------------------------------------------------------------------
   !> subroutine initStateVector
   ! initialize the state vector
   !
   ! INPUT:
   ! @param n - number of state elements own by current thread
   ! @param y - the state vector (segment owned by currect thread)
   !
   !----------------------------------------------------------------------
   SUBROUTINE initStateVector(n,y,y0)

     IMPLICIT NONE

     INTEGER, INTENT(IN)    :: n
     REAL*8, INTENT(OUT)    :: y(n)
     REAL*8, INTENT(IN)    :: y0(n)

     y(1)=y0(1)
     y(2)=y0(2)
     y(3)=y0(3)
     y(4)=y0(4)

   END SUBROUTINE initStateVector

   !-----------------------------------------------------------------------
   !> subroutine odefun 
   !  Implements the 0d Nonlinear Burgers rheology
   !
   ! \author : Sagar Masuti
   !----------------------------------------------------------------------

   SUBROUTINE odefun(n,time,y,dydt,config,ind,dims,m)

      IMPLICIT NONE 
  
      INTEGER, INTENT(IN)   :: n
      REAL*8, INTENT(IN)    :: time
      REAL*8, INTENT(IN)    :: y(n)
      REAL*8, INTENT(INOUT) :: dydt(n)
      TYPE(CONFIG_STRUCT), INTENT(IN) :: config
      INTEGER, INTENT(IN)   :: ind
      INTEGER, INTENT(IN) :: dims
      REAL*8, DIMENSION(dims) :: m
      REAL*8 :: c1,c2,p,sigma,T,R,epsdot,G

      REAL*8 :: Ak,Ek,Gk,nk,A,nexp,E
      REAL*8 :: epsilonidot,epsilonik,epsilonikdot,epsilonimdot
      
      p=config%exps(ind)%pressure
      T=config%exps(ind)%T
      epsdot=config%exps(ind)%strain_rate
      R=8.3144

      Ak=m(1)
      Gk=m(2)
      nk=m(3)
      Ek=m(4)
      G=m(5)
      IF (config%ss_a .GT. 0) THEN
         A=10**config%ss_a
      ELSE
         A=m(6)
      END IF

      IF (config%ss_q .GT. 0) THEN
         E=config%ss_q
      ELSE
         E=m(7)
      END IF
      
      IF (config%ss_n .GT. 0) THEN
         nexp=config%ss_n
      ELSE
         nexp=m(8)
      END IF

      c1=exp(-(E)/(R*T));
      c2=exp(-(Ek)/(R*T));

      ! Implements the masuti et al, 2016 rheology
      sigma=y(STATE_VECTOR_STRESS)
      epsilonik=y(STATE_VECTOR_KELVIN)
      epsilonimdot=A*(sigma*(abs(sigma)**(nexp-1)))*c1
      epsilonikdot=Ak*((sigma - Gk*epsilonik)*(abs(sigma - Gk*epsilonik))**(nk-1))*c2
 
      dydt(STATE_VECTOR_STRESS)=G*(epsdot-(epsilonimdot+epsilonikdot))
      dydt(STATE_VECTOR_MAXWELL)=epsilonimdot
      dydt(STATE_VECTOR_KELVIN)=epsilonikdot
      dydt(STATE_VECTOR_TSTRAIN)=epsdot
 
   END SUBROUTINE 
END MODULE 
