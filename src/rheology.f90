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
! ---------------------------------------------------------------------
!                              RHEOLOGY
! ---------------------------------------------------------------------
!> program Rheology simulates evolution strain for a point source for 
!! olivine using a nonlinear Burgers type rheology. The flow law paramters
!! are obtained using the MCMC/Bayesian technique.
!!
!! \mainpage
!!
!! The time evolution can be evaluated numerically using either the 
!! predictor-corrector method with fixed step size or the 4th/5th order 
!! Runge-Kutta method with adaptive time steps. The state vector is
!! as follows:
!!
!!    / \sigma     \
!!   |  \epsilon_m  |
!!   |  \epsilon_k  |
!!    \ \epsilon   /
!!
!! where \sigma is the total stress in the system, \epsilon_m, \epsilon_k, 
!! and \epsilon are the strains in Maxwell element, Kelvin element, 
!! total strain. 
!! 
!! The MCMC sampling technique used here is Gibbs sampling. The code is 
!! implemented using MPI (Message Passing Interface) and highly scalable. 
!! However, scalability is limited to the number of conditional evaluations 
!! in each dimension. 
!!
!! \author Sagar Masuti (2018).
!----------------------------------------------------------------------

PROGRAM rheology

#include "macros.f90"
#include "include.f90"

   USE ode
   USE types
   USE sampler
   USE forward

   IMPLICIT NONE

   TYPE(CONFIG_STRUCT) :: config 
   REAL*8 :: mis
   INTEGER :: ierr,i
   REAL*8, DIMENSION(:), ALLOCATABLE :: m
   INTEGER :: dims

#ifdef USE_MPI
   INTEGER :: iproc,nproc

   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
#endif

   ! allocate buffer from ode module
   ALLOCATE(AK2(STATE_VECTOR_DGF), &
            AK3(STATE_VECTOR_DGF), &
            AK4(STATE_VECTOR_DGF), &
            AK5(STATE_VECTOR_DGF), &
            AK6(STATE_VECTOR_DGF), &
            yrkck(STATE_VECTOR_DGF), STAT=ierr)
   IF (ierr>0) STOP "could not allocate the AK1-6 work space"

   config%interval=1e5
   CALL read_config(config)
   CALL FLUSH(STDOUT)
   
   ALLOCATE (m(config%sample_param%nd),STAT=ierr) 
   IF (ierr>0) STOP "could not allocate to m"
   
#ifdef USE_MPI
   IF (0 .EQ. iproc) THEN
#endif 
   PRINT '("# starting gibbs sampler ...")'
   CALL FLUSH(STDOUT)
#ifdef USE_MPI
   END IF
#endif 

   CALL gibbs(config)

#ifdef USE_MPI
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr); 
#endif
   ! free memory
   IF (ALLOCATED(config%exps)) DEALLOCATE(config%exps)
   DO i=1,config%n
      IF (ALLOCATED(config%d(i)%strain)) DEALLOCATE(config%d(i)%strain)
      IF (ALLOCATED(config%d(i)%stress)) DEALLOCATE(config%d(i)%stress)
      IF (ALLOCATED(config%d(i)%error)) DEALLOCATE(config%d(i)%error)
   END DO
   IF (ALLOCATED(config%d)) DEALLOCATE(config%d)

   DEALLOCATE(m,AK2,AK3,AK4,AK5,AK6,config%sample_param%m0,config%sample_param%bounds)

#ifdef USE_MPI   
   CALL MPI_FINALIZE(ierr)
#endif

CONTAINS 

   !-----------------------------------------------------------------------
   !> subroutine read_config 
   ! Reads the configuration file in data directory and the associated data 
   !
   ! INPUT:
   ! @param config - pointer to the config file which needs to be filled. 
   !
   !----------------------------------------------------------------------
   SUBROUTINE read_config(config)
   
      USE types
      USE sampler

      IMPLICIT NONE 

      TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config

      CHARACTER :: ch
      CHARACTER(256) :: dataline,filename,configfile
      INTEGER :: iunit,noptions

      INTEGER :: i,ierr,k,iostatus,nolines
      INTEGER :: dummy

#ifdef USE_MPI      
      INTEGER :: nproc,iproc,position
      INTEGER, PARAMETER :: psize=1024
      CHARACTER, DIMENSION(psize) :: packed

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
#endif 

      config%debug=.FALSE. 
 
#ifdef USE_MPI
      IF (0 .EQ. iproc) THEN
#endif
      noptions=0
      IF (noptions .LT. COMMAND_ARGUMENT_COUNT()) THEN
          ! read from input file
          iunit=newunit()
          CALL GET_COMMAND_ARGUMENT(noptions+1,filename)
          OPEN (UNIT=iunit,FILE=filename,IOSTAT=ierr)
      ELSE
       ! get input parameters from standard input
       iunit=5
      END IF      

      PRINT '("# number of parameters to optimize")'
      CALL getdata(iunit,dataline)
      READ (dataline,*) config%sample_param%nd 
      PRINT '(I02.2)', config%sample_param%nd

      ALLOCATE(config%sample_param%bounds(config%sample_param%nd,2), &
               config%sample_param%m0(config%sample_param%nd),STAT=ierr)
      IF (ierr>0) STOP "could not allocate bounds and m0"

      PRINT '("# i      init     lower     upper")'
      DO k=1,config%sample_param%nd
         CALL getdata(iunit,dataline)
         READ (dataline,*) i,config%sample_param%m0(k),       &
                             config%sample_param%bounds(k,1), &
                             config%sample_param%bounds(k,2)
         PRINT '(I3,3ES10.2)', k,config%sample_param%m0(k),      &
                                 config%sample_param%bounds(k,1),&
                                 config%sample_param%bounds(k,2)
      END DO

      PRINT '("# Number of samples")'
      CALL getdata(iunit,dataline)
      READ (dataline,*) config%sample_param%nsamples 
      PRINT '(I7)', config%sample_param%nsamples 

      PRINT '("# Number of conditionals for gibbs sampling")'
      CALL getdata(iunit,dataline)
      READ (dataline,*) config%sample_param%ncond 
      PRINT '(I4)', config%sample_param%ncond 

      PRINT '("# config file name")'
      CALL getdata(iunit,dataline)
      READ (dataline,'(a)') configfile 
      PRINT '(a)', trim(configfile)

      PRINT '("# location of the stress vs. strain data")'
      CALL getdata(iunit,dataline)
      READ (dataline,'(a)') config%datadir 
      PRINT '(a)', trim(config%datadir)

      PRINT '("# steady-state parameters (A, E, and n)")'
      CALL getdata(iunit,dataline)
      READ (dataline,*) config%ss_a,config%ss_q,config%ss_n
      PRINT '(3ES10.2)', config%ss_a,config%ss_q,config%ss_n 

      CLOSE(iunit)
         
      ! Get the number of lines in the file 
      CALL fileLineCount(configfile,nolines)
      config%n=nolines
      iunit=newunit()
      OPEN (UNIT=iunit,FILE=trim(configfile),IOSTAT=iostatus)

      ALLOCATE(config%exps(nolines),STAT=iostatus)
      IF (ierr>0) STOP "could not allocate exps in config"
   
      PRINT '("#n         T   stress pressure  e_d(1/s)    in_e filename")'
      DO i=1,nolines
         CALL getdata(iunit,dataline)
         READ (dataline,*,IOSTAT=ierr) k,config%exps(i)%T,config%exps(i)%sigma, &
                                       config%exps(i)%pressure,config%exps(i)%filename, &
                                       config%exps(i)%strain_rate,config%exps(i)%instrain

         PRINT '(I2.2," ",5ES9.2E1," ",A40)',i,config%exps(i)%T,config%exps(i)%sigma, &
                                        config%exps(i)%pressure,config%exps(i)%strain_rate, &
                                        config%exps(i)%instrain,config%exps(i)%filename
      END DO

#ifdef USE_MPI
         CALL FLUSH(STDOUT)
         position=0
         CALL MPI_PACK(config%interval,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%n,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%ss_a,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%ss_q,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%ss_n,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%sample_param%nd,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%sample_param%ncond,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%sample_param%nsamples,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_PACK(config%datadir,256,MPI_CHAR,packed,psize,position,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

         ! send sampler configuration 
         DO i=1,config%sample_param%nd
            position=0
            CALL MPI_PACK(config%sample_param%m0(i),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%sample_param%bounds(i,1),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%sample_param%bounds(i,2),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
         END DO

         ! send experimental settings
         DO i=1,config%n
            position=0
            CALL MPI_PACK(config%exps(i)%T,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%exps(i)%sigma,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%exps(i)%filename,256,MPI_CHAR,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%exps(i)%pressure,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%exps(i)%strain_rate,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(config%exps(i)%instrain,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
         END DO

      ELSE ! receive in slaves

         position=0
         CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%interval,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%n,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%ss_a,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%ss_q,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%ss_n,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%sample_param%nd,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%sample_param%ncond,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%sample_param%nsamples,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         CALL MPI_UNPACK(packed,psize,position,config%datadir,256,MPI_CHAR,MPI_COMM_WORLD,ierr)


         IF (config%sample_param%nd .gt. 0) THEN
            ALLOCATE(config%sample_param%bounds(config%sample_param%nd,2), &
                     config%sample_param%m0(config%sample_param%nd),STAT=ierr)
            IF (ierr>0) STOP "could not allocate bounds and m0"

            DO i=1,config%sample_param%nd
               position=0
               CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%sample_param%m0(i),1,MPI_REAL8,MPI_COMM_WORLD,ierr) 
               CALL MPI_UNPACK(packed,psize,position,config%sample_param%bounds(i,1),1,MPI_REAL8,MPI_COMM_WORLD,ierr) 
               CALL MPI_UNPACK(packed,psize,position,config%sample_param%bounds(i,2),1,MPI_REAL8,MPI_COMM_WORLD,ierr) 
            END DO
         END IF

         IF (config%n .gt. 0) THEN
            ALLOCATE(config%exps(config%n),STAT=iostatus)
            DO i=1,config%n
               position=0
               CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%T,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%sigma,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%filename,256,MPI_CHAR,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%pressure,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%strain_rate,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
               CALL MPI_UNPACK(packed,psize,position,config%exps(i)%instrain,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
            END DO 
         END IF
      END IF
#endif 
      ! reading stress vs strain data 
      CALL readobs(config)

   END SUBROUTINE

   !---------------------------------------------------------------------
   !> function newunit
   !! This is a simple function to search for an available unit.
   !! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
   !! The UNIT value is returned by the function, and also by the optional
   !! argument. This allows the function to be used directly in an OPEN
   !! statement, and optionally save the result in a local variable.
   !! If no units are available, -1 is returned.
   !---------------------------------------------------------------------
   INTEGER FUNCTION newunit(unit)
     INTEGER, INTENT(OUT), OPTIONAL :: unit

     INTEGER, PARAMETER :: LUN_MIN=10, LUN_MAX=1000
     LOGICAL :: opened
     INTEGER :: lun

     newunit=-1
     DO lun=LUN_MIN,LUN_MAX
       INQUIRE(UNIT=lun,OPENED=opened)
       IF (.NOT. opened) THEN
         newunit=lun
         EXIT
       END IF
     END DO
     IF (PRESENT(unit)) unit=newunit

   END FUNCTION newunit

   !--------------------------------------------------------------------------!
   ! fileLineCount
   ! count the number of lines in ASCII file ignoring
   ! commented lines and white space
   ! Input:
   ! @param    filename            name of ASCII file
   !--------------------------------------------------------------------------!
   SUBROUTINE fileLineCount(filename,count)
     IMPLICIT NONE
     CHARACTER(256), INTENT(IN) :: filename
     INTEGER, INTENT(OUT) :: count

     INTEGER :: iunit,i
     CHARACTER(256) :: line
     CHARACTER :: char

     count=0

     iunit=newunit()
     OPEN (iunit,FILE=filename)

     DO
        char='#'
100  CONTINUE
        IF (char .EQ. '#') THEN
           READ(iunit,'(a)',END=10) line
           i=1
           char=line(1:1)
200     CONTINUE
           IF (char .EQ. ' ') THEN
              i=i+1
              char=line(i:i)
              GOTO 200
           ENDIF
           GOTO 100
        ENDIF
        count=count+1
     END DO

   10 CLOSE (iunit)

   END SUBROUTINE fileLineCount

   !--------------------------------------------------------------------------!
   ! readobs : Reads observations
   ! Input : Input file with the data.
   ! Output: Reads the data into config%d
   !--------------------------------------------------------------------------!
   SUBROUTINE readobs(config)

     USE types
     IMPLICIT NONE

     TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config 

     INTEGER :: iunit,i,j,n,iostatus
     CHARACTER(256) :: filename,dataline
     LOGICAL :: file_exists

#ifdef USE_MPI
     INTEGER :: iproc,nproc

     iunit=newunit()
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
#else
     iunit=1000
#endif

     ALLOCATE(config%d(config%n),STAT=iostatus)
     IF (iostatus>0) STOP "could not allocate data"

     DO j=1,config%n
        filename=TRIM(config%datadir)// "/" //config%exps(j)%filename

        file_exists=.FALSE.
        INQUIRE (FILE=filename,EXIST=file_exists)

        IF (file_exists) THEN
           CALL fileLineCount(filename,n)
           OPEN (UNIT=iunit,FILE=filename,STATUS="old",FORM="FORMATTED")

           ALLOCATE(config%d(j)%strain(n),config%d(j)%stress(n), &
                    config%d(j)%error(n),STAT=iostatus)
           IF (iostatus>0) STOP "could not allocate time series"


           IF (.TRUE. .EQV. config%debug) THEN
              PRINT '(A)',trim(config%exps(j)%filename)
           ENDIF

           config%d(j)%n=n
           DO i=1,n
              CALL getdata(iunit,dataline)
              READ (dataline,*) config%d(j)%strain(i), &
                                config%d(j)%stress(i)
              config%d(j)%error(i)=2
              IF (.TRUE. .EQV. config%debug) THEN
                 PRINT '(3ES10.2E2)',config%d(j)%strain(i),config%d(j)%stress(i),config%d(j)%error(i)
              ENDIF
           END DO
           CLOSE (iunit)

        ELSE
           PRINT *, filename
           STOP "could not find data (stress vs strain) file"
     END IF
   END DO
END SUBROUTINE 

END PROGRAM
