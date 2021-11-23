
PROGRAM bestfit

#include "macros.f90"
#include "include.f90"

   USE heap 
   USE ode
   USE types
   USE forward

   REAL*8 :: chisquare
   INTEGER :: ierr
   REAL*8, DIMENSION(5) :: m
   CHARACTER(256) :: dataline,filename
   INTEGER :: iunit,nd

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
   
   config%debug=.TRUE.
   CALL forward_constant_strain_rate(config,chisquare,nd,m) 
   PRINT '("m: ",5ES9.2E2)',m
   PRINT '("misfit: ",1ES9.2E2)',chisquare
 
   IF (ALLOCATED(config%exps)) DEALLOCATE(config%exps)
   DO i=1,config%n
      IF (ALLOCATED(config%d(i)%strain)) DEALLOCATE(config%d(i)%strain)
      IF (ALLOCATED(config%d(i)%stress)) DEALLOCATE(config%d(i)%stress)
   END DO
   IF (ALLOCATED(config%d)) DEALLOCATE(config%d)

   DEALLOCATE(AK2,AK3,AK4,AK5,AK6)

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
      USE getopt_m

      IMPLICIT NONE

      TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config

      CHARACTER :: ch
      CHARACTER(256) :: dataline,filename,configfile
      INTEGER :: iunit,noptions
      TYPE(OPTION_S) :: opts(2)

      INTEGER :: i,ierr,k,iostatus,nolines
      INTEGER :: dummy

      config%debug=.TRUE.
      
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

      PRINT '("# number of bestfit model parameters to read")'
      CALL getdata(iunit,dataline)
      READ (dataline,*) nd 
      PRINT '(I02.2)', nd

      PRINT '("# m")'
      DO k=1,nd
         CALL getdata(iunit,dataline)
         READ (dataline,*) m(k)
         PRINT '(I02.2,1ES10.2)', k,m(k)
      END DO

      PRINT '("# config file name")'
      CALL getdata(iunit,dataline)
      READ (dataline,'(a)') configfile
      PRINT '(a)', trim(configfile)

      CLOSE(iunit)

      ! Get the number of lines in the file
      CALL fileLineCount(configfile,nolines)
      config%n=nolines
      iunit=newunit()
      OPEN (UNIT=iunit,FILE=trim(configfile),IOSTAT=iostatus)

      ALLOCATE(config%exps(nolines),STAT=iostatus)
   
      PRINT '("#n         T   stress pressure    strain_rate   instrain filename")'
      DO i=1,nolines
         CALL getdata(iunit,dataline)
         READ (dataline,*,IOSTAT=ierr) k,config%exps(i)%T,config%exps(i)%sigma, &
                                       config%exps(i)%pressure,config%exps(i)%filename, &
                                       config%exps(i)%strain_rate,config%exps(i)%instrain

         PRINT '(I2.2," ",5ES9.2E1," ",A40)',i,config%exps(i)%T,config%exps(i)%sigma, &
                                        config%exps(i)%pressure,config%exps(i)%strain_rate, &
                                        config%exps(i)%instrain,config%exps(i)%filename
      END DO

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
      !                               count the number of lines in ASCII file ignoring
      !                               commented lines and white space
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
100        CONTINUE
           IF (char .EQ. '#') THEN
              READ(iunit,'(a)',END=10) line
              i=1
              char=line(1:1)
200           CONTINUE
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
   ! readobs
   !                               Reads observations
   ! Input:
   !                               Input file with the data.
   ! Output:
   !                               Reads the data into actd
   !--------------------------------------------------------------------------!
   SUBROUTINE readobs(config)

     USE types
     IMPLICIT NONE

     TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config

     INTEGER :: iunit,i,j,n,iostatus
     CHARACTER(256) :: filename,dataline
     LOGICAL :: file_exists

     iunit=1000

     ALLOCATE(config%d(config%n),STAT=iostatus)
     IF (iostatus>0) STOP "could not allocate data"

     DO j=1,config%n
        filename="./"//config%exps(j)%filename

        file_exists=.FALSE.
        INQUIRE (FILE=filename,EXIST=file_exists)

        IF (file_exists) THEN
           CALL fileLineCount(filename,n)
           OPEN (UNIT=iunit,FILE=filename,STATUS="old",FORM="FORMATTED")

           ALLOCATE(config%d(j)%strain(n),config%d(j)%stress(n), &
                    config%d(j)%error(n),STAT=iostatus)
           IF (iostatus>0) STOP "could not allocate time series"


           config%d(j)%n=n
           DO i=1,n
              CALL getdata(iunit,dataline)
              READ (dataline,*) config%d(j)%strain(i), &
                                config%d(j)%stress(i)
              config%d(j)%error(i)=2
           END DO
           CLOSE (iunit)

        ELSE
           PRINT *, filename
           STOP "could not find data (stress vs strain) file"
     END IF
   END DO
END SUBROUTINE

END PROGRAM bestfit
