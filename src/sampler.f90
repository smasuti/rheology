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

MODULE sampler

#include "include.f90"
#include "macros.f90"
   
   USE types
   USE forward

   IMPLICIT NONE 
   
#include "mpif.h"

CONTAINS 
   !-----------------------------------------------------------------------
   !> subroutine Gibbs 
   !  Implements the Gibbs sampling method with specified conditional evaluation
   !  and calls the user specified forward subroutine.
   !
   ! INPUT:
   ! @param config - pointer to the config file.
   !  
   ! OUTPUT: 
   ! @param samples and its misfit are output to standard display 
   !
   ! \author : Sagar Masuti
   !----------------------------------------------------------------------
   SUBROUTINE gibbs(config)

      IMPLICIT NONE 
     
      TYPE(CONFIG_STRUCT), INTENT(INOUT) :: config

      INTEGER :: i,j,k,ierr,ind,iseed,t,nd,ncond,nsamples,nproc
      REAL*8 :: a,b,ms(config%sample_param%nd),chisquare,mcond

      REAL*8, DIMENSION(:), ALLOCATABLE :: residual 
      REAL*8, DIMENSION(:), ALLOCATABLE :: ncondsamples
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: bounds
      REAL*8 :: m(config%sample_param%nd)

#ifdef USE_MPI
      INTEGER :: ii,iproc
      REAL*8, DIMENSION(:), ALLOCATABLE :: residualproc 
      REAL*8, DIMENSION(:), ALLOCATABLE :: ncondsamplesproc

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

      ALLOCATE(residualproc(nproc),STAT=ierr)
      IF (ierr>0) STOP "could not allocate residualproc"
   
      ALLOCATE(ncondsamplesproc(nproc),STAT=ierr)
      IF (ierr>0) STOP "could not allocate ncondsamplesproc"

      iseed=12345+iproc
#else
      nproc=1
      iseed=12345
#endif
      config%iseed=iseed
      nd=config%sample_param%nd
      ncond=config%sample_param%ncond
      nsamples=config%sample_param%nsamples
      bounds=config%sample_param%bounds
      
      ALLOCATE(residual(ncond),STAT=ierr)
      IF (ierr>0) STOP "could not allocate residual"
   
      ALLOCATE(ncondsamples(ncond),STAT=ierr)
      IF (ierr>0) STOP "could not allocate ncondsamples"

      DO j=1,nd
         m(j)=10**config%sample_param%m0(j)
      END DO

      DO i=1,nsamples
         DO j=1,nd
            DO k=1,ncond,nproc
               a=ran3(iseed)
               b=1-a
               m(j)=10**(b*bounds(j,1)+a*bounds(j,2))
               CALL forward_constant_strain_rate(config,chisquare,nd,m) 
#ifdef USE_MPI
               CALL MPI_ALLGATHER(chisquare,1,MPI_REAL8,residualproc,1,MPI_REAL8, &
                                  MPI_COMM_WORLD,ierr)
               CALL MPI_ALLGATHER(m(j),1,MPI_REAL8,ncondsamplesproc,1,MPI_REAL8, &
                                  MPI_COMM_WORLD,ierr)
                
               DO ii=1,nproc
                  residual(k+(ii-1))=residualproc(ii) 
                  ncondsamples(k+(ii-1))=ncondsamplesproc(ii) 
               END DO
#else
               residual(k)=chisquare
               ncondsamples(k)=m(j)
#endif
            END DO 
            ind=MINLOC(residual,1) 
            m(j)=ncondsamples(ind)
         END DO
#ifdef USE_MPI
         IF (0 .EQ. iproc) THEN
            CALL forward_constant_strain_rate(config,chisquare,nd,m)
#else
            CALL forward_constant_strain_rate(config,chisquare,nd,m)
#endif
         PRINT '("Sample and misfit: ",I7," ",6ES10.3E2)',i,m,chisquare
         CALL FLUSH(STDOUT)
#ifdef USE_MPI
         END IF
#endif 
      END DO
      
      DEALLOCATE(residual,ncondsamples)
#ifdef USE_MPI
      DEALLOCATE(residualproc,ncondsamplesproc)
#endif
   END SUBROUTINE gibbs 

   !-----------------------------------------------------------------------
   !> subroutine ran3 
   !-----------------------------------------------------------------------
   FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      IF (idum .LT. 0.or.iff.eq.0) THEN
       iff=1
       mj=MSEED-iabs(idum)
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       DO 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          IF (mk .LT. MZ)mk=mk+MBIG
          mj=ma(ii)
      11     CONTINUE
       DO 13 k=1,4
          DO 12 i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             IF (ma(i) .LT. MZ)ma(i)=ma(i)+MBIG
      12        CONTINUE
      13     CONTINUE
       inext=0
       inextp=31
       idum=1
      END IF
      inext=inext+1
      IF (inext.eq.56)inext=1
      inextp=inextp+1
      IF (inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      IF (mj .LT. MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      RETURN
   END FUNCTION ran3 

END MODULE sampler
