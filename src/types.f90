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

#include "macros.f90"

MODULE types

   IMPLICIT NONE

   

   TYPE EXPERIMENT_STRUCT
      ! temperature
      REAL*8 :: T

      ! stress 
      REAL*8 :: sigma 

      ! filename of data 
      CHARACTER(256) :: filename 

      ! pressure
      REAL*8 :: pressure

      ! strain rate of experiment
      REAL*8 :: strain_rate

      ! water content 
      REAL*8 :: coh

      ! initial strain 
      REAL*8 :: instrain 

      ! TEMP
      REAL*8 :: temp1
      
      ! TEMP
      REAL*8 :: temp2

      ! sigmabmax
      REAL*8 :: sigmabmax

   END TYPE EXPERIMENT_STRUCT 
   
   TYPE DATA_STRUCT
      ! Number of data points
      INTEGER :: n

      ! Time for which the experiment/simulation was performed.
      REAL*8, DIMENSION(:), ALLOCATABLE :: time

      ! Measured/computed strain
      REAL*8, DIMENSION(:), ALLOCATABLE :: strain 

      REAL*8, DIMENSION(:), ALLOCATABLE :: stress 

      REAL*8, DIMENSION(:), ALLOCATABLE :: error 
      
   END TYPE DATA_STRUCT

   TYPE SAMPLER_STRUCT

      ! number of dimensions 
      INTEGER :: nd

      ! Initial value
      REAL*8, DIMENSION(:), ALLOCATABLE :: m0

      ! number of samples to be generated 
      INTEGER :: nsamples 

      ! conditionals to be evaluated in each dimension
      INTEGER :: ncond

      ! bounds on the model parameters
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: bounds
   END TYPE SAMPLER_STRUCT

   TYPE CONFIG_STRUCT 
     INTEGER :: n
 
     ! simulation time
     REAL*8 :: interval

     ! random seed
     INTEGER :: iseed

     ! covariance matrix flag
     INTEGER :: cm_flag

     ! covariance matrix
     REAL*8 :: cm(6)

     TYPE(EXPERIMENT_STRUCT), DIMENSION(:), ALLOCATABLE :: exps 

     TYPE(DATA_STRUCT), DIMENSION(:), ALLOCATABLE :: d

     TYPE(SAMPLER_STRUCT) :: sample_param
 
     ! other options
     LOGICAL :: isdryrun=.FALSE.
     LOGICAL :: ishelp=.FALSE.
     LOGICAL :: isversion=.FALSE.
     LOGICAL :: debug=.FALSE.

  END TYPE CONFIG_STRUCT 

END MODULE types

