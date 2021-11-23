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

#ifndef RHEOLOGY_MACROS

#define RHEOLOGY_MACROS 1

#define WRITE_DEBUG_INFO(e) WRITE (0,'("error ",I3.3," at line ",I5.5,", rank ",I4.4)') e,__LINE__,rank
#define WRITE_DEBUG_INFO_SERIAL(e) WRITE (0,'("error ",I3.3," at line ",I5.5)') e,__LINE__

! for the Runge-Kutta method
#define TINY 1.d-30

! constant variables
#define ZERO 0.000000000000000000000000000000000000000000000000d0
#define  ONE 1.000000000000000000000000000000000000000000000000d0

! standard output
#define STDOUT 6

! standard error
#define STDERR 0

! time and time step file
#define FPTIME 14

!-----------------------
! index in state vector
!-----------------------
#define STATE_VECTOR_STRESS         1
#define STATE_VECTOR_KELVIN         2
#define STATE_VECTOR_MAXWELL        3
#define STATE_VECTOR_TSTRAIN        4

#define STATE_VECTOR_DGF            4

#endif

