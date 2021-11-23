!-----------------------------------------------------------------------
! Copyright 2018, Sagar Masuti 
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
! \author : Sylvain Barbot & Sagar Masuti
!-----------------------------------------------------------------------
SUBROUTINE interp1(m,x,y,n,xi,yi)
  INTEGER, INTENT(IN) :: n,m
  REAL*8, DIMENSION(m), INTENT(IN) :: x,y
  REAL*8, DIMENSION(n), INTENT(IN) :: xi
  REAL*8, DIMENSION(n), INTENT(OUT) :: yi

  INTEGER :: i,j

  DO i=1,n
     ! check for end points
     IF (xi(i) .LE. x(1)) THEN
        yi(i)=y(1)
        CYCLE
     END IF
     IF (xi(i) .GE. x(m)) THEN
        yi(i)=y(m)
        CYCLE
     END IF

     IF (2 .EQ. m) THEN
        ! linear interpolation
        yi(i)=y(1)+(y(2)-y(1))/(x(2)-x(1))*(xi(i)-x(1))
        CYCLE
     END IF

     ! left quadratic interpolation
     IF (xi(i) .GT. x(m-1) .AND. xi(i) .LT. x(m)) THEN
        yi(i)=((xi(i)-x(m-1))*(xi(i)-x(m-0)))/((x(m-2)-x(m-1))*(x(m-2)-x(m-0)))*y(m-2) + &
              ((xi(i)-x(m-2))*(xi(i)-x(m-0)))/((x(m-1)-x(m-2))*(x(m-1)-x(m-0)))*y(m-1) + &
              ((xi(i)-x(m-2))*(xi(i)-x(m-1)))/((x(m-0)-x(m-2))*(x(m-0)-x(m-1)))*y(m-0)
        CYCLE
     END IF

     ! right quadratic interpolation
     DO j=1,m-2
        ! check if interpolation is needed
        IF (xi(i).EQ.x(j)) THEN
           yi(i)=y(j)
           CYCLE
        ELSE IF (xi(i).EQ.x(j+1)) THEN
           yi(i)=y(j+1)
           CYCLE
        ELSE IF (xi(i) .EQ. x(j+2)) THEN
           yi(i)=y(j+2)
           CYCLE
        END IF

        ! quadratic interpolation
        IF (xi(i) .GT. x(j) .AND. xi(i) .LT. x(j+1)) THEN
           yi(i)=((xi(i)-x(j+1))*(xi(i)-x(j+2)))/((x(j+0)-x(j+1))*(x(j+0)-x(j+2)))*y(j+0) + &
                 ((xi(i)-x(j+0))*(xi(i)-x(j+2)))/((x(j+1)-x(j+0))*(x(j+1)-x(j+2)))*y(j+1) + &
                 ((xi(i)-x(j+0))*(xi(i)-x(j+1)))/((x(j+2)-x(j+0))*(x(j+2)-x(j+1)))*y(j+2)
        END IF
     END DO
  END DO

END SUBROUTINE interp1

