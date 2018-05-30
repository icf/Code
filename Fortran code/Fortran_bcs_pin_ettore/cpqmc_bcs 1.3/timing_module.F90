!    Copyright (C) 2007, 2008, 2009  D. Schirmer, M. L. Wall
!    This file is part of OpenSourceTEBD.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
MODULE timing_module
!
!Purpose :Contains routines used to time a simulation
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   2/24/09   M. L. Wall	v1.0 release
!

IMPLICIT NONE

!Global variables required to keep track of initial time
CHARACTER(LEN=8) :: init_date
CHARACTER(LEN=10) :: init_time 

CONTAINS

SUBROUTINE BeginTiming()
!
!Purpose: Begin the timing routine
!
IMPLICIT NONE

CALL DATE_AND_TIME(init_date, init_time)
END SUBROUTINE BeginTiming

!SUBROUTINE EndTiming(length, fileID)
SUBROUTINE EndTiming(fileID)
!
!Purpose: End the timing routine
!
IMPLICIT NONE
#ifdef MPI
include "mpif.h"
#endif
CHARACTER(LEN=8) :: end_date
CHARACTER(LEN=10)  :: end_time
!INTEGER, INTENT(OUT) :: length
INTEGER  :: length
INTEGER, INTENT(IN), OPTIONAL :: fileID
INTEGER :: days, hours, minutes, seconds
INTEGER :: end_time_comb, end_date_comb, end_year, end_month, end_day, end_hour, end_minute, end_second
INTEGER :: init_time_comb, init_date_comb, init_year, init_month, init_day, init_hour, init_minute, init_second

integer::ierr,rank,Nsize

#ifdef MPI
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,Nsize,ierr)
#else
 rank=0
 Nsize=1
#endif


length =0

CALL DATE_AND_TIME(end_date, end_time)
READ(end_date,'(I8)') end_date_comb
READ(end_time(1:6),'(I6)') end_time_comb
READ(init_date,'(I8)') init_date_comb
READ(init_time(1:6),'(I6)') init_time_comb

end_year = end_date_comb/10000
end_month = MOD(end_date_comb,10000)/100
end_day = MOD(end_date_comb,100)
end_hour = end_time_comb/10000
end_minute = MOD(end_time_comb,10000)/100
end_second = MOD(end_time_comb,100)

init_year = init_date_comb/10000
init_month = MOD(init_date_comb,10000)/100
init_day = MOD(init_date_comb,100)
init_hour = init_time_comb/10000
init_minute = MOD(init_time_comb,10000)/100
init_second = MOD(init_time_comb,100)

!The following are only needed for either:
! (1) Test purposes
! (2) Extremely detailed report

!PRINT *, 'init_year = ', init_year
!PRINT *, 'init_month = ', init_month
!PRINT *, 'init_day = ', init_day
!PRINT *, 'init_hour = ', init_hour
!PRINT *, 'init_minute = ', init_minute
!PRINT *, 'init_second = ', init_second

!PRINT *, 'end_year = ', end_year
!PRINT *, 'end_month = ', end_month
!PRINT *, 'end_day = ', end_day
!PRINT *, 'end_hour = ', end_hour
!PRINT *, 'end_minute = ', end_minute
!PRINT *, 'end_second = ', end_second

!PRINT *,'init date: ' , init_date,'  init time: ', init_time
!PRINT *,'end date: ', end_date, '  end time: ', end_time

length = end_second - init_second                      !seconds
length = length + 60*(end_minute-init_minute)          !minutes
length = length + 3600*(end_hour - init_hour)          !hours
length = length + 86400*(end_day - init_day)           !days
length = length + 2635200*(end_month - init_month)     !months
length = length + 31356000*(end_year - init_year)      !years

days = length/86400
hours = length/3600 - 24*days
minutes = length/60 - 60*hours - 60*24*days
seconds = length - 60*minutes - 3600*hours - 3600*24*days

IF(PRESENT(fileID)) THEN
if(rank.eq.0) WRITE(fileID,*) '----------Calculated Time-------------'
if(rank.eq.0) WRITE(fileID,*) ''
if(rank.eq.0) WRITE(fileID,*) '   *Total Seconds:  ', length
if(rank.eq.0) WRITE(fileID,*) '   *Total Hours:    ', length/3600
if(rank.eq.0) WRITE(fileID,*) '   *Readable Time:  ', days , ' Days, ', hours,  &
' Hours, ', minutes, ' Minutes, ', seconds, ' seconds'
if(rank.eq.0) WRITE(fileID,*) '______________________________________'

ELSE
if(rank.eq.0) PRINT *, '----------Calculated Time-------------'
if(rank.eq.0) Print *, ''
if(rank.eq.0) PRINT *, '   *Total Seconds:  ', length
if(rank.eq.0) PRINT *, '   *Total Hours:    ', length/3600
if(rank.eq.0) PRINT *, '   *Readable Time:  ', days , ' Days, ', hours,  &
 ' Hours, ', minutes, ' Minutes, ', seconds, ' seconds'
if(rank.eq.0) PRINT *, '______________________________________'
END IF

END SUBROUTINE EndTiming

SUBROUTINE PrintTiming()
IMPLICIT NONE
#ifdef MPI
include "mpif.h"
#endif
character(len=8)::date
character(len=10)::time
character(len=5)::zone
integer::values(8)
integer::ierr,rank,Nsize

#ifdef MPI
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,Nsize,ierr)
#else
 rank=0
 Nsize=1
#endif

call DATE_AND_TIME(date,time,zone,values)
if(rank.eq.0) write(*,*) "date=",date
if(rank.eq.0) write(*,*) "time=",values(5),":",values(6),":",values(7)
if(rank.eq.0) write(*,*) "--------------------------------------------------------"
if(rank.eq.0) write(*,*) " "
if(rank.eq.0) write(*,*) " "
if(rank.eq.0) write(*,*) " "
END SUBROUTINE PrintTiming
END MODULE timing_module
 
