!    Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009  J. Williams, I. Danshita, R. Mishmash, D. Schirmer, M. L. Wall
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
MODULE io_module
!
! Purpose: Module to perform basic i/o operations
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
IMPLICIT NONE

INTERFACE appendBaseName
MODULE PROCEDURE appendBaseName_r, appendBaseName_i, appendBaseName_c
END INTERFACE appendBaseName

CONTAINS

SUBROUTINE createFileName(basename,diRectory)
!
!Purpose: Begin a file name in the directory diRectory
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

baseName=diRectory

END SUBROUTINE createFileName

SUBROUTINE appendBaseName_r(basename,partName,partDigs,partValue)
!
!Purpose: Append to a file name the character string partName followed by 
!the value partValue to partDigs digits for real numbers
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname
INTEGER, INTENT(IN) :: partDigs
REAL(KIND=8), INTENT(IN) :: partValue
CHARACTER(16) :: iWstring, specString

!Write the number of decimal places wanted into a string
	WRITE(iWstring,'(I4)') partDigs
!Generate the format type for a float string with partDigs decimal places	
	specString="(1f16."//TRIM(ADJUSTL(iWstring))//")"
!Write the value partValue to a string using partDigs decimal places
	WRITE(iWstring,specString) partValue
!append partname and the value to the given basename
basename=TRIM(basename)//partName//TRIM(ADJUSTL(iWstring))

END SUBROUTINE appendBaseName_r

SUBROUTINE appendBaseName_i(basename,partName,partValue)
!
!Purpose: Append to a file name the character string partName followed by 
!the value partValue to partDigs digits for integers
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname
INTEGER, INTENT(IN) :: partValue
CHARACTER(16) :: iWstring, specString

!Write the number of digits wanted into a string
	WRITE(iWstring,'(I4)') 16
!Generate the format type for an integer string with partDigs digits	
	specString="(I"//TRIM(ADJUSTL(iWstring))//")"
!Write the value partValue to a string using partDigs digits
	WRITE(iWstring,specString) partValue
!append partname and the value to the given basename
basename=TRIM(basename)//partName//TRIM(ADJUSTL(iWstring))

END SUBROUTINE appendBaseName_i

SUBROUTINE appendBaseName_c(basename,partName)
!
!Purpose: Append to a file name the character string partName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname

basename=TRIM(basename)//TRIM(partName)

END SUBROUTINE appendBaseName_c

SUBROUTINE copyName(name1,name2)
!
!Purpose:Copy name1 to name2
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: name1
CHARACTER(len=*), INTENT(OUT) :: name2

name2=TRIM(name1)

END SUBROUTINE copyName

LOGICAL FUNCTION CheckName(baseName)
!
!Purpose: Returns TRUE if file name baseName exists and FALSE otherwise
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: baseName

INQUIRE(FILE=baseName, EXIST=CheckName)

END FUNCTION CheckName

SUBROUTINE openUnit(fileName,myUnit,openKind)
!
!Purpose: Open a file 'filename' as UNIT=myUnit
!An error message is returned if the specified file cannot be opened
!
IMPLICIT NONE
INTEGER, INTENT(IN) ::  myUnit
CHARACTER(len=*), INTENT(IN) :: fileName
CHARACTER(132) ::  stopname
CHARACTER, INTENT(IN), OPTIONAL :: openKind
INTEGER :: fileStatus !Status integer to ensure file open happens properly

stopname='*** Cannot open file named '//filename//'***'

IF(PRESENT(openKind)) THEN

IF(openKind=='N') THEN
OPEN(UNIT=myUnit, FILE=filename, STATUS='NEW', ACTION='WRITE',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF

ELSE IF(openKind=='R') THEN
OPEN(UNIT=myUnit, FILE=filename, STATUS='REPLACE',ACTION='READWRITE',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF


ELSE IF(openKind=='O') THEN
OPEN(UNIT=myUnit, FILE=filename, STATUS='OLD', ACTION='READWRITE',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF

ELSE IF(openKind=='A') THEN
OPEN(UNIT=myUnit, FILE=filename, ACTION='WRITE',POSITION='APPEND',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF

ELSE IF(openKind=='B') Then 

OPEN(UNIT=myUnit, FILE=filename, STATUS='OLD',Form='Unformatted',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF


ELSE IF(openKind=='C') Then

OPEN(UNIT=myUnit, FILE=filename,STATUS='REPLACE',Form='Unformatted',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF


else 
STOP "Unknown option in openUnit!"
END IF

ELSE
OPEN(UNIT=myUnit, FILE=filename, STATUS='UNKNOWN', ACTION='READWRITE',IOSTAT=fileStatus)
IF(fileStatus>0) THEN
PRINT *,stopname
STOP
END IF
END IF

END SUBROUTINE openUnit

END MODULE io_module
