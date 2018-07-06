
      SUBROUTINE ASSERT(COND,MSG)
      implicit none

      LOGICAL COND

      CHARACTER(LEN=*) MSG

      IF (.NOT.COND) THEN

         WRITE(0,*) 'ASSERT:',MSG

         STOP

      END IF

      END SUBROUTINE ASSERT

