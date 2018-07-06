        program datdinf
        implicit none

        integer, parameter:: Maxstep=4000
        integer, parameter:: Datset=8

        real(kind=8):: datin(Datset), ddat(Datset)
        real(kind=8):: time(Maxstep), nbin(Maxstep)
        real(kind=8):: dat(Maxstep,Datset), datq(Maxstep,Datset)
       
        integer:: il,i
        character(len=132):: filename, str

        integer::  len, status, ios

        call get_command_argument(1, str, len, status)
        read(str,'(a132)') filename


        print *, 'filename', filename

          dat  = 0d0
          datq = 0d0
          nbin = 0d0

        open (11,file=trim(filename))
         do while(ios/=-1) 
              read(11,*,iostat=ios) il,time(il),datin(1:Datset)
              do i=1,Datset
              dat(il,i)  = dat(il,i)+datin(i)
              datq(il,i) = datq(il,i)+datin(i)**2
              enddo 
              nbin(il) = nbin(il)+1d0
                   if (il>Maxstep) then 
                      print *, "il exceed Maxstep", il
                      stop
                    endif 
         enddo
         close(11)



         OPEN (UNIT=12,FILE=trim(filename)//'.stat')

         do il=1,Maxstep
          do i=1,Datset
       ddat(i)=dsqrt((datq(il,i)/nbin(il)-(dat(il,i)/nbin(il))**2)/(nbin(il)-1.d0))
          enddo 
         write (12,1000) il,time(il),(dat(il,i)/nbin(il) ,ddat(i),i=1,Datset),nbin(il)
         enddo

        1000  format(i6,40f17.6)

       stop

       end program 
       
