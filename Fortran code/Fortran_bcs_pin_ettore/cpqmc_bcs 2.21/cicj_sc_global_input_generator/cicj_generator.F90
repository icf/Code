!main program of code
program main

implicit none
integer:: Nsite,Nsitex,Nsitey
complex(kind=8)::cicj_sc_global(2*16,2*16)
complex(kind=8)::rho01,rho02,rho1,rho2,rho3,rho4,rho5
real(kind=8)::distance

integer::sitex,sitey
integer::i,j

Nsitex=4
Nsitey=4
Nsite=16

!U=4 77  
rho01=4.3750134113590200e-01 !may have some problem 
rho02=4.3750134113590200e-01 !may have some problem
rho1= 1.6716243704626116e-01    9.5618305263811248e-03
rho2= 5.4219127339802040e-02    2.3747915469839599e-04
rho3=-5.0582689894488900e-02   -9.2844283438057174e-03
rho4=-5.0467122247422232e-02   -2.2228959533793100e-02
rho5=-0.051

!U=4 55
rho01=0.3125 !may have some problem 
rho02=0.3125 !may have some problem
rho1=0.1759
rho2=0.05716
rho3=0.05716
rho4=-0.05616
rho5=-0.1634

!U=4 33 u=8
rho01=0.1875 !may have some problem 
rho02=0.1875 !may have some problem
rho1=0.1196
rho2=0.0573
rho3=0.0573
rho4=-0.003
rho5=-0.0627

!U=8
!rho01=0.4375 !may have some problem 
!rho02=0.4375 !may have some problem
!rho1=0.1323
!rho2=0.0398
!rho3=-0.0444
!rho4=-0.033
!rho5=-0.0272

cicj_sc_global=0

do sitex=1,Nsitex
do sitey=1,Nsitey
   if( mod((sitex+sitey),2) .EQ. 1 )cicj_sc_global(sitex+(sitey-1)*Nsitex,sitex+(sitey-1)*Nsitex)=rho02
   if( mod((sitex+sitey),2) .EQ. 0 )cicj_sc_global(sitex+(sitey-1)*Nsitex,sitex+(sitey-1)*Nsitex)=rho01
   if( mod((sitex+sitey),2) .EQ. 1 )cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+sitex+(sitey-1)*Nsitex)=rho01
   if( mod((sitex+sitey),2) .EQ. 0 )cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+sitex+(sitey-1)*Nsitex)=rho02
   do i=1,Nsitex
      do j=1,Nsitey
         distance=(sitex-i)**2+(sitey-j)**2
         !spin up
         if(distance .EQ. 1)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho1   
         if(distance .EQ. 9)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho1  
         if(distance .EQ. 2)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho2 
         if(distance .EQ. 10)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho2
         if(distance .EQ. 18)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho2 
         if(distance .EQ. 4)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho3 
         if(distance .EQ. 5)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho4
         if(distance .EQ. 13)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho4 
         if(distance .EQ. 8)cicj_sc_global(sitex+(sitey-1)*Nsitex,i+(j-1)*Nsitex)=rho5
         !spin dn
         if(distance .EQ. 1)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho1   
         if(distance .EQ. 9)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho1  
         if(distance .EQ. 2)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho2 
         if(distance .EQ. 10)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho2
         if(distance .EQ. 18)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho2 
         if(distance .EQ. 4)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho3 
         if(distance .EQ. 5)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho4
         if(distance .EQ. 13)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho4 
         if(distance .EQ. 8)cicj_sc_global(Nsite+sitex+(sitey-1)*Nsitex,Nsite+i+(j-1)*Nsitex)=rho5
      enddo
   enddo
enddo
enddo

!generator: generate "2*Nsite by 2*Nsite" input matrix
open(unit=10,file='cicj_sc_global_input.inputdat',status='UNKNOWN')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        write(10,*)cicj_sc_global(i,j)
  end do
end do
close(10)
 
end program main
