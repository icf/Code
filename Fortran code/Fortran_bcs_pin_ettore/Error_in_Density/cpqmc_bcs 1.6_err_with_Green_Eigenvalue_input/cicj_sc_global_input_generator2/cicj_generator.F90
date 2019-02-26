!main program of code
program main

implicit none
integer:: Nsite,Nsitex,Nsitey
complex(kind=8)::cicj_sc_global(2*16,2*16)
complex(kind=8)::rho(16)
real(kind=8)::distance

integer::sitex,sitey
integer::i,j

Nsitex=4
Nsitey=4
Nsite=16

!U=4 4433
rho(0+1)= ( 0.18750000020155699     ,  0.0000000000000000     )
rho(1+1)= ( 0.12271853821724019     ,  1.1755088683879263E-011)
rho(2+1)= (  6.0334959582000129E-002,  5.8094461726641541E-012)
rho(3+1)= ( 0.12271853806358329     ,  9.4730242866907625E-012)
rho(4+1)= ( 0.12271853808508064     ,  1.3979911446247975E-011)
rho(5+1)= (  6.0334959587875832E-002,  8.5846840708145110E-012)
rho(6+1)= ( -1.2602079723448621E-003,  6.8864408232737778E-012)
rho(7+1)= (  6.0334959469512388E-002,  2.0773465204281888E-012)
rho(8+1)= (  6.0334959498272396E-002, -1.3592766201367117E-011)
rho(9+1)= ( -1.2602079331338061E-003, -1.2285458554820805E-011)
rho(10+1)= ( -6.2628351858797851E-002, -4.0466772555104405E-011)
rho(11+1)= ( -1.2602080902464317E-003, -4.1973399915229984E-011)
rho(12+1)= ( 0.12271853812043160     ,  2.0544047247073587E-012)
rho(13+1)= (  6.0334959654799493E-002, -9.0684603103231261E-013)
rho(14+1)= ( -1.2602079401649211E-003, -2.8174011062035492E-011)
rho(15+1)= (  6.0334959473238269E-002, -4.1557135719734882E-012)

!U=8 77  
rho(0+1)= 4.3750000613530587e-01
rho(1+1)= 1.4085383576049862e-01
rho(2+1)= -9.4824216729012829e-03
rho(3+1)=1.2765463297796770e-01
rho(4+1)=1.2995041768833560e-01 
rho(5+1)= 2.2484782411706507e-02
rho(6+1)=  -3.3746564741084652e-02
rho(7+1)=2.2856439734942939e-02
rho(8+1)=-1.2345418579995084e-02
rho(9+1)=  -3.5464855046015466e-02
rho(10+1)= -2.7101453135364592e-02
rho(11+1)=-2.8694844186302960e-02
rho(12+1)=1.3100473191349937e-01 
rho(13+1)= 2.8258905172913379e-02  
rho(14+1)= -3.4287339607711131e-02
rho(15+1)=1.7881540445067623e-02

!U=8 77  twist(0.01,0.02)
rho(0+1)=( 4.3749998476193230e-01  ,  0.0000000000000000e+00)
rho(1+1)=(  1.3231482118476523e-01  ,  4.4732388586367442e-03)
rho(2+1)=( -4.4473144481448248e-02  , -1.0820942057112801e-08)
rho(3+1)=(1.3231411706440971e-01  , -4.4732240600597285e-03)
rho(4+1)=(1.3242516278577826e-01  ,  9.3275584253363058e-03)
rho(5+1)=( 3.9799937186233147e-02  ,  5.6438824381261916e-04)
rho(6+1)=( -3.2993451762258173e-02 ,  -8.5515506192273237e-03)
rho(7+1)=( 3.9417312367152728e-02  ,  2.0050499927566742e-04)
rho(8+1)=(-4.4157178193562664e-02  , -2.0630720005618275e-10)
rho(9+1)=( -3.3157370164776802e-02 ,  -4.0612429121947976e-03)
rho(10+1)=( -2.7182278236518906e-02  , -2.6640607481533476e-09)
rho(11+1)=(-3.3157016285184859e-02  ,  4.0612491879768765e-03)
rho(12+1)=(1.3242617803379175e-01  , -9.3275899894553690e-03)
rho(13+1)=( 3.9417962548860042e-02 ,  -2.0056721052966555e-04 ) 
rho(14+1)=(  -3.2993970730726213e-02 ,   8.5515443892756997e-03)
rho(15+1)=(3.9800059213273994e-02  , -5.6444509146797513e-04)

    


cicj_sc_global=0
do i=1,Nsitex
do j=1,Nsitey  !the standard point
   do sitex=1,Nsitex
      do sitey=1,Nsitey  

         if(sitex-i .GE. 0 .AND. sitey-j .GE. 0)cicj_sc_global(i+(j-1)*Nsitex,sitex+(sitey-1)*Nsitex) &
&=rho(sitex-i+(sitey-j)*Nsitex +1) 
         if(sitex-i .LT. 0 .AND. sitey-j .GE. 0)cicj_sc_global(i+(j-1)*Nsitex,sitex+(sitey-1)*Nsitex) &
&=rho(sitex-i+4+(sitey-j)*Nsitex +1) 
         if(sitex-i .GE. 0 .AND. sitey-j .LT. 0)cicj_sc_global(i+(j-1)*Nsitex,sitex+(sitey-1)*Nsitex) &
&=rho(sitex-i+(sitey-j+4)*Nsitex +1) 
         if(sitex-i .LT. 0 .AND. sitey-j .LT. 0)cicj_sc_global(i+(j-1)*Nsitex,sitex+(sitey-1)*Nsitex) &
&=rho(sitex-i+4+(sitey-j+4)*Nsitex +1) 
      enddo
   enddo
enddo
enddo

!spin symmetry
do i=1,Nsite
do j=1,Nsite
   cicj_sc_global(i+Nsite,j+Nsite)=cicj_sc_global(i,j)
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

!real number generator: generate "2*Nsite by 2*Nsite" input matrix
open(unit=10,file='real_number_cicj_sc_global_input.inputdat',status='UNKNOWN')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        write(10,*)dble(cicj_sc_global(i,j))
  end do
end do
close(10)

!imag number generator: generate "2*Nsite by 2*Nsite" input matrix
open(unit=10,file='imag_number_cicj_sc_global_input.inputdat',status='UNKNOWN')
do i=1,2*Nsite,1
  do j=1,2*Nsite,1
        write(10,*)aimag(cicj_sc_global(i,j))
  end do
end do
close(10)
 
end program main
