!This subroutine should compile with the module caldet_module and matrixcal subroutine
!This page contains some matrix manipulate in cpqmc paticularly.
!Contains: 
!1.wavefunction overlap:m*n and n*m wave function overlap
!2.Green function matrix:The  <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> matrix
!3.Fast Green Function matrix:when we get inverse of the overlap, we can use it
!  straightly,do no need to calculate it again.
!4.Green function element:Only calculate i,j element,sometimes,we do no need the
!  whole matrix, only need some element.
!5.Fast Green Function element.Only calculate i,j element. We use the known
!  inverse overlap matrix information.
!6.Green function element,only calculate  <phi_l|ci^+ ci|phi_r>, it is a
!derivation method, which might have small error due to numerical
!derivation. It can be used when <phi_l|phi_r>=0
!7.Green function element,only calculate  <phi_l|ci^+ ci|phi_r>, it is a
!multi determinates method, which computational time scale to m slatter determinates.
!It can be used when <phi_l|phi_r>=0
!8.Green function element,only calculate <phi_l|ci^+ ci|phi_r>, it is a exact
!expand method, which calculate <phi_l|exp(ni)|phi_r>, can get the green
!function element by transfer, it can be used when <phi_l|phi_r>=0 

!This subroutine calculate <phi_l|phi_r>,Tanspose[Nsite*N(sigma)].Nsite*N(sigma)
subroutine deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: ovlpinv_temp(m,m)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)

 call zgemm('C','N',m,m,n,one,phi_l,n,phi_r,n,zero,ovlpinv_temp,m)
 !ovlpinv_temp= matmul(transpose(conjg(phi_l)), phi_r)
 !write(*,*) ovlpinv_temp
end subroutine



!This subroutine calculate Det[<phi_l|phi_r>]
subroutine deter_overlap_imp(n,m,phi_l,phi_r,imp)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT)::imp
complex(kind=8):: ovlpinv_temp(m,m)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)

 call zgemm('C','N',m,m,n,one,phi_l,n,phi_r,n,zero,ovlpinv_temp,m)
 call caldet(m,ovlpinv_temp(1:m,1:m),imp)
end subroutine deter_overlap_imp




!This subroutine calcuate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r>,all the i,j matrix
!is wrote into Amax(n,n)
subroutine cal_Amat(n,m,phi_l,phi_r,Amat)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat(n,n)

complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)


complex(kind=8):: tempmat(n,m)
complex(kind=8):: ovlpinv_temp(m,m)
real(kind=8)::rdummy

  !ovlpinv_temp= matmul(transpose(conjg(phi_l)), phi_r)
  call deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)

  call  inverse(ovlpinv_temp,m)

  call ZGEMM('N','N',n,m,m,one,phi_r,n,ovlpinv_temp,m,zero,tempmat,n)
  call ZGEMM('N','c',n,n,m,one,tempmat,n,phi_l,n,zero,Amat,n)

  Amat=transpose(Amat)

end subroutine cal_Amat


!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r>,here we use
!the information of ovlpinv, which is the inverse of <phi_l|phi_r>
subroutine cal_Amat_withovlpinv(n,m,phi_l,phi_r,ovlpinv,Amat)
implicit none

integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m), ovlpinv(m,m)
complex(kind=8),intent(OUT):: Amat(n,n)

complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)


integer:: i,j
integer:: k,l

complex(kind=8):: tempmat(n,m)


   call ZGEMM('N','N',n,m,m,one,phi_r,n,ovlpinv,m,zero,tempmat,n)
   call ZGEMM('N','c',n,n,m,one,tempmat,n,phi_l,n,zero,Amat,n)

   Amat=transpose(Amat)

  ! G=RO^(-1)L^{dagger}
end subroutine cal_Amat_withovlpinv



!This subroutine also calculate <phi_l|ci cj^+|phi_r>/<phi_l|phi_r>,here we use
!the information of ovlpinv, which is the inverse of <phi_l|phi_r>
subroutine cal_Amat_withovlpinv2(n,m,phi_l,phi_r,ovlpinv,Amat)
implicit none

integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m), ovlpinv(m,m)
complex(kind=8),intent(OUT):: Amat(n,n)

complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)


integer:: i,j
integer:: k,l

complex(kind=8):: tempmat(n,m)


   call ZGEMM('N','N',n,m,m,one,phi_r,n,ovlpinv,m,zero,tempmat,n)
   call ZGEMM('N','c',n,n,m,one,tempmat,n,phi_l,n,zero,Amat,n)

   Amat=-Amat
   do i=1,n
     Amat(i,i)=Amat(i,i)+one
   enddo

  ! G=RO^(-1)L^{dagger}
end subroutine cal_Amat_withovlpinv2


!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> element.
subroutine cal_cidcj(n,m,phi_l,phi_r,Amat,i,j)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat
integer,intent(IN)::i,j
complex(kind=8):: ovlpinv_temp(m,m)
integer::k,l
 !Get the inverse of phi_l over lap phi_r
  call deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)
  call inverse(ovlpinv_temp,m)

 !Get the element.
  Amat=dcmplx(0.d0,0.d0)
  do k=1,m,1
     do l=1,m,1
        Amat=Amat+phi_r(j,k)*ovlpinv_temp(k,l)*conjg(phi_l(i,l))
     end do
  end do
end subroutine cal_cidcj


!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> element.
!We use the ovlpinv to faster the simulation
subroutine cal_cidcj_withovlpinv(n,m,phi_l,phi_r,ovlpinv,Amat,i,j)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m),ovlpinv(m,m)
complex(kind=8),intent(OUT):: Amat
integer,intent(IN)::i,j
integer::k,l

 !Get the element.
  Amat=dcmplx(0.d0,0.d0)
  do k=1,m,1
     do l=1,m,1
        Amat=Amat+phi_r(j,k)*ovlpinv(k,l)*conjg(phi_l(i,l))
     end do
  end do
end subroutine cal_cidcj_withovlpinv


!Green function element,only calculate  <phi_l|ci^+ ci|phi_r>, it is a
!derivation method, which might have small error due to numerical
!derivation. It can be used when <phi_l|phi_r>=0
subroutine cal_ni_derivation(n,m,phi_l,phi_r,det,i,Amat)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
character(len=1),intent(IN)::det
integer,intent(IN)::i
complex(kind=8),intent(OUT):: Amat
complex(kind=8)::phi_rn(n,m),ovp,ovpeta,ovlp(m,m) !(ovpeta-ovp)/(eta-zero)
real(kind=8)::eta
integer::j,k

!set the derivation parameter
eta=1d-5

!Get the overlap of phi_l.phi_r
if(det.eq.'Y') then !zero overlap
  ovp=dcmplx(0.d0,0.d0)
else if(det.eq.'N') then !None zero overlap
  call deter_overlap(n,m,phi_l,phi_r,ovlp)
  call caldet(m,ovlp,ovp)
else
write(*,*) "Something is wrong with det input in cal_ni_derivation"
stop
end if

!Get exp(eta.ni).phi_r
phi_rn(1:n,1:m)=phi_r(1:n,1:m)
do j=1,m,1
   phi_rn(i,j)=phi_r(i,j)*exp(eta)
end do

!Get det(phi_l.exp(eta.ni).phi_r)
call deter_overlap(n,m,phi_l,phi_rn,ovlp)
call caldet(m,ovlp,ovpeta)

!Get the <phi_l|ci^+ ci|phi_r>
Amat=(ovpeta-ovp)/eta
end subroutine cal_ni_derivation

!Green function element,only calculate  <phi_l|ci^+ ci|phi_r>, it is a
!multi determinates method, which computational time scale to m slatter
!determinates. It can be used when <phi_l|phi_r>=0
subroutine cal_ni_multi(n,m,phi_l,phi_r,i,Amat)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
integer,intent(IN)::i
complex(kind=8),intent(OUT):: Amat
complex(kind=8)::ovlp(m,m),ovp
complex(kind=8),allocatable::phi_rn(:,:,:)
integer::j,k

allocate(phi_rn(n,m,m))

!set the multi determinate
do j=1,m,1
   phi_rn(1:n,1:m,j)=phi_r(1:n,1:m)
   do k=1,n,1
      if(k.ne.i) then
        phi_rn(k,j,j)=dcmplx(0.d0,0.d0)
      end if
   end do
end do

!Get the <phi_l|ci^+ ci|phi_r>
Amat=dcmplx(0.d0,0.d0)
do j=1,m,1
   call deter_overlap(n,m,phi_l(1:n,1:m),phi_rn(1:n,1:m,j),ovlp(1:m,1:m))
   call caldet(m,ovlp(1:m,1:m),ovp)
   Amat=Amat+ovp
end do

deallocate(phi_rn)
end subroutine cal_ni_multi


!Green function element,only calculate <phi_l|ci^+ ci|phi_r>, it is a exact
!expand method, which calculate <phi_l|exp(ni)|phi_r>, can get the green
!function element by transfer, it can be used when <phi_l|phi_r>=0
subroutine cal_ni_exp(n,m,phi_l,phi_r,det,i,Amat)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
character(len=1),intent(IN)::det
integer,intent(IN)::i
complex(kind=8),intent(OUT):: Amat
complex(kind=8)::phi_rn(n,m),ovp,ovpexp,ovlp(m,m) !ovp: <phi_l|phi_r> 
                                                  !ovpexp:<phi_l|exp(ni)|phi_r>
integer::j,k

!Get the overlap of phi_l.phi_r
if(det.eq.'Y') then !zero overlap
  ovp=dcmplx(0.d0,0.d0)
else if(det.eq.'N') then !None zero overlap
  call deter_overlap(n,m,phi_l,phi_r,ovlp)
  call caldet(m,ovlp,ovp)
else
write(*,*) "Something is wrong with det input in cal_ni_derivation"
stop
end if

!Get exp(eta.ni).phi_r
phi_rn(1:n,1:m)=phi_r(1:n,1:m)
do j=1,m,1
   phi_rn(i,j)=phi_r(i,j)*exp(1.d0)
end do

!Get det(phi_l.exp(ni).phi_r)
call deter_overlap(n,m,phi_l,phi_rn,ovlp)
call caldet(m,ovlp,ovpexp)

!Get the <phi_l|ci^+ ci|phi_r>
Amat=dcmplx(1.d0/(exp(1.d0)-1.d0))*(ovpexp-ovp)
end subroutine cal_ni_exp



!Green function element, calculate <phi_l|ci^+ cj|phi_r>, it is a exact
!expand method, which calculate <phi_l|exp(ci^+cj)|phi_r>, can get the green
!function element by transfer, it can be used when <phi_l|phi_r>=0
subroutine cal_ij_exp(n,m,phi_l,phi_r,det,i,j,Amat)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
character(len=1),intent(IN)::det
integer,intent(IN)::i,j
complex(kind=8),intent(OUT):: Amat
complex(kind=8)::phi_rn(n,m),ovp,ovpexp,ovlp(m,m) !ovp: <phi_l|phi_r> 
                                                  !ovpexp:<phi_l|exp(ni)|phi_r>
integer::k,p,q

!Get the overlap of phi_l.phi_r
if(det.eq.'Y') then !zero overlap
  ovp=dcmplx(0.d0,0.d0)
else if(det.eq.'N') then !None zero overlap
  call deter_overlap(n,m,phi_l,phi_r,ovlp)
  call caldet(m,ovlp,ovp)
else
write(*,*) "Something is wrong with det input in cal_ni_derivation"
stop
end if


if(i.eq.j) then
  !Get exp(ni).phi_r
  phi_rn(1:n,1:m)=phi_r(1:n,1:m)
  do p=1,m,1
     phi_rn(i,p)=phi_r(i,p)*exp(1.d0)
  end do
  !Get det(phi_l.exp(ni).phi_r)
  call deter_overlap(n,m,phi_l,phi_rn,ovlp)
  call caldet(m,ovlp,ovpexp)
  
  !Get the <phi_l|ci^+ ci|phi_r>
  Amat=dcmplx(1.d0/(exp(1.d0)-1.d0))*(ovpexp-ovp)
else
  !Get exp(ci^+cj).phi_r
  phi_rn(1:n,1:m)=phi_r(1:n,1:m)
  do p=1,m,1
     phi_rn(i,p)=phi_rn(i,p)+phi_r(j,p)
  end do
  !Get det(phi_l.exp(ci^+cj).phi_r)
  call deter_overlap(n,m,phi_l,phi_rn,ovlp)
  call caldet(m,ovlp,ovpexp)

  !Get the <phi_l|ci^+ cj|phi_r>
  Amat=ovpexp-ovp
end if
end subroutine cal_ij_exp
