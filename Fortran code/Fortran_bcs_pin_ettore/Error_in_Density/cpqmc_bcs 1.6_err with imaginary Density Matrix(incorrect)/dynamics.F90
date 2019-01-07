subroutine measure_green_holes(Gd,cicj_t_local)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none

complex(kind=8),intent(IN) ::Gd(2*Nsite,2*Nsite,Nbeta,Dtot)
complex(kind=8),intent(OUT)::cicj_t_local(2*Nsite,Nbeta,Dtot)

integer         ::i_beta,k,j,j2,j3,i,m,n
integer,external::latt_label
integer,external::bound
integer         ::cc(1:Dimen),ctmp
complex(kind=8) ::GG(2*Nsite,2*Nsite,Nbeta)

integer ::i_deb,j_deb !DEB

cicj_t_local=zero

do k=1,Dtot
  do i_beta=1,Nbeta
    if(i_beta.eq.1)then

      call copy_G_dc(Gd(1,1,i_beta,k),GG(1,1,i_beta))

    else


      if(dtype.EQ.'c') then
        call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,GG(1,1,i_beta-1),2*Nsite              &
              &   ,Gd(1,1,i_beta,k)                                                          &
              &   ,2*Nsite,zero,GG(1,1,i_beta),2*Nsite)
      else if(dtype.EQ.'d') then
        call zgemm('N','N',Nsite,Nsite,Nsite,one,GG(1,1,i_beta-1),2*Nsite                    &
              &   ,Gd(1,1,i_beta,k)                                                          &
              &   ,2*Nsite,zero,GG(1,1,i_beta),2*Nsite)
        call zgemm('N','N',Nsite,Nsite,Nsite,one,GG(Nsite+1,Nsite+1,i_beta-1),2*Nsite        &
              &   ,Gd(Nsite+1,Nsite+1,i_beta,k)                                              &
                  ,2* Nsite,zero,GG(Nsite+1,Nsite+1,i_beta),2*Nsite)
      endif

    endif


    do i=1,Nsite,1
      do j=1,Nsite,1

        do m=1,Dimen,1
          ctmp=coor(i,m)-coor(j,m)+1
          cc(m)=bound(ctmp,Nl(m))
        enddo

        n=latt_label(cc(1:Dimen))
        cicj_t_local(n,i_beta,k)=cicj_t_local(n,i_beta,k)+GG(i,j,i_beta)
        cicj_t_local(n+Nsite,i_beta,k)=cicj_t_local(n+Nsite,i_beta,k)+GG(i+Nsite,j+Nsite,i_beta)

      enddo
    enddo

    do n=1,2*Nsite,1
      cicj_t_local(n,i_beta,k)=cicj_t_local(n,i_beta,k)/dcmplx(Nsite)
    enddo

  enddo
enddo


end subroutine measure_green_holes










subroutine measure_green_particles(Gd,cicj_t_local)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none

complex(kind=8),intent(IN) ::Gd(2*Nsite,2*Nsite,Nbeta,Dtot)
complex(kind=8),intent(OUT)::cicj_t_local(2*Nsite,Nbeta,Dtot)

integer         ::i_beta,k,j,j2,j3,i,m,n
integer,external::latt_label
integer,external::bound
integer         ::cc(1:Dimen),ctmp
complex(kind=8) ::GG(2*Nsite,2*Nsite,Nbeta)

integer ::i_deb,j_deb !DEB

cicj_t_local=zero

do k=1,Dtot
  do i_beta=1,Nbeta
    if(i_beta.eq.1)then

      call copy_G_dc(Gd(1,1,i_beta,k),GG(1,1,i_beta))


!      if(dtype.EQ.'c') then
!        do j=1,2*Nsite,1
!          do j2=1,2*Nsite,1
!            GG(j2,j,i_beta)            =Gd(j2,j,i_beta,k)
!          enddo
!        enddo
!      elseif(dtype.EQ.'d') then
!        do j=1,Nsite,1
!          do j2=1,Nsite,1
!            GG(j2,j,i_beta,k)            =Gd(j2,j,i_beta,k)
!            GG(j2+Nsite,j+Nsite,i_beta,k)=Gd(j2+Nsite,j+Nsite,i_beta,k)
!          enddo
!        enddo
!      endif

    else


      if(dtype.EQ.'c') then
        call zgemm('N','N',2*Nsite,2*Nsite,2*Nsite,one,Gd(1,1,i_beta,k),2*Nsite              &
              &   ,GG(1,1,i_beta-1)                                                          &
              &   ,2*Nsite,zero,GG(1,1,i_beta),2*Nsite)
      else if(dtype.EQ.'d') then
        call zgemm('N','N',Nsite,Nsite,Nsite,one,Gd(1,1,i_beta,k),2*Nsite                  &
              &   ,GG(1,1,i_beta-1)                                                         &
              &   ,2*Nsite,zero,GG(1,1,i_beta),2*Nsite)
        call zgemm('N','N',Nsite,Nsite,Nsite,one,Gd(Nsite+1,Nsite+1,i_beta,k),2*Nsite       &
              &   ,GG(Nsite+1,Nsite+1,i_beta-1)                                             & 
                  ,2* Nsite,zero,GG(Nsite+1,Nsite+1,i_beta),2*Nsite)
      endif

!      if(dtype.EQ.'c') then
!        do j=1,2*Nsite,1
!          do j2=1,2*Nsite,1
!            GG(j2,j,i_beta,k)=zero
!            do j3=1,2*Nsite,1
!              GG(j2,j,i_beta,k)=GG(j2,j,i_beta,k)+Gd(j2,j3,i_beta,k)*GG(j3,j,i_beta-1,k)
!            enddo
!          enddo
!        enddo
!      elseif(dtype.EQ.'d') then
!        do j=1,Nsite,1
!          do j2=1,Nsite,1
!
!            GG(j2,j,i_beta,k)=zero
!            do j3=1,Nsite,1
!              GG(j2,j,i_beta,k)=GG(j2,j,i_beta,k)+Gd(j2,j3,i_beta,k)*GG(j3,j,i_beta-1,k)
!            enddo
!
!            GG(j2+Nsite,j+Nsite,i_beta,k)=zero
!            do j3=1,Nsite,1
!              GG(j2+Nsite,j+Nsite,i_beta,k)=GG(j2+Nsite,j+Nsite,i_beta,k)      &
!                           +Gd(j2+Nsite,j3+Nsite,i_beta,k)*GG(j3+Nsite,j+Nsite,i_beta-1,k)
!            enddo
!          enddo
!        enddo
!      endif

    endif



    do i=1,Nsite,1
      do j=1,Nsite,1
    
        do m=1,Dimen,1
          ctmp=coor(i,m)-coor(j,m)+1
          cc(m)=bound(ctmp,Nl(m))
        enddo 

        n=latt_label(cc(1:Dimen))
        cicj_t_local(n,i_beta,k)=cicj_t_local(n,i_beta,k)+GG(i,j,i_beta)
        cicj_t_local(n+Nsite,i_beta,k)=cicj_t_local(n+Nsite,i_beta,k)+GG(i+Nsite,j+Nsite,i_beta)
     
      enddo
    enddo

    do n=1,2*Nsite,1
      cicj_t_local(n,i_beta,k)=cicj_t_local(n,i_beta,k)/dcmplx(Nsite)
    enddo

  enddo
enddo


end subroutine measure_green_particles









subroutine calculate_Gp(ph1,ph2,G)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none

complex(kind=8),intent(IN) ::ph1(2*Nsite,Ntot),ph2(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::G(2*Nsite,2*Nsite)

integer        ::alpha,beta
complex(kind=8)::imp_local
complex(kind=8)::ovp_local(Ntot,Ntot)
!complex(kind=8)::Amat_local(2*Nsite,2*Nsite)

integer        :: i_deb,j_deb  !DEBUG

!write(*,*)'******** IN calculate_Gp '
!write(*,*)
!write(*,*)'ph1  up   '
!write(*,*)
!do i_deb=1,Nsite
!  write(*,*)(ph1(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo
!write(*,*)

!write(*,*)
!write(*,*)'Phi1  dn    '
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!   write(*,*)(ph1(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
!enddo
!write(*,*)

!write(*,*)'****'

!write(*,*)
!write(*,*)'ph2  up   '
!write(*,*)
!do i_deb=1,Nsite
!  write(*,*)(ph2(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo
!write(*,*)

!write(*,*)
!write(*,*)'Phi2  dn    '
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!   write(*,*)(ph2(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
!enddo
!write(*,*)



call over_lap_dc(ph1,ph2,ovp_local)
!ovp_local = ph1^dagger ph2


!write(*,*)'**** ovp_local'
!do i_deb=1,Nspin(1)
!  write(*,*)(ovp_local(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo
!write(*,*)
!do i_deb=Nspin(1)+1,Ntot
!  write(*,*)(ovp_local(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
!enddo
!write(*,*)


call inverse_d_dc(ovp_local,imp_local)
!ovp_local = Inverse of ph1^dagger ph2, imp_local = det(ph1^dagger ph2)

!write(*,*)'**** inverse of ovp_local'
!do i_deb=1,Nspin(1)
!  write(*,*)(ovp_local(i_deb,j_deb),j_deb=1,Nspin(1))
!enddo
!write(*,*)
!do i_deb=Nspin(1)+1,Ntot
!  write(*,*)(ovp_local(i_deb,j_deb),j_deb=Nspin(1)+1,Ntot)
!enddo
!write(*,*)
!write(*,*)
!write(*,*)'imp_local = ',imp_local



call cal_Amat_withovlpinv_dc2(ph1,ph2,ovp_local,G)


!write(*,*)
!write(*,*)
!write(*,*)'Amat_local'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!     write(*,*)i_deb,j_deb,Amat_local(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!     write(*,*)i_deb,j_deb,Amat_local(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)

!if(dtype.EQ.'c') then
!  do alpha=1,2*Nsite
!    do beta=1,2*Nsite
!      G(alpha,beta)=-Amat_local(beta,alpha)
!    enddo
!    G(alpha,alpha)=G(alpha,alpha)+one
!  enddo
!else if(dtype.EQ.'d') then
!  do alpha=1,Nsite
!    do beta=1,Nsite
!      G(alpha,beta)=-Amat_local(beta,alpha)
!      G(alpha+Nsite,beta+Nsite)=-Amat_local(beta+Nsite,alpha+Nsite)
!    enddo
!    G(alpha,alpha)=G(alpha,alpha)+one
!    G(alpha+Nsite,alpha+Nsite)=G(alpha+Nsite,alpha+Nsite)+one
!  enddo
!endif

end subroutine calculate_Gp

!----------------------------------------------------------------------

subroutine calculate_Bh(i_beta,i_walker,B)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none

integer, intent(IN)         ::i_beta,i_walker
complex(kind=8), intent(OUT)::B(2*Nsite,2*Nsite)

integer        ::j,j2,j3,icall
integer        :: i_deb,j_deb  !DEBUG
integer        ::aux
real(kind=8)   ::x
complex(kind=8)::explr_up,explr_dn
complex(kind=8)::Temp(2*Nsite,2*Nsite)
!data icall/0/
!save 

!DEBUG
!write(*,*)
!write(*,*)
!write(*,*)'exp_halfK'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)

if(i_beta.eq.Nbeta)then
B=zero
Temp=zero
!icall=1
endif

call copy_B_dc(exp_mhalfK,B)

!write(*,*)
!write(*,*)
!write(*,*)'B new initial'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)

!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    B(j2,j)=exp_halfK(j2,j)
!    B(j2+Nsite,j+Nsite)=exp_halfK(j2+Nsite,j+Nsite)
!  enddo
!enddo


!write(*,*)
!write(*,*)
!write(*,*)'B old initial'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)

do j=1,Nsite,1

  x=back_store(j,i_beta,i_walker)
  aux=NINT(x)
 ! explr_up=expln_up(aux)
 ! explr_dn=expln_dn(aux)
  if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

  elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

  endif

  do j2=1,Nsite,1
    B(j2,j)=B(j2,j)                        /(explr_up+one)
    B(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)/(explr_dn+one)
  enddo

enddo


call copy_B_dc(B,Temp)

!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    Temp(j2,j)=B(j2,j)
!    Temp(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)
!  enddo
!enddo


call zgemm('N','N',Nsite,Nsite,Nsite,one,Temp(1,1),2*Nsite,exp_mhalfK(1,1),2*Nsite,zero,B,2*Nsite)
call zgemm('N','N',Nsite,Nsite,Nsite,one,Temp(Nsite+1,Nsite+1),2*Nsite,exp_mhalfK(Nsite+1,Nsite+1)&
     &    ,2*Nsite,zero,B(Nsite+1,Nsite+1),2*Nsite)


!write(*,*)
!write(*,*)
!write(*,*)'B new final'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)


!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    B(j2,j)=dcmplx(0.d0,0.d0)
!    do j3=1,Nsite,1
!      B(j2,j)=B(j2,j)+Temp(j2,j3)*exp_halfK(j3,j)
!    enddo
!    B(j2+Nsite,j+Nsite)=dcmplx(0.d0,0.d0)
!    do j3=1,Nsite,1
!      B(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)+Temp(j2+Nsite,j3+Nsite)*exp_halfK(j3+Nsite,j+Nsite)
!    enddo
!  enddo
!enddo

!write(*,*)
!write(*,*)
!write(*,*)'B old final'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)
!
!stop !DEB


!write(*,*)
!write(*,*)
!write(*,*)'exp_K'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,exp_K(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,exp_K(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!stop 'DEB'


end subroutine calculate_Bh


!----------------------------------------------------------------------

subroutine calculate_B(i_beta,i_walker,B)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
implicit none

integer, intent(IN)         ::i_beta,i_walker
complex(kind=8), intent(OUT)::B(2*Nsite,2*Nsite)

integer        ::j,j2,j3,icall
integer        :: i_deb,j_deb  !DEBUG
integer        ::aux
real(kind=8)   ::x
complex(kind=8)::explr_up,explr_dn
complex(kind=8)::Temp(2*Nsite,2*Nsite)
!data icall/0/
!save 

!DEBUG
!write(*,*)
!write(*,*)
!write(*,*)'exp_halfK'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,exp_halfK(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)

if(i_beta.eq.Nbeta)then
B=zero
Temp=zero
!icall=1
endif

call copy_B_dc(exp_halfK,B)

!write(*,*)
!write(*,*)
!write(*,*)'B new initial'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)

!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    B(j2,j)=exp_halfK(j2,j)
!    B(j2+Nsite,j+Nsite)=exp_halfK(j2+Nsite,j+Nsite)
!  enddo
!enddo


!write(*,*)
!write(*,*)
!write(*,*)'B old initial'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)




do j=1,Nsite,1
  
  x=back_store(j,i_beta,i_walker)
  aux=NINT(x)
!  explr_up=expln_up(aux)
!  explr_dn=expln_dn(aux)

  if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

  elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

  endif

  do j2=1,Nsite,1
    B(j2,j)=B(j2,j)                        *(explr_up+one)
    B(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)*(explr_dn+one)
  enddo

enddo


call copy_B_dc(B,Temp)

!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    Temp(j2,j)=B(j2,j)
!    Temp(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)
!  enddo
!enddo


call zgemm('N','N',Nsite,Nsite,Nsite,one,Temp(1,1),2*Nsite,exp_halfK(1,1),2*Nsite,zero,B,2*Nsite)
call zgemm('N','N',Nsite,Nsite,Nsite,one,Temp(Nsite+1,Nsite+1),2*Nsite,exp_halfK(Nsite+1,Nsite+1)&
     &    ,2*Nsite,zero,B(Nsite+1,Nsite+1),2*Nsite)


!write(*,*)
!write(*,*)
!write(*,*)'B new final'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)


!do j=1,Nsite,1
!  do j2=1,Nsite,1
!    B(j2,j)=dcmplx(0.d0,0.d0)
!    do j3=1,Nsite,1
!      B(j2,j)=B(j2,j)+Temp(j2,j3)*exp_halfK(j3,j)
!    enddo
!    B(j2+Nsite,j+Nsite)=dcmplx(0.d0,0.d0)
!    do j3=1,Nsite,1
!      B(j2+Nsite,j+Nsite)=B(j2+Nsite,j+Nsite)+Temp(j2+Nsite,j3+Nsite)*exp_halfK(j3+Nsite,j+Nsite)
!    enddo
!  enddo
!enddo

!write(*,*)
!write(*,*)
!write(*,*)'B old final'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,B(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!write(*,*)
!
!stop !DEB


!write(*,*)
!write(*,*)
!write(*,*)'exp_K'
!do i_deb=1,Nsite
!  do j_deb=1,Nsite
!    write(*,*)i_deb,j_deb,exp_K(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!do i_deb=Nsite+1,2*Nsite
!  do j_deb=Nsite+1,2*Nsite
!    write(*,*)i_deb,j_deb,exp_K(i_deb,j_deb)
!  enddo
!enddo
!write(*,*)
!stop 'DEB'
  
    
end subroutine calculate_B



!----------------------------------------------------------------------

subroutine one_step_left(phR,i_beta,i_walker)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
use mc_loop_param
implicit none
integer,intent(IN)::i_beta,i_walker
complex(kind=8),intent(INOUT)::phR(2*Nsite,Ntot)

integer        ::j,k
integer        ::aux
real(kind=8)   ::x
real(kind=8)   ::anm,anm1,anm2
complex(kind=8)::explr_up,explr_dn
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::phtmp(2*Nsite,Ntot)

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

do j=1,Nsite,1

  x=back_store(j,i_beta,i_walker)
  aux=NINT(x)
!  explr_up=expln_up(aux)
!  explr_dn=expln_dn(aux)
  if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

  elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

  endif
  
  if(dtype.EQ.'c') then
    do k=1,Ntot,1
      phR(j,k)      =phR(j,k)      *(explr_up+one)
      phR(j+Nsite,k)=phR(j+Nsite,k)*(explr_dn+one)
    end do
  else if(dtype.EQ.'d') then
    do k=1,Nspin(1),1
      phR(j,k)=phR(j,k)*(explr_up+one)
    end do
    do k=Nspin(1)+1,Ntot,1
      phR(j,k)=phR(j,k)*(explr_dn+one)
    end do
  end if

enddo

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

if(mod(i_beta,StepforGram).eq.0)then
if(dtype.EQ.'c') then
  call modGS(phR(1,1),2*Nsite,Ntot,anm,Rmat)
else if(dtype.EQ.'d') then
  call modGS(phR(1:Nsite,1:Nspin(1)),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
  call modGS(phR((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))
  anm=anm1*anm2
end if
endif

end subroutine one_step_left


!----------------------------------------------------------------------

subroutine one_step_left_new(phR,i_beta,anm,i_walker)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
use mc_loop_param
implicit none
integer,intent(IN)::i_beta,i_walker
complex(kind=8),intent(INOUT)::phR(2*Nsite,Ntot)
real(kind=8),intent(OUT)::anm

integer        ::j,k,i_b
integer        ::aux
real(kind=8)   ::x
real(kind=8)   ::anm1,anm2
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::phtmp(2*Nsite,Ntot)

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

i_b=i_beta

do j=1,Nsite,1

  x=back_store(j,i_beta,i_walker)

   if(i_b.LE.Nbk(1).AND.i_b.GE.1) then
     if(kcrn.eq.1) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.2) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.3) then
       lr_up=x*sqrt(dcmplx(dt*onsitU))
       lr_up=lr_up-(dt*onsitU/2.d0)

       lr_dn=-x*sqrt(dcmplx(dt*onsitU))
       lr_dn=lr_dn-(dt*onsitU/2.d0)

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.4) then
       lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else
       write(*,*) "Something is wrong with kcrn input:",kcrn
       call mystop
     end if
   else if(i_b.GT.Nbk(1).AND.i_b.LE.Nstps_fwd) then
     aux=NINT(x)
!     explr_up=expln_up(aux)
!     explr_dn=expln_dn(aux)
     if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

     elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

     endif
   else
     write(*,*) "Something is wrong of i_b",i_b
     call mystop
   end if



  if(dtype.EQ.'c') then
    do k=1,Ntot,1
      phR(j,k)      =phR(j,k)      *(explr_up+one)
      phR(j+Nsite,k)=phR(j+Nsite,k)*(explr_dn+one)
    end do
  else if(dtype.EQ.'d') then
    do k=1,Nspin(1),1
      phR(j,k)=phR(j,k)*(explr_up+one)
    end do
    do k=Nspin(1)+1,Ntot,1
      phR(j,k)=phR(j,k)*(explr_dn+one)
    end do
  end if

enddo

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

anm=1.d0
if(mod(i_beta,StepforGram).eq.0)then
if(dtype.EQ.'c') then
  call modGS(phR(1,1),2*Nsite,Ntot,anm,Rmat)
else if(dtype.EQ.'d') then
  call modGS(phR(1:Nsite,1:Nspin(1)),Nsite,Nspin(1),anm1,Rmat(1:Nspin(1),1:Nspin(1)))
  call modGS(phR((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot),Nsite,Nspin(2),anm2,Rmat(1:Nspin(2),1:Nspin(2)))
  anm=anm1*anm2
end if
endif

end subroutine one_step_left_new

!----------------------------------------------------------------------



!----------------------------------------------------------------------

subroutine one_step_left_bare(phR,i_beta,i_walker)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
use mc_loop_param
implicit none
integer,intent(IN)::i_beta,i_walker
complex(kind=8),intent(INOUT)::phR(2*Nsite,Ntot)

integer        ::j,k,i_b
integer        ::aux
real(kind=8)   ::x
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
complex(kind=8)::phtmp(2*Nsite,Ntot)

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

do j=1,Nsite,1
  x=back_store(j,i_beta,i_walker)
  i_b=i_beta
  if(i_b.LE.Nbk(1).AND.i_b.GE.1) then
     if(kcrn.eq.1) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.2) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.3) then
       lr_up=x*sqrt(dcmplx(dt*onsitU))
       lr_up=lr_up-(dt*onsitU/2.d0)

       lr_dn=-x*sqrt(dcmplx(dt*onsitU))
       lr_dn=lr_dn-(dt*onsitU/2.d0)

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.4) then
       lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else
       write(*,*) "Something is wrong with kcrn input:",kcrn
       call mystop
     end if
  else if(i_b.GT.Nbk(1).AND.i_b.LE.Nstps_fwd) then 
     aux=NINT(x)
!     explr_up=expln_up(aux)
!     explr_dn=expln_dn(aux)
     if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

     elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

     endif 
  else
     write(*,*) "Something is wrong of i_b one_step_left_bare",i_b
     call mystop
  end if
   

  if(dtype.EQ.'c') then
    do k=1,Ntot,1
      phR(j,k)      =phR(j,k)      *(explr_up+one)
      phR(j+Nsite,k)=phR(j+Nsite,k)*(explr_dn+one)
    end do
  else if(dtype.EQ.'d') then
    do k=1,Nspin(1),1
      phR(j,k)=phR(j,k)*(explr_up+one)
    end do
    do k=Nspin(1)+1,Ntot,1
      phR(j,k)=phR(j,k)*(explr_dn+one)
    end do
  end if

enddo

call copy_wf_dc(phR,phtmp)
call k_to_ph_dc(exp_halfK,phtmp,phR)

end subroutine one_step_left_bare

!----------------------------------------------------------------------


subroutine one_step_right(phL,i_beta,i_walker,coe)

use param
use lattice_param
use model_param
use method_param
use phiT_param
use phi_param
use project_param
use mc_loop_param
implicit none
integer,intent(IN)::i_beta,i_walker
complex(kind=8),intent(INOUT)::coe(Dtot)
complex(kind=8),intent(INOUT)::phL(2*Nsite,Ntot,Dtot)

integer        ::j,k,l,i_b
integer        ::aux
real(kind=8)   ::x
real(kind=8)   ::anm,anm1,anm2
complex(kind=8)::explr_up,explr_dn,lr_up,lr_dn
complex(kind=8)::Rmat(Ntot,Ntot)
complex(kind=8)::phtmp(2*Nsite,Ntot)

do k=1,Dtot,1
   call copy_wf_dc(phL(1,1,k),phtmp(1,1))
   call dk_to_ph_dc(exp_halfK,phtmp(1,1),phL(1,1,k))  !dagger
end do

i_b=i_beta

do j=1,Nsite,1

  x=back_store(j,i_beta,i_walker)

   if(i_b.LE.Nbk(1).AND.i_b.GE.1) then
     if(kcrn.eq.1) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.2) then
       if(x.LT.0.5d0) then
         lr_up=-1.d0*(dt*onsitU*0.5d0+gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0+gamaf)
       else
         lr_up=-1.d0*(dt*onsitU*0.5d0-gamaf)
         lr_dn=-1.d0*(dt*onsitU*0.5d0-gamaf)
       end if
       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.3) then
       lr_up=x*sqrt(dcmplx(dt*onsitU))
       lr_up=lr_up-(dt*onsitU/2.d0)

       lr_dn=-x*sqrt(dcmplx(dt*onsitU))
       lr_dn=lr_dn-(dt*onsitU/2.d0)

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else if(kcrn.eq.4) then
       lr_up=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_up=lr_up+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       lr_dn=x*sqrt(dcmplx(-1.d0*dt*onsitU))
       lr_dn=lr_dn+(dt*onsitU/2.d0)*(1.d0-2.d0*ng(j))

       explr_up=exp(lr_up)-one
       explr_dn=exp(lr_dn)-one
     else
       write(*,*) "Something is wrong with kcrn input:",kcrn
       call mystop
     end if
   else if(i_b.GT.Nbk(1).AND.i_b.LE.Nstps_fwd) then
     aux=NINT(x)
!     explr_up=expln_up(aux)
!     explr_dn=expln_dn(aux)

     if(Nbands.eq.1)then

       explr_up=expln_up(aux)
       explr_dn=expln_dn(aux)

     elseif(Nbands.eq.3)then

       if(j.le.Nbravais)then
         explr_up=explnd_up(aux)
         explr_dn=explnd_dn(aux)
       elseif(j.le.2*Nbravais)then
         explr_up=explnx_up(aux)
         explr_dn=explnx_dn(aux)
       else
         explr_up=explny_up(aux)
         explr_dn=explny_dn(aux)
       endif

     endif

   else
     write(*,*) "Something is wrong of i_b one_step_right",i_b
     call mystop
   end if

  

  if(dtype.EQ.'c') then
     do k=1,Dtot,1
       do l=1,Ntot,1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   else if(dtype.EQ.'d') then
     do k=1,Dtot,1
       do l=1,Nspin(1),1
          phL(j,l,k)=phL(j,l,k)*conjg(explr_up+one)
       end do
       do l=Nspin(1)+1,Ntot,1
          phL(j+Nsite,l,k)=phL(j+Nsite,l,k)*conjg(explr_dn+one)
       end do
     end do
   end if

enddo

do k=1,Dtot,1
   call copy_wf_dc(phL(1,1,k),phtmp(1,1))
   call dk_to_ph_dc(exp_halfK,phtmp(1,1),phL(1,1,k))  !dagger
end do

if(mod(i_beta,StepforGram).eq.0)call gs_phL(phL,coe)

end subroutine one_step_right
