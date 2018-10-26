!------------------------------
!normalization the wavefunction
!------------------------------
subroutine norm_wave(wav)
use param
implicit none
complex(kind=8),intent(inout)::wav(Nhilbert)
complex(kind=8)::norm,sqrtnorm,temp
real(kind=8),external::dznrm2

norm=dcmplx(dznrm2(Nhilbert,wav(1),1))
sqrtnorm=1.d0/norm
call zscal(Nhilbert,sqrtnorm,wav(1),1)


!write(*,*) dznrm2(Nhilbert,wav(1),1);pause
end subroutine norm_wave

!---------------------
!S^2 wav - S(S+1) wave
!---------------------
subroutine S_con(wav)
use param
implicit none
complex(kind=8),intent(inout)::wav(Nhilbert)
complex(kind=8)::tmp
complex(kind=8),allocatable::tmpwav(:)
integer(kind=8)::i,j
integer::iS
real(kind=8)::SS
real(kind=8)::renorm,arpha

allocate(tmpwav(Nhilbert))
do iS=two_SS,0,-2
   if(iS.eq.two_rSS) cycle

   SS=dble(iS)/2.d0
   !if(rank.eq.0) write(*,*) SS

   call Square_towf(wav(1),tmpwav(1))

   renorm=dble(two_rSS)/2.d0*(dble(two_rSS)/2.d0+1.d0)
   renorm=renorm-SS*(SS+1.d0)
   arpha=-SS*(SS+1.d0)
   call zaxpy(Nhilbert,dcmplx(arpha),wav,1,tmpwav,1)
   arpha=1.d0/renorm
   call zscal(Nhilbert,dcmplx(arpha),tmpwav,1)
   call zcopy(Nhilbert,tmpwav,1,wav,1) 
end do

deallocate(tmpwav)
end subroutine S_con

!----------------------
!Ti wav -exp(i*ki) wave
!----------------------
subroutine k_con(wav)
use param
implicit none
complex(kind=8),intent(inout)::wav(Nhilbert)
complex(kind=8)::tmp
complex(kind=8),allocatable::tmpwav(:)
integer(kind=8)::i,j
integer::dd,ki
complex(kind=8)::renorm,arpha


allocate(tmpwav(Nhilbert))
do dd=1,Dimen,1
   do ki=0,Nl(dd)-1,1
      !if(rank.eq.0) write(*,*) ki,k_keep(dd)
      if(ki.eq.k_keep(dd)) cycle

      call k_towf(wav,tmpwav,dd)

      renorm=exp(-Xi*2.d0*Pi*k_keep(dd)/dble(Nl(dd)))
      renorm=renorm-exp(-Xi*2.d0*Pi*ki/dble(Nl(dd)))
      arpha=-exp(-Xi*2.d0*Pi*ki/dble(Nl(dd)))
      call zaxpy(Nhilbert,arpha,wav,1,tmpwav,1)
      arpha=1.d0/renorm
      call zscal(Nhilbert,arpha,tmpwav,1)
      call zcopy(Nhilbert,tmpwav,1,wav,1)      
   end do
end do
deallocate(tmpwav)
end subroutine k_con

