!This page is used to do population control

subroutine PopControl()
use param
use adET_param
use mc_loop_param
use phi_param
use mpi_serial_param
implicit none
logical::reconfig
integer::i
#ifdef MPI
 if(rank.eq.0) allocate(Mweight(MNwalkers))
#endif


!We set m_w in assess_walkers,if reconfig m_w=mean reweight
!else m_w=1.d0, we do not do any population control,no reweight.
m_w_old=m_w

!DEBUG
!write(*,*)
!write(*,*)'In popcontrol '
!write(*,*)'m_w_old, m_w', m_w_old, m_w
!write(*,*)

call assess_walkers(reconfig)

if(reconfig) then
  call population_control()
end if

#ifdef MPI
 if(rank.eq.0) deallocate(Mweight)
#endif



!DEBUG
!write(*,*)
!write(*,*)'Exit from  popcontrol '
!write(*,*)'m_w_old, m_w', m_w_old, m_w
!write(*,*)
end subroutine PopControl


!Make the population control while the population size is constant.
subroutine population_control()
use param
use phi_param
use phiT_param
use mc_loop_param
use mpi_serial_param
use method_param
use model_param
use lattice_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer:: table(Nwalkers)
integer:: i,j,k,l,m,n
#ifdef MPI
integer,allocatable::Mtable(:)
integer,allocatable::Nsend(:) !The number send out in each tread, and the directions.
integer,allocatable::Nscatter(:) !It is the Nsend(i)*2, for the scatterv
integer::Nsedth !The total number need to be send through different thread
integer,allocatable::msdtable(:,:) !The total send out table: contains parents and children
integer,allocatable::sdtable(:,:) !Send out table in each thread:contains parents and children
integer,allocatable::disp(:) ! For scatterv msdtable to sdtable
integer::nbl,nrk  !Nlabel and Nrank
integer,allocatable::SREQ(:),SSTA(:,:)
integer,allocatable::RREQ(:),RSTA(:,:)
integer::is,ir,Ns,Nr
character,allocatable::sbuf(:,:),rbuf(:,:)
integer::Nbuf
integer::posit
integer::is_tmp,ir_tmp  !is_tmp is the send number i tmp
                        !ir_tmp is the receive rank i_rank tmp
                        !ir_tmp also represent receive number i tmp
#endif

real(kind=8)::cc
real(kind=8)::wt(Nwalkers)


!If cpmc then only use weight
!if fp,rcp we use the ccoe
if(crn.GT.0) then
  cc=ccoe
else
  cc=0.d0
end if



#ifdef MPI


if(rank.eq.0) allocate(Mtable(MNwalkers))
if(rank.eq.0) call reconfiguration(MNwalkers,Mweight,Mtable)


!Scatter Mtable,Nreceive
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_SCATTER(Mtable(1),Nwalkers,MPI_INTEGER,table(1),Nwalkers,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


!Nsend: The numbers need to be sent from different thread
allocate(Nsend(0:(Nsize-1)))
if(rank.eq.0) then
  Nsend=0
  do m=0,Nsize-1,1
     do n=1,Nwalkers,1
        i=n+m*Nwalkers
        j=Mtable(i)
        call nthd(j,Nwalkers,k)
        if(k.ne.m) then 
          Nsend(k)=Nsend(k)+1
        end if
     end do
  end do

  Nsedth=0
  do i=0,Nsize-1,1
     Nsedth=Nsedth+Nsend(i)
  end do
  if(Nsedth.GT.0) then
    allocate(msdtable(2,Nsedth))
  else
    allocate(msdtable(1,1))
  end if

  l=1
  do m=0,Nsize-1,1
     do n=1,Nwalkers,1
        i=n+m*Nwalkers
        j=Mtable(i)
        call nthd(j,Nwalkers,k)
        if(k.ne.m) then
          msdtable(1,l)=Mtable(i) !parent
          msdtable(2,l)=i   !child
          l=l+1
        end if
     end do
  end do
  if((l-1).NE.Nsedth) then
    write(*,*) "Something is wrong with Nsend:",Nsedth,l-1
    call mystop
  end if

end if
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_BCAST(Nsend(0),Nsize,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)



!Scatter the sendtable
!disp: The first number to be sent in different thread.
!sdtable(rank): the ones need to be sent from rank
allocate(disp(0:(Nsize-1)))
allocate(Nscatter(0:(Nsize-1)))
do i=0,Nsize-1,1
   Nscatter(i)=Nsend(i)*2
end do
disp(0)=0
do j=1,Nsize-1,1
   disp(j)=disp(j-1)+Nscatter(j-1)
end do
if(Nsend(rank).GT.0) then 
  allocate(sdtable(2,Nsend(rank)))
else if(Nsend(rank).EQ.0) then
  allocate(sdtable(1,1))
else
  write(*,*) "Nsend rank is lt 0",rank,Nsend(rank)
  call mystop
end if
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_SCATTERV(msdtable(1,1),Nscatter(0),disp(0),MPI_INTEGER,sdtable(1,1),Nscatter(rank),MPI_INTEGER,0,MPI_COMM_WORLD,IERR)



!The population control--message expression
Nbuf=0
if(dtype.EQ.'c') then
  Nbuf=Nbuf+2*Nsite*Ntot*16+16
else if(dtype.EQ.'d') then
  Nbuf=Nbuf+Nsite*Ntot*16+16
end if
!if(crn.LT.0.d0) then
  Nbuf=Nbuf+(Dtot+1)*16
!else if(crn.GT.0.d0) then
!end if
if(back_pro) then
  Nbuf=Nbuf+Nsite*i_back*8
  if(dtype.EQ.'c') then
    Nbuf=Nbuf+2*Nsite*Ntot*16
  else if(dtype.EQ.'d') then
    Nbuf=Nbuf+Nsite*Ntot*16
  end if
end if


!Packed and Send message.
is=0;is_tmp=0;ir_tmp=-1
do i=1,Nsend(rank),1
   !Receive array
   l=sdtable(2,i)
   nbl=l
   call nthd(nbl,Nwalkers,nrk)
    

   !if(l.ne.sdtable(1,i)) then
     !if(rank.ne.nrk) then
       if(sdtable(1,i).EQ.is_tmp.AND.nrk.eq.ir_tmp) then
       else
         is=is+1
         is_tmp=sdtable(1,i)
         ir_tmp=nrk
       end if
     !end if
   !end if
end do
Ns=is
if(Ns.GT.0) then
  allocate(sbuf(Nbuf,Ns))
  allocate(SREQ(Ns),SSTA(MPI_STATUS_SIZE,Ns))
end if
is=0;is_tmp=0;ir_tmp=-1
do i=1,Nsend(rank),1
   !Send array
   k=sdtable(1,i)
   call nthd(k,Nwalkers,nrk)
   if(rank.NE.nrk) then
     write(*,*) "Something is wrong with population send.",rank
     call mystop
   end if
   !Receive array
   l=sdtable(2,i)
   nbl=l
   call nthd(nbl,Nwalkers,nrk)


   !if(rank.ne.nrk) then
     if(sdtable(1,i).EQ.is_tmp.AND.nrk.eq.ir_tmp) then
     else
       is=is+1
       is_tmp=sdtable(1,i)
       ir_tmp=nrk

       call pack_pop(Nbuf,k,sbuf(1,is))

       call MPI_ISEND(sbuf(1,is),Nbuf,MPI_BYTE,nrk,nbl,MPI_COMM_WORLD,SREQ(is),IERR)
     end if
   !end if
end do


!Recieve message and change weight dlogw
ir=0;is_tmp=0
do i=1,Nwalkers,1
   l=i+rank*Nwalkers !The number labeled in total MPI
   j=table(i) !j is the parent,the number labeled in total MPI
   nbl=j
   call nthd(nbl,Nwalkers,nrk) !nbl labeled in each thread, nrk means rank

   !if (l<j) then
      if(nrk.ne.rank) then
        if(j.eq.is_tmp) then
        else
          ir=ir+1
          is_tmp=j
        end if
      end if
   !end if 

enddo
Nr=ir
if(Nr.GT.0) then
  allocate(rbuf(Nbuf,Nr))
  allocate(RREQ(Nr),RSTA(MPI_STATUS_SIZE,Nr))
end if
ir=0;is_tmp=0;ir_tmp=0
do i=1,Nwalkers,1
   l=i+rank*Nwalkers !The number labeled in total MPI
   j=table(i) !j is the parent,the number labeled in total MPI
   nbl=j
   call nthd(nbl,Nwalkers,nrk) !nbl labeled in each thread, nrk means rank 


   if (l.ne.j) then
      if(nrk.eq.rank) then
        call change_i_to_j(i,nbl)
      else
        if(j.eq.is_tmp) then
          call change_i_to_j(i,ir_tmp)
        else
          ir=ir+1
          is_tmp=j
          ir_tmp=i
          !receive a pack
          call MPI_RECV(rbuf(1,ir),Nbuf,MPI_BYTE,nrk,i,MPI_COMM_WORLD,RSTA(1,ir),IERR)
          call unpack_pop(Nbuf,i,rbuf(1,ir))
        end if!j compare with is_tmp
      end if!nrk compare with rank
   end if

   weight(i)=1.d0/(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc)
   dlogw(i)=dlog(1.d0/(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc))
enddo


if(Ns.GT.0) call MPI_WAITALL(Ns,SREQ,SSTA,IERR)
if(allocated(sbuf)) deallocate(sbuf)
if(allocated(rbuf)) deallocate(rbuf)
if(allocated(SREQ)) deallocate(SREQ)
if(allocated(SSTA)) deallocate(SSTA)
if(allocated(RREQ)) deallocate(RREQ)
if(allocated(RSTA)) deallocate(RSTA)




deallocate(disp)
deallocate(sdtable)
deallocate(Nsend)
deallocate(Nscatter)
if(rank.eq.0) deallocate(Mtable)
if(rank.eq.0) deallocate(msdtable)
call MPI_BARRIER(MPI_COMM_WORLD,IERR)


#else

do i=1,Nwalkers,1
   wt(i)=weight(i)*(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc)
end do

call reconfiguration(Nwalkers,wt,table)

do i=1,Nwalkers,1
   j=table(i) !j is the parent  

   if (i.NE.j) then
      call change_i_to_j(i,j)
   end if

   weight(i)=1.d0/(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc)
   dlogw(i)=dlog(1.d0/(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc))
enddo

#endif
end subroutine population_control




!In the population control,change walker i from walker j
subroutine change_i_to_j(i,j)
use model_param
use phiT_param
use phi_param
use lattice_param
use method_param
implicit none
integer,intent(IN)::i,j

if(dtype.EQ.'c') then
  phi(:,:,i)=phi(:,:,j)
  rx(i)=rx(j)

  !if(crn.LT.0.d0) then
    impfunc(:,i)=impfunc(:,j)
    tot_imp(i)=tot_imp(j)
  !else if(crn.GT.0.d0) then
  !end if


  if(back_pro) then
     back_store(1:Nsite,1:i_back,i)=back_store(1:Nsite,1:i_back,j)
     phi0(:,:,i)=phi0(:,:,j)
  end if
else if(dtype.EQ.'d') then
  phi(1:Nsite,1:Nspin(1),i)=phi(1:Nsite,1:Nspin(1),j)
  phi((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,i)=phi((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,j)
  rx(i)=rx(j)

  !if(crn.LT.0.d0) then
    impfunc(:,i)=impfunc(:,j)
    tot_imp(i)=tot_imp(j)
  !else if(crn.GT.0.d0) then
  !end if


  if(back_pro) then
     back_store(1:Nsite,1:i_back,i)=back_store(1:Nsite,1:i_back,j)
     phi0(1:Nsite,1:Nspin(1),i)=phi0(1:Nsite,1:Nspin(1),j)
     phi0((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,i)=phi0((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,j)
  end if
end if

end subroutine change_i_to_j



!pack the information to sbuf in population control
subroutine pack_pop(Nbuf,k,sbuf)
use phi_param
use phiT_param
use mpi_serial_param
use method_param
use lattice_param
implicit none
#ifdef MPI
include "mpif.h"
integer,intent(IN)::Nbuf
integer,intent(IN)::k
character,intent(OUT)::sbuf(Nbuf)
integer::posit

posit=0
!call MPI_PACK(phi(1,1,k),2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
call MPI_PACK(phi(1,1,k),1,phtype,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
call MPI_PACK(rx(k),1,MPI_DOUBLE_COMPLEX,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
!if(crn.LT.0.d0) then
  call MPI_PACK(impfunc(1,k),Dtot,MPI_DOUBLE_COMPLEX,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
  call MPI_PACK(tot_imp(k),1,MPI_DOUBLE_COMPLEX,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
!else if(crn.GT.0.d0) then
!end if

if(back_pro) then
  call MPI_PACK(back_store(1,1,k),Nsite*i_back,MPI_DOUBLE_PRECISION,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
  !call MPI_PACK(phi0(1,1,k),2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
  call MPI_PACK(phi0(1,1,k),1,phtype,sbuf(1),Nbuf,posit,MPI_COMM_WORLD,IERR)
end if

if(posit.NE.Nbuf) then
  write(*,*) "Something is wrong with packed",rank,posit,Nbuf
  call mystop
end if
#else
integer,intent(IN)::Nbuf
integer,intent(IN)::k
character,intent(IN)::sbuf(Nbuf)
return
#endif
end subroutine pack_pop


!unpack the information from rbuf in population control
subroutine unpack_pop(Nbuf,i,rbuf)
use phi_param
use phiT_param
use mpi_serial_param
use method_param
use lattice_param
implicit none
integer,intent(IN)::Nbuf
integer,intent(IN)::i
character,intent(IN)::rbuf(Nbuf)
integer::posit
#ifdef MPI
include "mpif.h"

posit=0
!call MPI_UNPACK(rbuf(1),Nbuf,posit,phi(1,1,i),2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERR)
call MPI_UNPACK(rbuf(1),Nbuf,posit,phi(1,1,i),1,phtype,MPI_COMM_WORLD,IERR)
call MPI_UNPACK(rbuf(1),Nbuf,posit,rx(i),1,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERR)
!if(crn.LT.0.d0) then
  call MPI_UNPACK(rbuf(1),Nbuf,posit,impfunc(1,i),Dtot,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERR)
  call MPI_UNPACK(rbuf(1),Nbuf,posit,tot_imp(i),1,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERR)
!else if(crn.GT.0.d0) then
!end if 

if(back_pro) then
  call MPI_UNPACK(rbuf(1),Nbuf,posit,back_store(1,1,i),Nsite*i_back,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
  !call MPI_UNPACK(rbuf(1),Nbuf,posit,phi0(1,1,i),2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERR)
  call MPI_UNPACK(rbuf(1),Nbuf,posit,phi0(1,1,i),1,phtype,MPI_COMM_WORLD,IERR)
end if

if(posit.NE.Nbuf) then
  write(*,*) "Something is wrong with packed",rank,posit,Nbuf
  call mystop
end if
#else
return
#endif
end subroutine unpack_pop


!change to a new way, try to move less walkers
subroutine reconfiguration(n,w,table)
use rand_num
use mc_loop_param
implicit none
integer,intent(IN):: n
real(kind=8),intent(IN):: w(n)
integer,intent(OUT):: table(n)
real(kind=8),allocatable::p(:)
integer,allocatable::num_n(:)
real(kind=8):: wsum,z,cap
integer::n_s
integer:: i,j,k,l,m


allocate(p(n),num_n(n))

wsum=0d0
do j=1,n
   wsum=wsum+w(j)
enddo

!---Prepare the probability
do j=1,n
   p(j)=w(j)/wsum
enddo


!wsum=0.d0
!cap=1.d0  !0.1d0
!do j=1,n,1
!   if(p(j).GT.cap) then
!     write(*,*) "cap the p(j)",j,p(j),cap
!     p(j)=cap
!   end if
!   wsum=wsum+p(j)
!enddo
!do j=1,n
!   p(j)=p(j)/wsum
!enddo



!---initialize-----------
j=1
wsum=p(1)
num_n=0
do i=1,n
   z=dble((i-1)+rndm())/dble(n)
   do while(z>wsum)
      j=j+1
      wsum=wsum+p(j)
   enddo
   !table(i)=j
   num_n(j)=num_n(j)+1
   table(j)=j !if num_n>0, it is itself.
enddo


!test
!do i=1,n,1
!   table(i)=i
!end do
!num_n=0
!num_n(2)=22
!num_n(20)=2
!write(*,*) "----------------"
!write(*,*) num_n


if(mod(n,Nwalkers).NE.0) then
  write(*,*) "Warning!The total population should be devided by Nwalkers:",n,Nwalkers
  call mystop
end if
n_s=n/Nwalkers


!We first move the l in rank(i-1) to k in rank(i-1)
do i=1,n_s,1
   m=(i-1)*Nwalkers+1 
   do j=1,Nwalkers,1
      k=j+(i-1)*Nwalkers
      if(num_n(k).eq.0) then
        do l=m,i*Nwalkers,1
           if(num_n(l).GT.1) then
              num_n(l)=num_n(l)-1
              num_n(k)=num_n(k)+1
              table(k)=l
              exit
           end if
        end do
        m=l
      end if
   end do  
end do

!The population change in different core
if(n_s.GT.1) then
  m=1
  do k=1,n,1
     if(num_n(k).eq.0) then
       do l=m,n,1
          if(num_n(l).GT.1) then
             num_n(l)=num_n(l)-1
             num_n(k)=num_n(k)+1
             table(k)=l
             exit
          end if
       end do  
       m=l
     end if
  end do
end if

!Check
do k=1,n,1
   if(num_n(k).NE.1) then
     write(*,*) "Something is wrong in population control, num_n(k) is:",num_n(k),k
     call mystop
   end if
end do

deallocate(p,num_n)

!write(*,*) "----------------"
!write(*,*) table
end subroutine reconfiguration



!old code table(i)=j
!subroutine reconfiguration(n,w,table)
!use rand_num
!implicit none
!integer,intent(IN):: n
!real(kind=8),intent(IN):: w(n)
!integer,intent(OUT):: table(n)
!real(kind=8):: wsum,z
!real(kind=8),allocatable::p(:)

!real(kind=8)::cap
!real(kind=8):: ran
!
!
!allocate(p(n))
!ran=rndm()
!
!wsum=0d0
!do j=1,n
!   wsum=wsum+w(j)
!enddo
!
!!---Prepare the probability
!do j=1,n
!   p(j)=w(j)/wsum
!!print *, 'j,p', j,p(j)
!enddo
!
!
!!wsum=0.d0
!!cap=1.d0  !0.1d0
!!do j=1,n,1
!!   if(p(j).GT.cap) then
!!     write(*,*) "cap the p(j)",j,p(j),cap
!!     p(j)=cap
!!   end if
!!   wsum=wsum+p(j)
!!enddo
!!do j=1,n
!!   p(j)=p(j)/wsum
!!enddo
!
!
!
!!---initialize-----------
!j=1
!wsum=p(1)
!
!do i=1,n
!   z=dble((i-1)+rndm())/dble(n)
!   !z=dble((i-1)+ran)/dble(n)
!   do while(z>wsum)
!      !print *, i,j,z,wsum
!      !pause
!      j=j+1
!      wsum=wsum+p(j)
!   enddo
!   table(i)=j
!enddo
!
!deallocate(p)
!end subroutine reconfiguration



!-------------------------------------------------------------
!This subroutine decides whether to do the pop control or not.
!-------------------------------------------------------------
subroutine assess_walkers(reconfig)
use param
use mc_loop_param
use phi_param
use adET_param
use mpi_serial_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
logical,intent(OUT):: reconfig
real(kind=8)::  maxw, minw, meanw
integer:: i
real(kind=8)::cc
real(kind=8)::wt(Nwalkers)


!If cpmc then only use weight
!if fp,rcp we use the ccoe
if(crn.GT.0) then
  cc=ccoe
else
  cc=0.d0
end if
do i=1,Nwalkers,1
   wt(i)=weight(i)*(cc*abs(tot_imp(i)/abs_imp_avg)+1.d0-cc)
end do



#ifdef MPI

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(wt(1),Nwalkers,MPI_DOUBLE_PRECISION,Mweight(1), &
               & Nwalkers,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
 if(rank.eq.0) then

   meanw=0d0
   do i=1,MNwalkers,1
      if (Mweight(i)<0d0) then
         print *, 'negative weight', Mweight(i),i
         call mystop !This stop is not safe, it only stop the rank.eq.0 thread.
      endif
      meanw=meanw+Mweight(i)
   enddo
   meanw=meanw/dble(MNwalkers)

   maxw=maxval(Mweight(:),1)
   minw=minval(Mweight(:),1)

   if (maxw>2d0 .or. minw<0.5d0)  then
      reconfig=.true.
      m_w=meanw
   else
      reconfig=.false.
      m_w=1.d0
   endif
 end if

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_BCAST(m_w,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(reconfig,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)


#else



 meanw=0d0
 do i=1,Nwalkers,1
    if (wt(i)<0d0) then
       print *, 'negative weight', weight(i),wt(i)
       call mystop
    endif
    meanw=meanw+wt(i)
 enddo
 meanw=meanw/dble(Nwalkers)


 maxw=maxval(wt(:),1)
 minw=minval(wt(:),1)

 if (maxw>2d0 .or. minw<0.5d0)  then
    reconfig=.true.
    m_w=meanw
 else
    reconfig=.false.
    m_w=1.d0
 endif

#endif

!DEBUG
!if(rank.eq.0) write(*,*) "The mean value of the total weight*c*imp/imp_avg+weight*(1.d0-c) before pop control:"
!if(rank.eq.0) write(*,*) "It should be 1 if after the pop contrl."
!if(rank.eq.0) write(*,*) meanw,reconfig
!if(rank.eq.0) write(*,*) ""
!if(rank.eq.0) write(*,*) ""
!if(rank.eq.0) write(*,*) ""

end subroutine assess_walkers


!----------------------------------------------------------------------------------
!Give in l between 1~MNwalkers, give out k between 0~(Nsize-1),l is the number in k
!----------------------------------------------------------------------------------
subroutine nthd(l,Nwalkers,k)
integer,intent(INOUT)::l
integer,intent(OUT)::k
integer,intent(IN)::Nwalkers
integer::j
j=l-1
k=j/Nwalkers
l=l-k*Nwalkers
end subroutine nthd

