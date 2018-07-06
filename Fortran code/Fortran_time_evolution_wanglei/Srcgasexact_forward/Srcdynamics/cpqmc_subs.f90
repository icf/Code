subroutine cal_overlap(n,m,phi_l,phi_r,overlap)
implicit none

integer,intent(IN):: n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: overlap

complex(kind=8):: tempmat(m,m)
integer:: i

tempmat= matmul(transpose(conjg(phi_l)), phi_r)
call caldet_c(m,tempmat,overlap)

end subroutine 

subroutine cal_green(n,m,phi_l,phi_r,gf)
implicit none

integer,intent(IN):: n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: gf(n,n)

integer:: i,j
integer:: k,l
complex(kind=8):: tempmat(m,m)

tempmat= matmul(transpose(conjg(phi_l)), phi_r)
call  inverse(tempmat,m)


gf=0d0
    do i=1,n
      do j=1,n
         if (i==j) then 
            gf(i,j)=1d0
           else
            gf(i,j)=0d0
         endif

          do k=1,m
            do l=1,m
          gf(i,j)=gf(i,j) - phi_r(i,k)*tempmat(k,l) *conjg(phi_l(j,l))
            enddo 
         enddo 

       enddo 
      enddo 

end subroutine 


subroutine cal_occ(n,gf,occ)
implicit none
integer,intent(IN):: n
complex(kind=8),intent(IN):: gf(n,n)
real(kind=8), intent(OUT):: occ(n)

integer:: i

do i=1,n
occ(i)=1d0-dble(gf(i,i))
enddo 

end subroutine 

