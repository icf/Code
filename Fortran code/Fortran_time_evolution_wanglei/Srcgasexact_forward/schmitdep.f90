
!--- Schmit orthogonalization---------------------------------      
do l=1,m
           ortho=0d0
            do j=1, l-1
               do k=1,m
               ortho(j)=ortho(j)+phi_out(n, k)*conjg(phi_(j,k))
               enddo 
           enddo

    do k=1,m
      do j=1,n
      phi_out(j,k)=phi_out(j,k)-phi_out(n,k)*ortho(j)
     enddo
    enddo 

enddo 


