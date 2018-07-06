module ising
contains
 function viscosity(temperature)
 real:: viscosity,temperature
    viscosity=temperature
 end function

subroutine initial(state) 
 integer:: state(5,5)
 state=1
end subroutine initial

subroutine get_energy(state,energy)
integer:: state(5,5)
real:: energy,U
U=1
energy=0.0
   do i=1,5
      do j=1,5
         energy=energy+U*state(i,j)*state(i,j)
      enddo
   enddo
end subroutine


end module
