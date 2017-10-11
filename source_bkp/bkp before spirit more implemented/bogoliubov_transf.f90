subroutine bogoliubov_transf(kpointindex,hamiltonian,evalues,leftevector,rightevector)
   use mod_global
   implicit none
   integer, intent(in) :: kpointindex
   complex(kind=pc), intent(inout) :: hamiltonian(twonaucell,twonaucell), evalues(twonaucell), leftevector(twonaucell,twonaucell), rightevector(twonaucell,twonaucell)
   complex(kind=pc) :: auxevector(twonaucell,twonaucell), diagonal(twonaucell,twonaucell), number
   real(kind=pc) :: summation
   character(len=10) :: format = 'f14.10'
   integer :: i, j
   
   leftevector = trimatmul( g , trp_cjg(rightevector), g )

   return
 
   ! Renormalizing the eigenvectors
   diagonal = matmul( leftevector , rightevector )

   do j=1, twonaucell
      if( abs(diagonal(j,j)) .GE. zero_toler ) then
         number = 1.0d0/sqrt(diagonal(j,j))
      else
         number = 1.0d0
      end if
      rightevector(:,j) = number*rightevector(:,j)
      leftevector(j,:) = number*leftevector(j,:)
   end do
   
   hamiltonian(naucell+1:twonaucell,:) = -hamiltonian(naucell+1:twonaucell,:)

   if(toprint) then
      auxevector(:,1) = real(evalues)
      auxevector(:,2) = aimag(evalues)
      call printmatrix("Eigenvalues of Dmatrix"//new_line('a')//"  Re          Im", real(auxevector(:,1:2)), format)
      call printmatrix("Rightevector", real(rightevector), format)

      diagonal = matmul( leftevector, rightevector )
      call printmatrix("L.R", real(diagonal), format)

      diagonal = trimatmul( trp_cjg(rightevector) , hamiltonian, rightevector )
      call printmatrix("R+.H.R",  real(diagonal), format)

      diagonal = trimatmul( trp_cjg(rightevector), g, rightevector)
      call printmatrix("R+.g.R",  real(diagonal), format)
      diagonal = trimatmul(rightevector , g, trp_cjg(rightevector))
      call printmatrix("R.g.R+",  real(diagonal), format)

      pause
      return
   end if

   ! R+.H.R
   diagonal = trimatmul( trp_cjg(rightevector) , hamiltonian, rightevector )
   do i=1, twonaucell
      if( abs( diagonal(i,i) - abs(evalues(i)) ) > zero_toler ) then
         print *, "Bogoliubov Transformation failed. Diagonal of R+.H.R is not equal to the eigenvalues:" 
         print *, "i: ", i, " diagonal(i,i)= ", diagonal(i,i),  " |evalues(i)|= ", abs(evalues(i)), " k point = ", kpointindex, ". Stopping..."
         stop
      end if
      diagonal(i,i) = 0.d0
   end do

   summation = sum(abs(diagonal))
   if( summation > zero_toler ) then
      print '(a,f12.8,a,i0,a)', "Bogoliubov Transformation failed. R+.H.R is not diagonal. Sum of off-diagonal elements= ", summation, " > 0  |  k point = ", kpointindex, ". Stopping..."
      stop
   end if

   ! R.g.R+
   diagonal = trimatmul(rightevector , g, trp_cjg(rightevector))
   summation = sum(abs(diagonal-g))
   if( summation > zero_toler ) then
      print '(a,f12.8,a,i0,a)', "Bogoliubov Transformation failed: R^+.g.R - g = ", summation, " > 0  |  k point = ", kpointindex, ". Stopping..."
      stop
   end if

   return

end subroutine bogoliubov_transf