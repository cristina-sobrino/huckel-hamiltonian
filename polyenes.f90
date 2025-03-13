program diag
      implicit none
      real*8, allocatable  :: mat(:,:), eigen_values(:), eigenvector(:,:), beta_plus, beta_minus !mat is the huckel matrix
      real*8 :: theta, x, alpha_1, alpha_2
      real*8 :: pi = ACOS(-1.0)
      real*8, allocatable :: e_minus(:), e_plus(:)
      integer               :: i, j, n_at, contador, k
      character(len=100) :: mol_type, bond_type, atom_type

      beta_plus = 0.d0
      beta_minus = 0.d0

      open(1, file="input.inp")
      read(1,*) bond_type, beta_plus, beta_minus ! alternated bonds or equal bonds
      read(1,*) mol_type ! linear or cyclic
      read(1,*) atom_type ! all atoms are equal or not
      read(1,*) n_at
      close(1)



! CLASSIFICATION OF THE MOLECULE 

! First case: can be either linear/chain or a ring
      if (mol_type == "linear" .or. mol_type == "ring")  then
         write(*,*) "The molecule is ", mol_type
      end if

! Second case: all the bonds can be equal or they can be alternated with two different beta values
      if (bond_type == "equal") then
      write(*,*) "All the bonds of the molecule are equal"
      beta_plus = -1.d0
      beta_minus = -1.d0

      else if (bond_type == "alternated") then
              write(*,*) "Two different bond lengths"
           if (beta_plus == 0.d0 .and. beta_minus == 0.d0) then
                   write(*,*) "Please select a value for beta_plus and beta_minus"
           end if
           
      else
          write(*,*) "Unrecognised keyword. The two possible keywords are:"
          write(*,*) "Equal for systems with equal bonds"
          write(*,*) "Alternated for two different bond lengths"
        end if


! Third case: all the atoms can be equal or they can be alternated
       if (atom_type == "alt_atoms" ) then
            write(*,*) "The system has two different types of atoms"
            alpha_1 = 0.d0
            alpha_2 = 1.d0 
      else if (atom_type == "eq_atoms" ) then 
              write(*,*) "All the atoms are equal"
        alpha_1 = 0.d0
        alpha_2 = 0.d0

      else
        write(*,*) "Unrecognized keyword for atom type. "
        write(*,*) " Keywords: alt_atoms for two different types of atoms"
        write(*,*) "           eq_atoms for only type of atom"
      end if



! Filling the Hamiltonian matrix with the appropriate values   
        contador = 0
 ! dimension of matrix
      allocate(mat(n_at, n_at), eigen_values(n_at), eigenvector(n_at, n_at))
      do i = 1,n_at
        do j = 1,n_at
          if (i == j) then

                if (mod(contador,2) == 0) then
                  mat(i,j) = alpha_1 !the diagonal elements are the alpha values
                  mat(i,j+1) = beta_plus !the off-diagonal elements are the beta values
                  mat(j+1,i) = beta_plus
                else 
                  mat(i,j) = alpha_2 
                  mat(i,j+1) = beta_minus
                  mat(j+1,i) = beta_minus
                end if
                contador = contador + 1

          end if
        end do
      end do
      
      !For the rings, we need do the hamiltonian as if we had a chain and 
      !then add beta values at the extreme positions to "close" the chain
  
      if (mol_type == "ring" ) then
         mat(n_at,1) = beta_minus
         mat(1,n_at) = beta_minus


       ! ANOTHER WAY OF GETTING THE EIGENVALUES FOR PBC

       allocate(e_minus(n_at/2), e_plus(n_at/2))
       theta = 2*pi / (n_at/2)

       !do k = 1, n_at/2 - 1

       e_plus = 0.d0

       do k = 1, n_at/2

          e_minus(k) = alpha_1 - sqrt( beta_plus**2 + beta_minus**2 + 2*beta_plus*beta_minus*cos(k*theta) )
          e_plus(k)  = alpha_1 + sqrt( beta_plus**2 + beta_minus**2 + 2*beta_plus*beta_minus*cos(k*theta) )
       end do


        contador  = 1

       open (25, file="analytical_eigenval.out")

       do i =1,n_at/2
          write(25, '(I6, f16.8)') contador, e_minus(i)
          contador = contador + 1
       end do

       do i = 1, n_at/2
          write(25, '(I6, f16.8)') contador, e_plus(i)
          contador = contador + 1
       end do

      close(25)


     end if

     open (15, file="huckel_hamiltonian.out" )
      do i = 1,n_at
        write(15,'(10f16.4)') mat(i,:)
      end do
     close(15)

      eigenvector(:,:) = mat(:,:)  ! to store the matrix because it is modified when calling the function

      call diagonalize_matrix(n_at,mat,eigen_values)






      ! Writing the eigenvalues and eigenvectors
      open (2, file="eigen_values.out")
      do i =1,n_at
       write(2, '(I6, f16.8)') i, eigen_values(i)
      end do
      close(2)


        open(3, file="eigen_vectors.out")
        do i = 1,n_at
                write(3, '(I6, 100(f16.8))') i,mat(i,:)
        end do 
        close(3)


        deallocate(mat, eigen_values, eigenvector)
end program

subroutine diagonalize_matrix(N,A,e)


        implicit none

        ! Input variables

        integer, intent(in) :: N !dimension
        double precision,intent(inout) :: A(N,N) !matrix to diagonalize
        double precision,intent(out) :: e(N) !eigevalues

        ! local variables

        integer :: lwork, info
        integer :: i
        double precision, allocatable :: work(:)

        allocate(work(3*N))
        lwork = size(work)

        call dsyev('V','U',N,A,N,e,work,lwork,info)

        if(info /= 0) then
                write(*,'(a)') 'Problem in diagonalize_matrix (dsyev)!!'
         stop
         endif


do i = 1,N
        if (abs(e(i)) < 1e-10) e(i) = 0
end do

end subroutine diagonalize_matrix
