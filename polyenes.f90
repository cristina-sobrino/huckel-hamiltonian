program diag
      implicit none
      real*8, allocatable  :: mat(:,:), eigen_values(:), eigenvector(:,:), alpha_1, alpha_2, beta_plus, beta_minus !mat is the huckel matrix
      integer               :: i, j, n_at, contador
      character(len=100) :: mol_type, bond_type, atom_type

      open(1, file="input.inp")
      read(1,*) bond_type ! alternated bonds or equal bonds
      read(1,*) mol_type ! linear or cyclic
      read(1,*) atom_type ! all atoms are equal or not
      read(1,*) n_at
      close(1)

      
      if (mol_type == "linear" .or. mol_type == "ring")  then
         write(*,*) "The molecule is ", mol_type
      end if

      if (bond_type == "equal") then
      write(*,*) "All the bonds of the molecule are equal"
      beta_plus = -1.d0
      beta_minus = -1.d0

      else if (bond_type == "alternated") then
              write(*,*) "Two different bond lengths"
        beta_plus = -1.0d0
        beta_minus = -0.1d0

      else
          write(*,*) "Unrecognised keyword. The two possible keywords are:"
          write(*,*) "Equal for systems with equal bonds"
          write(*,*) "Alternated for two different bond lengths"
        end if

      if (atom_type == "alt_atoms" ) then
            write(*,*) "The system has two different types of atoms"
            alpha_1 = 0.d0
            alpha_2 = 2.d0 
      else if (atom_type == "eq_atoms" ) then 
              write(*,*) "All the atoms are equal"
        alpha_1 = 0.d0
        alpha_2 = 0.d0

      else
        write(*,*) "Unrecognized keyword for atom type. "
        write(*,*) " Keywords: alt_atoms for two different types of atoms"
        write(*,*) "           eq_atoms for only type of atom"
      end if


        contador = 0
      !n = 100  ! dimension of matrix
      allocate(mat(n_at, n_at), eigen_values(n_at), eigenvector(n_at, n_at))
      do i = 1,n_at-1
        do j = 1,n_at-1
          if (i == j) then
                if (mod(contador,2) == 0) then
                  mat(i,j) = alpha_1
                  mat(i,j+1) = beta_plus
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
      
      if (mol_type == "ring" ) then
         mat(n_at,1) = beta_minus
         mat(1,n_at) = beta_minus
      end if

      do i = 1,n_at
        write(*,'(10f16.4)') mat(i,:)
      end do

      eigenvector(:,:) = mat(:,:)  ! to store matrix because it is modified when calling the function

      call diagonalize_matrix(n_at,mat,eigen_values)

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
