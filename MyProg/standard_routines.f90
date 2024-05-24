program standardRoutines

  implicit none

  intrinsic                 ::  selected_real_kind, abs
  integer,  parameter       ::  wp = selected_real_kind(15)
  real(wp), allocatable     ::  mat(:,:) 
  integer                   ::  i, j, ndim
  real(wp)                  :: thresh_GS

  ndim          = 5
  thresh_GS     = 1.d-6

  allocate(mat(ndim, ndim),)


  ! get matrix
  do i = 1, ndim
    do j = i, ndim
      if (j == i) then
        mat(i,j) = i + 1.0d0
      else if (j > i) then
        mat(i,j) = 1.0d0 / (dble(i + j))
        mat(j,i) = mat(i,j)
      end if
    end do
  end do

  print *, '----------------'
  print *, 'Input Matrix A'
  print *, '----------------'
  print *
  call printMatrix(mat, ndim, ndim)
  print *

contains
  subroutine GramSchmidtOrth(mat1, ncol1, mat2, ncol2, nrow)

!   orthogonalizes the vectors of the matrix mat2 to a the matrix mat1 in a slow loop implementation
!   by projecting the each vector of mat2 with each vector in the matrix mat1.

    real(wp),  intent(in)       :: mat1(:,:) 
    real(wp),  intent(inout)    :: mat2(:,:)
    real(wp), allocatable       :: temp(:,:)
    real(wp)                    :: dnrm2
    integer,   intent(in)       :: ncol1, ncol2, nrow
    integer                     :: i, j

    allocate(temp(nrow, ncol1))
!     orthogonal projection
    do i = 1, ncol2
      temp = 0.0d0
      do j = 1, ncol1
        temp(:,i) = temp(:,i) + dot_product(mat1(:,j), mat2(:,i)) / dot_product(mat2(:,j), mat2(:,j)) * mat2(:,j)
      end do
      mat2(:,i) = mat2(:,i) - temp(:,i)

!     orthonormalization
      mat2(:,i) = mat2(:,i) / dnrm2(nrow, mat2(:,i), 1)
    end do 

!    print a warning if matrix is not orthogonal
      call checkOrth2mat(mat2, ncol2, 'Vectors to be orthogonalized', mat1, ncol1, 'Matrix with projection vectors', thresh_GS)

  end subroutine GramSchmidtOrth
!

  subroutine checkOrth2mat(mat1, nrows1, mat1_name, mat2, nrows2, mat2_name, thresh)

!   check if vectors of two matrices are orthogonal!

    integer,              intent(in)      ::  nrows1, nrows2
    real(wp),             intent(in)      ::  thresh, mat1(:,:), mat2(:,:)
    character(len=*),     intent(in)      ::  mat1_name, mat2_name
    real(wp)                              ::  dot_prod
    integer :: i,j !

    do i = 1, nrows1
      do j = 1, nrows2
        dot_prod = abs(dot_product( mat2(:,j), mat1(:,i) ))
        if ( dot_prod .GT. thresh ) then
          print *
          print *, '--- WARNING ---'
          print *, 'vector', i, 'of matrix', mat1_name, 'is not orthogonal to vector', j, 'of matrix', mat2_name
          print *, 'Result of dot product:', dot_prod
          print *
        end if
      end do
    end do
  end subroutine checkOrth2mat


  subroutine printMatrix(mat, nrows, ncols) 

    real(wp), intent(in)  :: mat(:,:)
    integer , intent(in)  :: nrows, ncols
    integer :: i,j 

    do i = 1, ncols
      print *, (mat(i,j), j= 1, nrows )
    end do

  end subroutine printMatrix
  
end program standardRoutines
