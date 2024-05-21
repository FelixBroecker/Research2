program davidson
 use diagonalization
  implicit none
!
!
!
!                        0====================0
!                   0===//   0============0   \\===0
!               0==//   0===//   0====0   \\===0   \\==0
!           0==//  0===//    0==//    \\==0    \\===0  \\==0
!       ::://  0==//   0====//            \\====0   \\===0 \\:::
!    =========// 0====//                        \\====0 \\=========
!   ============//                                    \\============ 
!  ============x                                       x============
!   ============\\                                    //============ 
!    =========\\ 0====\\                        //====0 //=========
!       :::\\  0==\\   0====\\            //====0   //==0  //:::
!           0==\\  0===\\    0==\\    //==0    //===0  //==0
!               0==\\   0===\\   0====0   //===0   //==0 
!                   0===\\   0============0   //===0 
!                        0====================0 
!
!
!
!
!
 

  intrinsic                 ::  selected_real_kind, abs
  integer,  parameter       ::  wp = selected_real_kind(15)
  real(wp)                  ::  zero, threshold_residual
  real(wp), allocatable     ::  mat(:,:)
  integer                   ::  i, j, ndim, eigen_in, verbose, max_iter
  real(wp)                  ::  dnrm2

  

  ndim                  = 5
  zero                  = 0.0d0
!
! get matrix
!
!
! allocate space for matrix
!
  allocate(mat(ndim, ndim))
  mat = zero 
!
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


!
!  call cholesky decomposition
!
  call testCholesky(ndim, ndim, mat, .false.)

!
!  settings for davidson and call routine
!
  eigen_in              = 4
  verbose               = 1
  max_iter              = 100
  threshold_residual    = 1.d-7

  call testDavidson(ndim, mat, eigen_in, threshold_residual, verbose, max_iter)


  deallocate(mat)


contains
  subroutine testDavidson(ndim, mat, eigen_in, threshold_residual, verbose, max_iter)
!    
    implicit none
    intrinsic                 ::  selected_real_kind
    integer,  parameter       ::  wp = selected_real_kind(15)
!
    integer, intent(in)       :: ndim, eigen_in, verbose, max_iter
    real(wp), intent(in)          :: mat(:,:), threshold_residual
!
    real(wp)                  ::  start_david, end_david, start_lapack, end_lapack, zero
    real(wp), allocatable     ::  eigenvals_lap(:), eigenvecs_lap(:,:), eigenvals_dav(:), eigenvecs_dav(:,:)
    real(wp)                  ::  dnrm2

    
!   allocate space for eigenvecs and eigenvals
    allocate(eigenvecs_lap(ndim, ndim), eigenvals_lap(ndim), eigenvecs_dav(ndim, eigen_in), eigenvals_dav(eigen_in))


    zero = 0.0d0
    eigenvecs_lap = zero
    eigenvals_lap = zero

    call cpu_time(start_david)
    call symmetricDavidson(mat, ndim, eigenvecs_dav, eigenvals_dav, eigen_in, verbose, threshold_residual, max_iter)
    call cpu_time(end_david)
    
    print *
    print *, '----------------------------------------------------------------'
    print *, '                        Results Davidson'
    print *, '----------------------------------------------------------------'

    print *
    print *, 'Eigenvalues:'
    print *
    call printVector(eigenvals_dav, eigen_in)

    eigenvecs_lap = mat
    call cpu_time(start_lapack)
    call lapackDiag(eigenvecs_lap, eigenvals_lap, ndim)
    call cpu_time(end_lapack)

    print *
    print *, '----------------------------------------------------------------'
    print *, '                         Results Lapack'
    print *, '----------------------------------------------------------------'
    print *
    print *, 'Eigenvalues:'
    print *
    call printVector(eigenvals_lap, eigen_in)

    print *
    print *, '----------------------------------------------------------------'
    print *, '                          Comparison'
    print *, '----------------------------------------------------------------'
    print *

    print*, 'Difference in eigenvalues:'
    do i = 1, eigen_in
      print '(E13.3)', eigenvals_dav(i) - eigenvals_lap(i)
    end do


    print*
    print*, '================================================================='
    print*
    print '(A, E13.3, A)' , 'overall wall time Davidson', end_david - start_david, ' s'
    print '(A, E13.3, A)' , 'overall wall time Lapack  ', end_lapack - start_lapack, ' s'
    print*
    print*, '================================================================='  
    print*

    deallocate(eigenvecs_lap, eigenvecs_dav, eigenvals_dav, eigenvals_lap)
!    
  end subroutine testDavidson
!
!
!
subroutine testCholesky(n_row, n_col, mat, do_test)
!    
    implicit none
    intrinsic                 :: selected_real_kind
    integer,  parameter       :: wp = selected_real_kind(15)
!
    real(wp), intent(inout)      :: mat(:,:)
    integer,  intent(in)      :: n_row, n_col
    logical,  intent(in)      :: do_test
    real(wp), allocatable     :: test(:,:)
!
if (do_test) then
!      
      allocate(test(3,3))
      test(1,:)=(/ 25,  15,  -5 /)
      test(2,:)=(/ 15,  18,   0 /)
      test(3,:)=(/ -5,   0,  11 /)
!
      print *
      print *, 'Cholesky Decomposition:' 
      print *
      print *, 'Input Matrix:'
      print *
      call printMatrix(test, 3, 3)
      print *
      print *, 'result Cholesky routine'
      print *
      call cholesky(3,3,test)
    else
!      
      print *
      print *, 'Cholesky Decomposition:' 
      print *
      print *, 'Input Matrix:'
      print *
      call printMatrix(mat, n_row, n_col)
      print *
      print *, 'result Cholesky routine'
      print *
      call cholesky(n_row, n_col, mat)
    
    end if 


  end subroutine

end program davidson
