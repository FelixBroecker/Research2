program davidson

  implicit none
 

  intrinsic                 ::  selected_real_kind, abs
  integer,  parameter       ::  wp = selected_real_kind(15)
  real(wp)                  ::  start_david, end_david, start_lapack, end_lapack, zero
  real(wp), allocatable     ::  mat(:,:), diagonal(:), eigenvals_lap(:), eigenvecs_lap(:,:), eigenvals_dav(:), eigenvecs_dav(:,:)
  integer                   ::  i, j, ndim, eigen_in, verbose
  real(wp)                  ::  dnrm2


  ndim          = 12
  eigen_in      = 4
  verbose       = 1
  
! allocate space for matrix
  allocate(mat(ndim, ndim), diagonal(ndim),)
! allocate space for eigenvecs and eigenvals
  allocate(eigenvecs_lap(ndim, ndim), eigenvals_lap(ndim), eigenvecs_dav(ndim, eigen_in), eigenvals_dav(eigen_in))


  zero = 0.0d0
  mat = zero 
  diagonal = zero
  eigenvecs_lap = zero
  eigenvals_lap = zero


! get matrix
  do i = 1, ndim
    do j = i, ndim
      if (j == i) then
        mat(i,j) = i + 1.0d0
        diagonal(i) = mat(i,j)
      else if (j > i) then
        mat(i,j) = 1.0d0 / (dble(i + j))
        mat(j,i) = mat(i,j)
      end if
    end do
  end do

  print *
  print *, '----------------'
  print *, 'Input Matrix A'
  print *, '----------------'
  print *
  call printMatrix(mat, ndim, ndim)
  print *

  call cpu_time(start_david)
  call symmetricDavidson(mat, ndim, eigenvecs_dav, eigenvals_dav, eigen_in, verbose)
  call cpu_time(end_david)
  
  print *
  print *, '----------------'
  print *, 'Results Davidson'
  print *, '----------------'

  print *
  print *, 'Eigenvalues:'
  print *
  call printVector(eigenvals_dav, eigen_in)
  print *
  print *, 'Eigenvectors'
  print *
  call printMatrix(eigenvecs_dav, eigen_in, ndim)

  eigenvecs_lap = mat
  call cpu_time(start_lapack)
  call lapackDiag(eigenvecs_lap, eigenvals_lap, ndim)
  call cpu_time(end_lapack)

  print *
  print *, '----------------'
  print *, 'Results Lapack'
  print *, '----------------'
  print *
  print *, 'Eigenvalues:'
  print *
  call printVector(eigenvals_lap, ndim)
  print *
  print *, 'Eigenvectors'
  print *
  call printMatrix(eigenvecs_lap, ndim, ndim)
  print * 

  print*, 'Difference in eigenvalues:'
  do i = 1, eigen_in
    print *, eigenvals_dav(i) - eigenvals_lap(i)
  end do
  print * 

  print *, 'Difference in eigenvalue'
  do i= 1, eigen_in
    print *, abs(dnrm2(eigen_in, eigenvecs_dav(:,i), 1) - dnrm2(eigen_in, eigenvecs_lap(:,i), 1))
  end do
  
  

  print*
  print*, '------------------------------------------'
  print * , 'overall wall time Davidson', end_david - start_david, 's'
  print * , 'overall wall time Lapack', end_lapack - start_lapack, 's'

  deallocate(mat, diagonal, eigenvecs_lap, eigenvecs_dav, eigenvals_dav, eigenvals_lap)


contains
  subroutine symmetricDavidson(mat_in, dim_mat_in, return_eigenvecs, return_eigenvals, eigen_in, verbose)

    real(wp), intent(in)      ::  mat_in(:,:)
    real(wp), intent(out)     ::  return_eigenvecs(:,:), return_eigenvals(:)
    integer, intent(in)       ::  dim_mat_in, eigen_in, verbose
    real(wp)                  ::  lw(1), threshold_residual, zero, check_GS, thresh_GS
    real(wp), allocatable     ::  matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
    real(wp), allocatable     ::  ritzVector(:,:), temp_mat(:,:), tau(:)
    real(wp), allocatable     ::  residual(:,:), temp_mat_prime(:,:),  temp(:,:)
    integer                   ::  i, j, it, ndimA, ndimV, maxiter, idxMaxVal(1), lwork, info,  n_grow,  max_orth
    real(wp)                  ::  dnrm2, diff
    logical, allocatable      ::  mask(:), converged(:)
    logical                   ::  matrix_not_full, GS_in_loop, not_orthogonal



    ndimA                 = dim_mat_in
    ndimV                 = 20
    maxiter               = 30
    max_orth              = 5
    threshold_residual    = 1.d-3
    thresh_GS             = 1.d-5
    GS_in_loop            = .false.



!   allocate space to create matrix A
    allocate(matA( ndimA, ndimA), mask(ndimA), diagonalA(ndimA))

!   allocate space to obtain reduced space
    allocate(matV( ndimA, ndimV ), matW(ndimA, ndimV), matP(ndimV, ndimV))

!   allocate space for diagonalization
    allocate(eigenvals(ndimV), eigenvecs(ndimV, ndimV))

!   allocate space for ritz vector and residual
    allocate(ritzVector(ndimA,ndimV), residual(ndimA, eigen_in))

!   allocate space for diagonalization
    allocate(temp_mat(ndimV, eigen_in), temp_mat_prime(ndimA, eigen_in), temp(ndimA, eigen_in))

!   allocate space for else
    allocate(converged(eigen_in))



    zero                  = 0.0d0
    
    matA            = zero
    diagonalA       = zero
    matW            = zero
    matP            = zero
    eigenvals       = zero
    eigenvecs       = zero
    ritzVector      = zero
    temp_mat        = zero
    temp_mat_prime  = zero
    residual        = zero
    temp            = zero
    matV            = zero
    mask            = .true.
    converged       = .false.


    if (verbose .ge. 1) then
      print *, '-- start of Davidson routine --'
    end if

!   get matrix diagonal
    do i = 1, dim_mat_in
      do j = i, dim_mat_in
        if (j .eq. i ) then
          diagonalA(i) = mat_in(i,j)
        end if
      end do
    end do

    matA = mat_in

    if (verbose .ge. 2) then
      write(*,*) 'Matrix A'
      call printMatrix(matA, ndimA, ndimA)
      print *
    end if


!   get initial vector
!   search for column with maximum value and use that column for initial vector

    do i = 1, eigen_in
      idxMaxVal = minloc(diagonalA, mask=mask)
      matV(idxMaxVal, i) = 1.0d0
      mask(idxMaxVal) = .false.
    end do 

    if (verbose .ge. 2) then
      write(*,*) 'Initial Vector'
      call printMatrix(matV, eigen_in, ndimA)
      print *
    end if


    n_grow = eigen_in

!   start loop

    outer: do it = 1, maxiter
      if (verbose .ge. 2) then
        print *, '------------------------------'
        print *, 'Iteration', it
        print *, '------------------------------'
        print *
      end if

      eigenvecs = 0.0d0
      eigenvals = 0.0d0

!     get projection matrix P = (W)^T * V;  W = A * V

      call dgemm('n', 'n', ndimA, n_grow, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
      call dgemm('t', 'n', ndimV, n_grow, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)


      if (verbose .ge. 3) then
        print *, 'Matrixproduct W (A V) :'
        call printMatrix(matW, n_grow, ndimA)
        print *
        write(*,*) 'Matrixproduct (W)T V = P'
        call printMatrix(matP, n_grow, n_grow)
      end if



!     diagonalize and obtain eigenvalues and eigenvectors 

      eigenvecs = matP
      call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, lw, -1, info)
      lwork = int(lw(1))
      allocate(work(lwork))
      call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, work, lwork, info)
      deallocate(work)

      if (verbose .ge. 2) then
        print *
        print *, 'Eigenvalues:'
        call printVector(eigenvals, n_grow)
        print *
        print *, 'Eigenvectors:'
        call printMatrix(eigenvecs, n_grow, n_grow)
      end if

!     get ritz vector 

      call dgemm('n', 'n', ndimA, n_grow, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)

      if (verbose .ge. 2) then
        print *
        print *, 'Ritzvector all'
        call printMatrix(ritzVector, n_grow, ndimA)
        print *
      end if

!     get residual

      call dgemm('n', 'n', ndimA, eigen_in, ndimV, 1.0d0, matW, ndimA, eigenvecs, ndimV, 0.0d0, residual, ndimA)
      do i = 1, eigen_in
        call daxpy(ndimA, -eigenvals(i), ritzVector(:,i), 1, residual(:,i), 1)
      end do

      if (verbose .ge. 3) then
        print *
        print *, 'Residual'
        call printMatrix(residual, eigen_in, ndimA)
        print *
      end if


!     compute norm and check convergency

      do i = 1, eigen_in
        if (verbose .ge. 2) then
          print *, 'residual norm of eigenvector', i, ':', dnrm2(ndimA, residual(:,i), 1)
          print *
        end if
        if(dnrm2(ndimA, residual(:,i), 1) .le. threshold_residual) then
          if (verbose .ge. 2) then
            print *, 'Converged for:', i
          end if
          converged(i) = .true.
        end if
      end do

      if (all(converged)) then
        if (verbose .ge. 1) then
          print *, 'Converged for all sought eigenpairs after', it, 'iterations.'
          print *
          print *, 'Eigenvalues:'
          print *
          call printVector(eigenvals, eigen_in)
          print *
          print *, 'Eigenvectors:'
          print *
          call printMatrix(eigenvecs, eigen_in, n_grow)
        end if


!       copy eigenpairs in output

        do i = 1, eigen_in 
          return_eigenvals(i) = eigenvals(i)
          return_eigenvecs(:,i) = ritzVector(:,i)
        end do
       
        exit outer
      end if
      

!     check if matrix is not full

      matrix_not_full = .true.
      do i = 1, ndimA
        if (matV(i, ndimV + 1 - eigen_in) /= 0.0d0) then
          matrix_not_full = .false.   
          exit
        end if
      end do


!     precondition y_i = r_i / (D - lambda_i)

      if (matrix_not_full) then
        do i = 1, eigen_in
          do j = 1, ndimA
            diff = diagonalA(j) - eigenvals(1)
            if (diff >= 1.d-1) then
              residual(j,i) = residual(j,i) / (diff)
            else
              residual(j,i) = 1.d-1
            end if
          end do
        end do
        

!       Gram Schmidt orthogonalization

!       loop implementation

        if (GS_in_loop) then
          do i = 1, eigen_in
            temp = 0.0d0
            do j = 1, n_grow
              temp(:,i) = temp(:,i) + dot_product(residual(:,i), matV(:,j)) / dot_product(matV(:,j), matV(:,j)) * matV(:,j)
            end do
            residual(:,i) = residual(:,i) - temp(:,i)  
            residual(:,i) = residual(:,i) / dnrm2(ndimA, residual(:,i), 1)
          end do 

!       check if orthogonal

          not_orthogonal = .false.
          call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS, not_orthogonal, .true.)
        
            
!       matrix implementation

        else
          not_orthogonal = .false.
          do j = 1, max_orth
            if (.not. not_orthogonal) then

!             matrix product V^T * y  
              call dgemm('t', 'n',  n_grow, eigen_in ,ndimA, 1.0d0, matV, ndimA, residual, ndimA, 0.0d0, temp_mat, ndimV)

!             matrix product V * (V^T * y)
              call dgemm('n', 'n', ndimA, eigen_in, n_grow, 1.0d0, matV, ndimA, temp_mat, ndimV, 0.0d0, temp_mat_prime, ndimA)

!             y_prime = y - V * (V^T * y) 
              do i = 1, eigen_in
                residual(:, i) = residual(:, i) - temp_mat_prime(:, i)
              end do

              if (verbose .ge. 3) then
                print *
                print *, 'y_i: Matrix product V^T * y'
                call printMatrix(temp_mat, eigen_in, n_grow)
                print *
                print *, 'y_i: Matrix product V * ( V^T * y )'
                call printMatrix(temp_mat_prime, eigen_in, ndimA)
                print *
                print *, 'Orthogonalized precondition y'
                call printMatrix(residual, eigen_in, ndimA)
              end if

!             print a warning if matrix is not orthogonal 

              call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS, not_orthogonal, .true.)
            else
              exit
            end if
          end do
            


!         orthogonalize residual with itself

          if (eigen_in .gt. 1) then

            allocate(tau(eigen_in))
            tau = zero
            call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, lw, -1, info)
            lwork = int(lw(1))
            allocate(work(lwork))
            work = zero
            call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, work, lwork, info)
            call checkInfo(info, 'Orthogonalization of residual step 1')


            call dorgqr(ndimA, eigen_in, eigen_in, residual, ndimA, tau, work, lwork, info)
            call checkInfo(info, 'Orthogonalization of residual step 2')
            deallocate(work)
            deallocate(tau)


!           check if orthogonal

            not_orthogonal = .false.
            call checkOrth1mat(residual, eigen_in, 'Precondition', thresh_GS, not_orthogonal, .true.)
            if (not_orthogonal) then
              print *, '---WARNING --- RESIDUAL NOT ORTHOGONAL'
            end if
          end if

          

!         orthonormalization
          do i = 1, eigen_in
            residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
          end do

          if (verbose .ge. 2) then
            print *
            print *, 'Orthonormalized precondition y'
            call printMatrix(residual, eigen_in, ndimA)
          end if

        end if  




!       add orthonormal residual to subspace

        do i = 1, eigen_in
          matV(:, n_grow  + i) = residual(:,i)
        end do

        if (verbose .ge. 3) then
          print *
          print *, 'New subspace V'
          call printMatrix(matV, n_grow + eigen_in, ndimA)
          print *
          print *
        end if 


!       handle case if not converged after maxiter steps

        if (it .eq. maxiter) then
          if (verbose .ge. 1) then
            print *
            print *, 'Not converged after', it, 'iterations'
          end if
          return_eigenvecs = zero
          return_eigenvals = zero
        end if


!       restart
        n_grow = n_grow + eigen_in 


!       treat case if matrix full
      else

        matV            = zero

        do i = 1, eigen_in
          matV(:, i) = ritzVector(:, i)
        end do

        if (verbose .ge. 2) then
          print *
          print *, 'Matrix V full'
          print *
          print *, 'New subspace'
          call printMatrix(matV, eigen_in, ndimA)
          print *
          print *
        end if

!       restart

        matW            = zero
        matP            = zero
        eigenvals       = zero
        eigenvecs       = zero
        ritzVector      = zero
        temp_mat        = zero
        temp_mat_prime  = zero
        residual        = zero
        temp            = zero

        n_grow          =  eigen_in
        
      end if


    end do outer

    deallocate(matA, mask, diagonalA, matV, matW, matP, eigenvals, eigenvecs, &
               ritzVector, residual, temp_mat, temp_mat_prime, temp, converged)


    if (verbose .ge. 1) then
      print *
      print *, '-- end of Davidson routine --'
      print *
    end if

  end subroutine symmetricDavidson


  subroutine lapackDiag(mat, eigenvals, ndimMat)
! get diagonalized matrix with lapack function

    real(wp),   intent(inout)     :: mat(:,:)
    real(wp),   intent(out)       :: eigenvals(:)
    integer,    intent(in)        :: ndimMat
    real(wp)                      :: lw(1)
    real(wp), allocatable         :: work(:)
    integer                       :: lwork, info

    call dsyev('V', 'u', ndimMat, mat, ndimMat, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'u', ndimMat, mat, ndimMat, eigenvals, work, lwork, info)
    call checkInfo(info, 'diagonalize whole mat A')
    deallocate(work)

  end subroutine lapackDiag


  subroutine printMatrix(mat, nrows, ncols) 
! print formatted matrix

    real(wp), intent(in)  :: mat(:,:)
    integer , intent(in)  :: nrows, ncols
    integer :: i,j 

    do i = 1, ncols
      print *, (mat(i,j), j= 1, nrows )
    end do

  end subroutine printMatrix


  subroutine printVector(vec, lenRow)
! print formatted vector

    real(wp), intent(in)  :: vec(:)
    integer,  intent(in)  :: lenRow
    integer               :: i

    do i = 1, lenRow
      print *, vec(i)
    end do

  end subroutine printVector


  subroutine checkInfo(info, occasion)
! check if info is zero and print error message if not

    integer               :: info, zero
    character(len=*)     :: occasion

    zero = 0
    if (info .NE. zero) then
      print *
      print *, '--- WARNING ---'
      print *, occasion
      print *, 'Process terminated with info not equal to 0'
      print*
    end if

  end subroutine checkInfo


  subroutine checkOrth2mat(mat1, nrows1, mat1_name, mat2, nrows2, mat2_name, thresh, not_orthogonal, verbose)
! check if vectors of two matrices are orthogonal

    integer,              intent(in)      ::  nrows1, nrows2
    real(wp),             intent(in)      ::  thresh, mat1(:,:), mat2(:,:)
    character(len=*),    intent(in)       ::  mat1_name, mat2_name
    real(wp)                              ::  dot_prod
    integer                               :: i,j 
    logical,              intent(inout)   :: not_orthogonal
    logical,              intent(in)      :: verbose

    do i = 1, nrows1
      do j = 1, nrows2
        dot_prod = abs(dot_product( mat2(:,j), mat1(:,i) ))
        if ( dot_prod .GT. thresh ) then
          not_orthogonal = .true.
          if (verbose) then
            print *
            print *, '--- WARNING ---'
            print *, 'vector ', i, ' of matrix', mat1_name, 'is not orthogonal to vector', j, 'of matrix', mat2_name
            print *, 'Result of dot product:', dot_prod
            print *
          end if 
        end if
      end do
    end do

  end subroutine checkOrth2mat


  subroutine checkOrth1mat(mat1, nrows1, mat1_name, thresh, not_orthogonal, verbose)
! check if the vectors within the matrix are orthogonal

    integer,              intent(in)      ::  nrows1
    real(wp),             intent(in)      ::  thresh, mat1(:,:)
    character(len=*),     intent(in)      ::  mat1_name
    logical,              intent(inout)   ::  not_orthogonal
    logical,              intent(in)      :: verbose
    real(wp)                              ::  dot_prod
    integer :: i,j 

    do i = 1, nrows1
      do j = 1, nrows1
        if (j .NE. i) then
          dot_prod  = abs(dot_product( mat1(:,j), mat1(:,i) ))
          if ( dot_prod .GT. thresh ) then
            not_orthogonal = .true.
            if (verbose) then 
              print *
              print *, '--- WARNING ---'
              print *, 'vector ', i, ' of matrix ', mat1_name, ' is not orthogonal to vector ', j, ' of matrix', mat1_name
              print *, 'Result of dot product: ', dot_prod
              print *
            end if
          end if
        end if
      end do
    end do

  end subroutine checkOrth1mat
 
end program davidson
