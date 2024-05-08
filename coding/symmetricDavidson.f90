program davidson

  implicit none
 

  intrinsic                 ::  selected_real_kind
  integer,  parameter       ::  wp = selected_real_kind(15)
  real(wp)                  ::  start_david, end_david, start_lapack, end_lapack, zero
  real(wp), allocatable     ::  eigenvecs(:,:), mat(:,:), diagonal(:), eigenvals(:)
  integer                   ::  i, j, ndim


  ndim = 5


  allocate(mat(ndim, ndim), diagonal(ndim), eigenvecs(ndim, ndim), eigenvals(ndim))
  zero = 0.0d0
  mat = zero 
  diagonal = zero
  eigenvecs = zero

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

  print *, 'matrix A'
  call printMatrix(mat, ndim, ndim)


  call cpu_time(start_david)
  call symmetricDavidson()
  call cpu_time(end_david)

  print *
  print *, 'Results Lapack'
  call cpu_time(start_lapack)
  eigenvecs = mat
  call lapackDiag(eigenvecs, eigenvals, ndim)
  print *, 'eigenvals'
  call printVector(eigenvals, ndim)
  call cpu_time(end_lapack)
  

  print*
  print*, '------------------------------------------'
  print * , 'overall wall time Davidson', end_david - start_david, 's'
  print * , 'overall wall time Lapack', end_lapack - start_lapack, 's'


contains
  subroutine symmetricDavidson()

  real(wp)                  ::  lw(1), threshold_residual, zero, check_GS, thresh_GS, tau
  real(wp), allocatable     ::  matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
  real(wp), allocatable     ::  ritzVector(:,:), temp_mat(:,:), ritzVectorTemp(:)
  real(wp), allocatable     ::  residual(:,:), temp_mat_prime(:,:), diff, temp(:,:)
  integer                   ::  i, j, it, ndimA, ndimV, maxiter, idxMaxVal(1), lwork, info, eigen_in,  n_grow
  real(wp)                  ::  dnrm2
  logical, allocatable      ::  mask(:), converged(:)
  logical                   ::  matrix_not_full, verbose, GS_in_loop



  ndimA                 = 5
  ndimV                 = 4
  maxiter               = 6
  eigen_in              = 1
  threshold_residual    = 1.d-4
  verbose               = .false.
  GS_in_loop            = .false.
  thresh_GS             = 1.d-6



! allocate space to create matrix A
  allocate(matA( ndimA, ndimA ), mask(ndimA), diagonalA(ndimA))

! allocate space to obtain reduced space
  allocate(matV( ndimA, ndimV ), matW(ndimA, ndimV), matP(ndimV, ndimV))

! allocate space for diagonalization
  allocate(eigenvals(ndimV), eigenvecs(ndimV, ndimV))

! allocate space for ritz vector and residual
  allocate(ritzVector(ndimA,ndimV), residual(ndimA, eigen_in))
  allocate(ritzVectorTemp(ndimA))

! allocate space for diagonalization
  allocate(temp_mat(ndimV, eigen_in), temp_mat_prime(ndimA, eigen_in), temp(ndimA, eigen_in))

! allocate space for else
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



! get matrix
  do i = 1, ndimA 
    do j = i, ndimA
      if (j == i) then
        matA(i,j) = i + 1.0d0
        diagonalA(i) = matA(i,j)
      else if (j > i) then
        matA(i,j) = 1.0d0 / (dble(i + j))
        matA(j,i) = matA(i,j)
      end if
    end do
  end do

  write(*,*) 'Matrix A'
  call printMatrix(matA, ndimA, ndimA)
  print *


! get initial vector
  ! search for colum with maximum value and use that column for initial vector

  do i = 1, eigen_in
    idxMaxVal = minloc(diagonalA, mask=mask)
    matV(idxMaxVal, i) = 1.0d0
    mask(idxMaxVal) = .false.
  end do 

  write(*,*) 'Initial Vector'
  call printMatrix(matV, eigen_in, ndimA)
  print *

  n_grow = eigen_in

! start loop
  outer: do it = 1, maxiter
      print *, '------------------------------'
      print *, 'Iteration', it
      print *, '------------------------------'
      print *

    eigenvecs = 0.0d0
    eigenvals = 0.0d0
  ! get projection matrix P = (W)^T * V;  W = A * V
    call dgemm('n', 'n', ndimA, n_grow, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
    call dgemm('t', 'n', ndimV, n_grow, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)


    if (verbose) then
      print *, 'Matrixproduct W (A V) :'
      call printMatrix(matW, n_grow, ndimA)
      print *
      write(*,*) 'Matrixproduct (W)T V = P'
      call printMatrix(matP, n_grow, n_grow)
    end if



  ! diagonalize and obtain eigenvalues and eigenvectors 
    eigenvecs = matP
    call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, work, lwork, info)
    deallocate(work)

    print *
    print *, 'Eigenvalues:'
    call printVector(eigenvals, n_grow)
    print *
    print *, 'Eigenvectors:'
    call printMatrix(eigenvecs, n_grow, n_grow)


  ! get ritz vector 
    call dgemm('n', 'n', ndimA, n_grow, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)

    print *
    print *, 'Ritzvector all'
    call printMatrix(ritzVector, n_grow, ndimA)
    print *


  ! get residual
    call dgemm('n', 'n', ndimA, eigen_in, ndimV, 1.0d0, matW, ndimA, eigenvecs, ndimV, 0.0d0, residual, ndimA)
    do i = 1, eigen_in
      call daxpy(ndimA, -eigenvals(i), ritzVector(:,i), 1, residual(:,i), 1)
    end do

    if (verbose) then
      print *
      print *, 'Residual'
      call printMatrix(residual, eigen_in, ndimA)
      print *
    end if


  ! compute norm and check convergency
    do i = 1, eigen_in
      print *, 'residual norm of eigenvector', i, ':', dnrm2(ndimA, residual(:,i), 1)
      print *
      if(dnrm2(ndimA, residual, 1) <= threshold_residual) then
        print *, 'Converged for:', i
        converged(i) = .true.
      end if
      if (all(converged)) then
        print *, 'Converged for all sought eigenpairs.'
        exit outer
      end if
    end do
    


  ! check if matrix is not full
    matrix_not_full = .true.
    do i = 1, ndimA
      if (matV(i, ndimV + 1 - eigen_in) /= 0.0d0) then
        matrix_not_full = .false.   
        exit
      end if
    end do


    if (matrix_not_full) then
    ! precondition y_i = r_i / (D - lambda_i)
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

      
  ! for eigen_in > 1 orthogonalize residual matrix

    if (eigen_in .GT. 1) then

      call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, lw, -1, info)
      lwork = int(lw(1))
      allocate(work(lwork))
      call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, work, lwork, info)
      call checkInfo(info, 'Orthogonalization of residual step 1')


      call dorgqr(ndimA, eigen_in, eigen_in, residual, ndimA, tau, work, lwork, info)
      call checkInfo(info, 'Orthogonalization of residual step 2')
      deallocate(work)

      call checkOrth1mat(residual, eigen_in, 'Residual', thresh_GS)
      do i = 1, eigen_in
        do j = 1, eigen_in
          if (j .NE. i) then
            check_GS = abs(dot_product( residual(:,j), residual(:,i) ))
          end if
          if ( check_GS .GT. thresh_GS ) then
            print *
            print *, '--- WARNING ---'
            print *, 'residual vector of index', i, 'is not orthogonal to residual vector', j, 'in matrix V'
            print *, 'Result of dot product:', check_GS
            print *
          end if
        end do
      end do
    end if


  ! Gram Schmidt orthogonalization

      if (GS_in_loop) then

  ! loop implementation
        do i = 1, eigen_in
          temp = 0.0d0
          do j = 1, n_grow
            temp(:,i) = temp(:,i) + dot_product(residual(:,i), matV(:,j)) / dot_product(matV(:,j), matV(:,j)) * matV(:,j)
          end do
          residual(:,i) = residual(:,i) - temp(:,i)  
          residual(:,i) = residual(:,i) / dnrm2(ndimA, residual(:,i), 1)
        end do 

      ! print a warning if matrix is not orthogonal
          call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS)

      else

  ! matrix implementation

      ! matrix product V^T * y  
        call dgemm('t', 'n',  n_grow, eigen_in ,ndimA, 1.0d0, matV, ndimA, residual, ndimA, 0.0d0, temp_mat, ndimV)

      ! matrix product V * (V^T * y) 
        call dgemm('n', 'n', ndimA, eigen_in, n_grow, 1.0d0, matV, ndimA, temp_mat, ndimV, 0.0d0, temp_mat_prime, ndimA)

      ! y_prime = y - V * (V^T * y) 
        do i = 1, eigen_in
          residual(:, i) = residual(:, i) - temp_mat_prime(:, i)
        end do


        if (verbose) then
          print *, 'y_i: Matrix product V^T * y'
          call printMatrix(temp_mat, eigen_in, n_grow)
          print *
          print *, 'y_i: Matrix product V * ( V^T * y )'
          call printMatrix(temp_mat_prime, eigen_in, ndimA)
          print *
          print *, 'Orthogonalized precondition y'
          call printMatrix(residual, eigen_in, ndimA)
        end if


    ! print a warning if matrix is not orthogonal
        call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS)
          

    !  orthonormalization
        do i = 1, eigen_in
          residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
        end do
        print *
        print *, 'Orthonormalized precondition y'
        call printMatrix(residual, eigen_in, ndimA)

        print *

      end if  

      ! add orthonormal residual to subspace

      do i = 1, eigen_in
        matV(:, n_grow  + i) = residual(:,i)
      end do

      if (verbose) then
        print *
        print *, 'New subspace V'
        call printMatrix(matV, n_grow + eigen_in, ndimA)
        print *
        print *
      end if 

    ! restart 
      n_grow = n_grow + eigen_in 


    else
    ! treat case if matrix full
      print *, 'Matrix V full'
      matV            = zero

      do i = 1, eigen_in
        matV(:, i) = ritzVector(:, i)
      end do

      print *, 'New subspace'
      call printMatrix(matV, eigen_in, ndimA)
      print *
      print *

    ! restart

      matW            = zero
      matP            = zero
      eigenvals       = zero
      eigenvecs       = zero
      ritzVector      = zero
      temp_mat        = zero
      temp_mat_prime  = zero
      residual        = zero
      temp            = zero

      n_grow          = eigen_in
    
    end if



  end do outer
  print *, 'test'
  deallocate(matA, mask, diagonalA, matV, matW, matP, eigenvals, eigenvecs, &
             ritzVector, residual, ritzVectorTemp, temp_mat, temp_mat_prime, temp, converged)
  print *, 'test'

  end subroutine symmetricDavidson


  subroutine lapackDiag(mat, eigenvals, ndimMat)
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


!   print formatted matrix
  subroutine printMatrix(mat, nrows, ncols) 

    real(wp), intent(in)  :: mat(:,:)
    integer , intent(in)  :: nrows, ncols
    integer :: i,j 

    do i = 1, ncols
      print *, (mat(i,j), j= 1, nrows )
    end do

  end subroutine printMatrix


    subroutine printVector(vec, lenRow)

      real(wp), intent(in)  :: vec(:)
      integer,  intent(in)  :: lenRow
      integer               :: i

      do i = 1, lenRow
        print *, vec(i)
      end do

    end subroutine printVector


!   check if info is zero and print error message if not
    subroutine checkInfo(info, occasion)
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



!   check if vectors of two matrices are orthogonal
    subroutine checkOrth2mat(mat1, nrows1, mat1_name, mat2, nrows2, mat2_name, thresh)
      integer,              intent(in)      ::  nrows1, nrows2
      real(wp),             intent(in)      ::  thresh, mat1(:,:), mat2(:,:)
      character(len=*),    intent(in)      ::  mat1_name, mat2_name
      real(wp)                              ::  dot_prod
      integer :: i,j 

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


!   check if the vectors within the matrix are orthogonal
    subroutine checkOrth1mat(mat1, nrows1, mat1_name, thresh)

      integer,              intent(in)      ::  nrows1
      real(wp),             intent(in)      ::  thresh, mat1(:,:)
      character(len=*),    intent(in)       ::  mat1_name
      real(wp)                              ::  dot_prod
      integer :: i,j 

      do i = 1, nrows1
        do j = 1, nrows1
          if (j .NE. i) then
            dot_prod  = abs(dot_product( mat1(:,j), mat1(:,i) ))
          end if
          if ( dot_prod .GT. thresh ) then
            print *
            print *, '--- WARNING ---'
            print *, 'vector', i, 'of matrix', mat1_name, 'is not orthogonal to vector', j, 'of matrix', mat1_name
            print *, 'Result of dot product:', dot_prod
            print *
          end if
        end do
      end do

    end subroutine checkOrth1mat
 
  end program davidson
      
