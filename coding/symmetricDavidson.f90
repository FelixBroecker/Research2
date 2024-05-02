program davidson
  implicit none
 
  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp)               :: lw(1), threshold_residual
  real(wp), allocatable  :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
  real(wp), allocatable  :: ritzVector(:,:), temp_mat(:,:), ritzVectorTemp(:)
  real(wp), allocatable  :: residual(:,:), temp_mat_prime(:,:)
  integer                :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(1), lwork, info, eigen_in, idx
  real(wp)               :: dnrm2
  logical                :: matrix_not_full, verbose


  ndimA                 = 5
  ndimV                 = 6
  maxiter               = 2
  eigen_in              = 1
  threshold_residual    = 0.001d0
  verbose               = .true.


  allocate(matA( ndimA, ndimA ))
  allocate(diagonalA(ndimA))
  allocate(matW(ndimA, ndimV))
  allocate(matP(ndimV, ndimV))
  allocate(eigenvals(ndimV))
  allocate(eigenvecs(ndimV, ndimV))
  allocate(ritzVector(ndimA,ndimV))
  allocate(temp_mat(ndimV, eigen_in))
  allocate(temp_mat_prime(ndimA, eigen_in))
  allocate(ritzVectorTemp(ndimA))
  allocate(residual(ndimA, eigen_in))


! get matrix
  do i = 1, ndimA 
    do j = i, ndimA
      if (j == i) then
        matA(i,j) = i + 1
        diagonalA(i) = matA(i,j)
      else if (j > i) then
        matA(i,j) = 1 / (dble(i + j))
        matA(j,i) = matA(i,j)
      end if
    end do
  end do


! get initial vector
  allocate(matV( ndimA, ndimV ))
  ! search for colum with maximum value and use that column for initial vector
  idxMaxVal = maxloc(diagonalA)
  matV(idxMaxVal, 1) = 1
  write(*,*) 'Initital Vector'
  call printMatrix(matV)


! start loop
  idx = 0 
  do it = 1, maxiter
      print *, '------------------------------'
      print *, 'Iteration', it
      print *, '------------------------------'

  ! get projection matrix P = (W)^T * V;  W = A * V
    call dgemm('n', 'n', ndimA, it * eigen_in, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
    call dgemm('t', 'n', ndimV, it * eigen_in, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)


  ! diagonalize and obtain eigenvalues and eigenvectors 
    eigenvecs = matP
    call dsyev('V', 'U', it, eigenvecs, it, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'U', it, eigenvecs, it, eigenvals, work, lwork, info)
    deallocate(work)


  ! get ritz vector 
    call dgemm('n', 'n', ndimA, it * eigen_in, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)


  ! get residual
  ! go column wise and do residual = matW_i * eigenvectors_i * matV_i
    do i = 1, eigen_in
      residual(:,i) = matW(:, idx + i) - eigenvals(idx + i) * matV(:, idx + i)
    end do

    print *
    print *, 'Residual'
    call printMatrix(residual)


  ! compute norm and check convergency
    print *, 'residual norm', dnrm2(ndimA, residual, 1)
    if(dnrm2(ndimA, residual, 1) <= threshold_residual) then
      print *, 'jippi, converged'
      exit 
    end if
    

  ! check if matrix is not full
    matrix_not_full = .true.
    do i = 1, ndimA
      if (matV(i, ndimV) /= 0.0d0) then
        print *, 'Matrix is full'
        matrix_not_full = .false.   
        exit
      end if
    end do


    if (matrix_not_full) then
    ! precondition y_i = r_i / (D - lambda_i)
      do i = 1, eigen_in
        residual(:,i) = residual(:,i) / (diagonalA(i) - eigenvals(i))
      end do

      print *
      print *, 'Preconditioned residual'
      call printMatrix(residual)

      
    ! Gram Schmidt orthogonalization
    ! matrix product V^T * y  
      call dgemm('t', 'n', ndimA, eigen_in * it, ndimA, 1.0d0, matV, ndimV, residual, ndimA, 0.0d0, temp_mat, ndimV)

    ! matrix product V * (V^T * y) 
      call dgemm('n', 'n', ndimA, eigen_in * it, ndimV, 1.0d0, matV, ndimV, temp_mat, ndimV, 0.0d0, temp_mat_prime, ndimA)

    ! y_prime = y - V * (V^T * y) 
      do i = 1, eigen_in
        residual(:, i) = residual(:, i) - temp_mat_prime(:, i)
      end do
      print *
      print *, 'Orthogonalized precondition y'
      call printMatrix(residual)

    ! orthonormalization
      do i = 1, eigen_in
        residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
      end do
      print *
      print *, 'Orthonormalized precondition y'
      call printMatrix(residual)


    ! add to V subspace
      do i = 1, eigen_in
        matV(:, idx + 1 + i) = residual(:,i)
      end do
      

    else
      print *, 'Matrix V full'
      deallocate(matV)
      allocate(matV(ndimA, ndimV))

      do i = 1, eigen_in
        matV(:, i) = ritzVector(:, idx + i)
      end do

    end if

    idx = idx + eigen_in


! generate output
    if (verbose) then
      print *
      print *, 'Input Matrix A:'
      call printMatrix(matA)
      print *
      print *, 'Matrix V'
      call printMatrix(matV)
      print *,
      print *, 'Matrixproduct W (A V) :'
      call printMatrix(matW)
      print *
      write(*,*) 'Matrixproduct (W)T V = P'
      call printMatrix(matP)
      print *, 'Info about diagonalization (0 = worked)', info
      print *
      print *, 'Eigenvalues:', eigenvals
      print *, 'Eigenvectors:'
      call printMatrix(eigenvecs)
      print *
      print *, 'Ritzvector all'
      call printMatrix(ritzVector)
      print *
      print *, 'y_i: Matrix product V^T * y'
      call printMatrix(temp_mat)
      print *
      print *, 'y_i: Matrix product V * ( V^T * y )'
      call printMatrix(temp_mat_prime)
      print *
      print *, 'New subspace V'
      call printMatrix(matV)
      print *
      print *

    end if

  end do


  contains
    subroutine printMatrix(mat) 
        integer :: i, length
        real(wp), allocatable :: mat(:,:)
        length = size(mat(:,i))
        do i = 1, length
          print *, mat(i, :)
        end do
    end subroutine printMatrix
      

  end program davidson
      
