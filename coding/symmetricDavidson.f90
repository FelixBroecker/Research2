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
  logical                :: matrix_not_full

  ndimA = 5
  ndimV = 6
  maxiter = 2
  eigen_in = 1
  threshold_residual = 0.001d0

! get matrix
  allocate(matA( ndimA, ndimA ))
  allocate(diagonalA(ndimA))

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
  write(*,*) 'Input Matrix A:'
  call printMatrix(matA)

! get initial vector
  allocate(matV( ndimA, ndimV ))
  ! search for colum with maximum value and use that column for initial vector
  idxMaxVal = maxloc(diagonalA)
  matV(idxMaxVal, 1) = 1
  write(*,*) 'Initital Vector'
  call printMatrix(matV)

! get projection of A
  allocate(matW(ndimA, ndimV))
  allocate(matP(ndimV, ndimV))
  allocate(eigenvals(ndimV))
  allocate(eigenvecs(ndimV, ndimV))
  allocate(ritzVector(ndimA,ndimV))
  allocate(temp_mat(ndimV, eigen_in))
  allocate(temp_mat_prime(ndimA, eigen_in))
  allocate(ritzVectorTemp(ndimA))
  allocate(residual(ndimA, eigen_in))

idx = 0 
  do it = 1, maxiter
    call dgemm('n', 'n', ndimA, it * eigen_in, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
    write(*,*) 'Matrixproduct A V :'
    call printMatrix(matW)

    call dgemm('t', 'n', ndimV, it * eigen_in, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)
    write(*,*) 'Matrixproduct (A V)T V'
    call printMatrix(matP)

  ! diagonalize and obtain eigenvalues and eigenvectors 

    ! store matP in eigenvecs because dsyev writes eigenvecs in input matrix
    eigenvecs = matP

    call dsyev('V', 'U', it, eigenvecs, it, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'U', it, eigenvecs, it, eigenvals, work, lwork, info)
    deallocate(work)
    print *, 'info', info

    write(*,*) 'Eigenvalues:', eigenvals
    write(*,*) 'Matrix P'
    call printMatrix(matP)
    write(6,*) 'eigenvecs'
    call printMatrix(eigenvecs)
      
  ! get ritz vector 
    write(*,*) 'Matrix V'
    call printMatrix(matV)
    call dgemm('n', 'n', ndimA, it * eigen_in, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)
    write(*,*) 'Ritzvector all'
    call printMatrix(ritzVector)

  ! get residual
    ! go column wise and do residual = matW_i * eigenvectors_i * matV_i

    do i = 1, eigen_in
      residual(:,i) = matW(:, idx + i) - eigenvals(idx + i) * matV(:, idx + i)
    end do

    print *, 'residual'
    call printMatrix(residual)

    ! compute norm
    print *, 'residual norm', dnrm2(ndimA, residual, 1)
    if(dnrm2(ndimA, residual, 1) <= threshold_residual) then
      print *, 'jippi, converged'
      exit 
    end if
    
    ! check if matrix is not full
    matrix_not_full = .true.
    do i = 1, ndimA
      if (matV(i, ndimV) /= 0.0d0) then
        matrix_not_full = .false.
        exit
      end if
    end do
    print *, 'matrix not full', matrix_not_full

    if (matrix_not_full) then

      ! precondition 
      do i = 1, eigen_in
        residual(:,i) = residual(:,i) / ( diagonalA(i) - eigenvals(i))
      end do

      print *, 'preconditioned residual'
      call printMatrix(residual)
      
    ! Gram Schmit orthogonalization

      ! matrix product V^T * y  
      call dgemm('t', 'n', ndimA, eigen_in * it, ndimA, 1.0d0, matV, ndimV, residual, ndimA, 0.0d0, temp_mat, ndimV)
      write(*,*) 'matrix temp_mat'
      call printMatrix(temp_mat)

      ! matrix product V * (V^T * y) 
      call dgemm('n', 'n', ndimA, eigen_in * it, ndimV, 1.0d0, matV, ndimV, temp_mat, ndimV, 0.0d0, temp_mat_prime, ndimA)
      write(*,*) 'matrix temp_mat_prime'
      call printMatrix(temp_mat_prime)

      ! y_prime = y - temp_mat_prime
      do i = 1, eigen_in
        residual(:, i) = residual(:, i) - temp_mat_prime(:, i)
      end do
      write(*,*) 'orthogonalized precondition y'
      call printMatrix(residual)

      ! orthonormalization
        do i = 1, eigen_in
          residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
        end do
      write(*,*) 'orthonormalized precondition y'
      call printMatrix(residual)
      ! add to V subspace
      do i = 1, eigen_in
        matV(:, idx + 1 + i) = residual(:,i)
      end do
      write(*,*) 'new subspace V'
      call printMatrix(matV)
      
    else
      print *, 'Matrix V full'
      deallocate(matV)
      allocate(matV(ndimA, ndimV))
      do i = 1, eigen_in
        matV(:, i) = ritzVector(:, idx + i)
      end do

    end if
    
  idx = idx + eigen_in
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
      
