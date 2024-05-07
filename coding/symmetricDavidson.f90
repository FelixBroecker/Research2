program davidson
  implicit none
 
  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp)               :: lw(1), threshold_residual, zero, check_GS, thresh_GS
  real(wp), allocatable  :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
  real(wp), allocatable  :: ritzVector(:,:), temp_mat(:,:), ritzVectorTemp(:)
  real(wp), allocatable  :: residual(:,:), temp_mat_prime(:,:), diff, temp(:,:)
  integer                :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(1), lwork, info, eigen_in, idx, n_grow
  real(wp)               :: dnrm2
  logical, allocatable   :: mask(:), converged(:)
  logical                :: matrix_not_full, verbose, GS_in_loop


  ndimA                 = 5
  ndimV                 = 6
  maxiter               = 6
  eigen_in              = 1
  threshold_residual    = 1.d-4
  verbose               = .false.
  GS_in_loop            = .false.
  thresh_GS             = 1.d-6

  zero                  = 0.0d0


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
  allocate(temp(ndimA, eigen_in))
  allocate(mask(ndimA))
  allocate(matV( ndimA, ndimV ))
  allocate(converged(eigen_in))



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
  call printMatrix(matA)
  print *


! get initial vector
  ! search for colum with maximum value and use that column for initial vector

  do i = 1, eigen_in
    idxMaxVal = minloc(diagonalA, mask=mask)
    matV(idxMaxVal, i) = 1.0d0
    mask(idxMaxVal) = .false.
  end do 

  write(*,*) 'Initial Vector'
  call printMatrix(matV)
  print *


  idx = 1
  n_grow = eigen_in

! start loop
  outer: do it = 1, maxiter
      print *, '------------------------------'
      print *, 'Iteration', it
      print *, '------------------------------'

  eigenvecs = 0.0d0
  eigenvals = 0.0d0
  ! get projection matrix P = (W)^T * V;  W = A * V
    call dgemm('n', 'n', ndimA, n_grow, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
    call dgemm('t', 'n', ndimV, n_grow, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)


    if (verbose) then
      print *, 'Matrixproduct W (A V) :'
      call printMatrix(matW)
      print *
      write(*,*) 'Matrixproduct (W)T V = P'
      call printMatrix(matP)
    end if



  ! diagonalize and obtain eigenvalues and eigenvectors 
    eigenvecs = matP
    call dsyev('V', 'u', idx, eigenvecs, ndimV, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'u', idx, eigenvecs, ndimV, eigenvals, work, lwork, info)
    deallocate(work)

    if (verbose) then
      print *, 'Info about diagonalization (0 = worked)', info
    end if

      print *
      print *, 'Eigenvalues:', eigenvals
      print *, 'Eigenvectors:'
      call printMatrix(eigenvecs)


  ! get ritz vector 
    call dgemm('n', 'n', ndimA, idx, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)

    print *
    print *, 'Ritzvector all'
    call printMatrix(ritzVector)
    print *


  ! get residual
    call dgemm('n', 'n', ndimA, n_grow, ndimV, 1.0d0, matW, ndimA, eigenvecs, ndimV, 0.0d0, residual, ndimA)

    do i = 1, eigen_in
      call daxpy(ndimA, -eigenvals(i), ritzVector(:,i), 1, residual(:,i), 1)
    end do

    if (verbose) then
      print *
      print *, 'Residual'
      call printMatrix(residual)
    end if


  ! compute norm and check convergency
  do i = 1, eigen_in
    print *, 'residual norm of eigenvector', i, ':', dnrm2(ndimA, residual(:,i), 1)
    if(dnrm2(ndimA, residual, 1) <= threshold_residual) then
      print *, 'Converged for:', i
      converged(i) = .true.
    end if
    if (all(converged)) then
      print *, 'Converged for all sought eigenvalues.'
      exit outer
    end if
  end do
    

  ! check if matrix is not full
    matrix_not_full = .true.
    do i = 1, ndimA
      if (matV(i, ndimV) /= 0.0d0) then
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


  ! Gram Schmit orthogonalization

      if (GS_in_loop) then

  ! loop implementation
        do i = 1, eigen_in
          temp = 0.0d0
          do j = 1, idx
            temp(:,i) = temp(:,i) + dot_product(residual(:,i), matV(:,j)) / dot_product(matV(:,j), matV(:,j)) * matV(:,j)
          end do
          residual(:,i) = residual(:,i) - temp(:,i)  
          residual(:,i) = residual(:,i) / dnrm2(ndimA, residual(:,i), 1)
        end do 

      ! check if orthonormal
        do i = 1, eigen_in
          do j = 1, idx
            check_GS = abs(dot_product( matV(:,j), residual(:,i) ))
            if ( check_GS .GT. thresh_GS ) then
              print *
              print *, '--- WARNING ---'
              print *, 'precondition of index', i, 'is not orthogonal to vector', j, 'in matrix V'
              print *, 'Result of dot product:', check_GS
              print *
            end if
          end do
        end do

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
          call printMatrix(temp_mat)
          print *
          print *, 'y_i: Matrix product V * ( V^T * y )'
          call printMatrix(temp_mat_prime)
          print *
          print *, 'Orthogonalized precondition y'
          call printMatrix(residual)
        end if


        do i = 1, eigen_in
          do j = 1, idx
            check_GS = abs(dot_product( matV(:,j), residual(:,i) ))
            if ( check_GS .GT. thresh_GS ) then
              print *
              print *, '--- WARNING ---'
              print *, 'precondition of index', i, 'is not orthogonal to vector', j, 'in matrix V'
              print *, 'Result of dot product:', check_GS
              print *
            end if
          end do
        end do
          

    !  orthonormalization
        do i = 1, eigen_in
          residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
        end do
        print *
        print *, 'Orthonormalized precondition y'
        call printMatrix(residual)

        print *

      end if  

      ! add orthonormal residual to subspace

      do i = 1, eigen_in
        matV(:, idx  + i) = residual(:,i)
      end do

      if (verbose) then
        print *
        print *, 'New subspace V'
        call printMatrix(matV)
        print *
        print *
      end if 

    ! restart 
      idx = idx + eigen_in
      n_grow = n_grow + eigen_in 


    else
    ! treat case if matrix full
      print *, 'Matrix V full'
      matV            = zero

      do i = 1, eigen_in
        matV(:, i) = ritzVector(:, idx - 1 + i)
      end do

      print *, 'New subspace'
      call printMatrix(matV)
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

      idx             = 1
      n_grow          = eigen_in
    
    end if



  end do outer


  contains
    subroutine printMatrix(mat) 
        integer :: i, length
        real(wp), intent(in) :: mat(:,:)
        length = size(mat(:,1))
        do i = 1, length
          print *, mat(i, :)
        end do
    end subroutine printMatrix
      

  end program davidson
      
