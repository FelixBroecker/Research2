program davidson
  implicit none
 
  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp)               :: lw(1)
  real(wp), allocatable  :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
  real(wp), allocatable  :: ritzVector(:,:), temp(:,:), ritzVectorTemp(:)
  real(wp), allocatable  :: residual(:,:)
  integer                :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(1), lwork, info, eigen_in
  real(wp)               :: dnrm2

  ndimA = 5
  ndimV = 6
  maxiter = 10
  eigen_in = 1

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
  allocate(temp(ndimA, ndimA))
  allocate(ritzVectorTemp(ndimA))
  allocate(residual(ndimA, ndimA))

  do it = 1, 1
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

    write(*,*) 'Eigenvalues:', eigenvals
    write(*,*) 'Matrix P'
    call printMatrix(matP)
    write(6,*) 'eigenvecs'
    call printMatrix(eigenvecs)
      
  ! get ritz vector with to matrices
    write(*,*) 'Matrix V'
    call printMatrix(matV)
    call dgemm('n', 'n', ndimA, it * eigen_in, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)
    write(*,*) 'Ritzvector all'
    call printMatrix(ritzVector)

  ! get residual
    ! go column wise and do residual = matW_i * eigenvectors_i * matV_i


    !print *, maxloc(eigenvals)
    ! dot(eigenvals,  ritzvector)
    !call dgemv('n', ndimA, ndimV, 1.0d0, ritzVector, ndimA, eigenvals, 1, 0.0d0, temp, 1)
    !call dgemm('n', 'n', ndimA, ndimV, ndimV, -1.0d0, matW, ndimA, eigenvals, ndimV, 1.0d0, temp, ndimA)
    !print *, 'temp mat'
    !call printMatrix(temp)
  
    do i = 1, eigen_in
      residual(:,i) = matW(:,i) - eigenvals( it * i) * matV(:,i)
    enddo
    print *, 'residual', residual
    !do i = 1, n_search
    !  res(:,i) = w(:,i) - lambda(it)*v(:,i)
    !enddo

    ! compute norm
! if( conv)

! elseif (matV not full)
    ! precondition

    !do i = i, n_add (n_add = it*n_search)
    !  res(:,i) = res(:,i) / (diag(i) - lambda(i))
    !enddo

    do i = 0, it - 1

     ! !compute ritzvector for index i
     !   call dgemv('n', ndimA, ndimV, 1.0d0, matV, ndimA, eigenvecs(:, size( eigenvecs(1,:) ) - i), 1, 0.0d0, ritzVectorTemp, 1)
     !   print *, 'ritz Vector index i'
     !   print *, ritzVectorTemp
     ! !compute A - eigenval of index i * Identity
     !   temp = matA - eigenvals(size(eigenvals) - i) * identity
     !   print *, 'temp mat'
     !   call printMatrix(temp)
     ! !dot Product of both
     !   call dgemv('n', ndimA, ndimA, 1.0d0, temp, ndimA, ritzVectorTemp, 1, 0.0d0, residual, 1)
     !   print *, 'Vector r_i', residual
     ! !get new directions tk
     !   direction = eigenvals(size(eigenvals) - i) * identityVec - diagonalA
     !   direction = residual / direction 
     !   print *, 'direction', direction
        
 !else

   !restart

 !endif
    enddo
    
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
      
