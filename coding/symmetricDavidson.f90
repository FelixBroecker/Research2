program davidson
  implicit none
 
  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp)               :: lw(1)
  real(wp), allocatable  :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
  real(wp), allocatable  :: ritzVector(:,:), identity(:,:), temp(:,:), ritzVectorTemp(:), residual(:)
  integer                :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(2), lwork, info

  ndimA = 5
  ndimV = 6
  maxiter = 10

!get matrix
  allocate(matA( ndimA, ndimA ))
  allocate(diagonalA(ndimA))

  do i = 1, ndimA 
    do j = 1, ndimA
      if (j == i) then
        matA(i,j) = i + 1
        diagonalA(i) = matA(i,j)
      else if (j > i) then
        matA(i,j) = 1 / (real(i + j))
        matA(j,i) = matA(i,j)
      end if
    end do
  end do
  write(*,*) 'Input Matrix A:'
  call printMatrix(matA)

!get initial vector
  allocate(matV( ndimA, ndimV ))
  !search for colum with maximum value and use that column for initial vector
  idxMaxVal = maxloc(matA)
  matV(idxMaxVal(1), 1) = 1
  write(*,*) 'Initital Vector'
  call printMatrix(matV)

!get projection of A
  allocate(matW(ndimA, ndimV))
  allocate(matP(ndimV, ndimV))

  do it = 1, 1
    call dgemm('n', 'n', ndimA, ndimV, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
    write(*,*) 'Matrixproduct A V :'
    call printMatrix(matW)

    call dgemm('t', 'n', ndimV, ndimV, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)
    write(*,*) 'Matrixproduct (A V)T V'
    call printMatrix(matP)

  !diagonalize and obtain eigenvalues and eigenvectors 
    allocate(eigenvals(ndimV))
    allocate(eigenvecs(ndimV, ndimV))

    !store matP in eigenvecs because dsyev writes eigenvecs in input matrix
    eigenvecs = matP

    call dsyev('V', 'U', ndimV, eigenvecs, ndimV, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'U', ndimV, eigenvecs, ndimV, eigenvals, work, lwork, info)

    write(*,*) 'Eigenvalues:', eigenvals
    write(*,*) 'Matrix P'
    call printMatrix(matP)
    write(6,*) 'eigenvecs'
    call printMatrix(eigenvecs)
      
  !get ritz vector with to matrices
    write(*,*) 'Matrix V'
    call printMatrix(matV)
    allocate(ritzVector(ndimA,ndimV))
    call dgemm('n', 'n', ndimA, ndimV, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)
    write(*,*) 'Ritzvector'
    call printMatrix(ritzVector)
  
  !get residual
    !get identity matrix and put eigenvalues on diagonal and substract from matA
    allocate(identity(ndimA, ndimA))
    allocate(temp(ndimA, ndimA))

    !print *, maxloc(eigenvals)
    ! dot(eigenvals,  ritzvector)
    !call dgemv('n', ndimA, ndimV, 1.0d0, ritzVector, ndimA, eigenvals, 1, 0.0d0, temp, 1)
    !call dgemm('n', 'n', ndimA, ndimV, ndimV, -1.0d0, matW, ndimA, eigenvals, ndimV, 1.0d0, temp, ndimA)
    !print *, 'temp mat'
    !call printMatrix(temp)

    identity(1:ndimA, 1:ndimA) = 0.0d0
    forall(i = 1:ndimA) identity(i, i) = 1.0d0
  
    do i = 0, it - 1
      !compute ritzvector for index i
        allocate(ritzVectorTemp(ndimA))
        allocate(residual(ndimA))
        call dgemv('n', ndimA, ndimV, 1.0d0, matV, ndimA, eigenvecs(:, size( eigenvecs(1,:) ) - i), 1, 0.0d0, ritzVectorTemp, 1)
        print *, 'ritz Vector Temp'
        print *, ritzVectorTemp
      !compute A - eigenval of index i * Identity
        temp = matA - eigenvals(size(eigenvals) - i) * identity
        print *, 'temp mat'
        call printMatrix(temp)
      !dot Product of both
        call dgemv('n', ndimA, ndimA, 1.0d0, temp, ndimA, ritzVectorTemp, 1, 0.0d0, residual, 1)
        print *, 'Vector r_i', residual


    !call printMatrix(identity) 
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
      
