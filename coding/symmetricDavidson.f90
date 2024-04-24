program davidson
  implicit none
 
  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp)               :: lw(1)
  real(wp), allocatable  :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), work(:) 
  integer                :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(2), lwork, info
  real(wp), allocatable  :: ritzVector(:,:)
  real(wp), allocatable  :: ritzVector2(:,:)

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
  end do

!diagonalize 
  allocate(eigenvals(ndimV))
  call dsyev('V', 'U', ndimV, matP, ndimV, eigenvals, lw, -1, info)
  lwork = int(lw(1))
  allocate(work(lwork))
  !lwork = max(1, 3 * ndimV - 1)
  call dsyev('V', 'U', ndimV, matP, ndimV, eigenvals, work, lwork, info)
  write(*,*) 'Eigenvalues:', eigenvals
  write(*,*) 'Matrix P'
  call printMatrix(matP)
  
  
!get ritz vector
  write(*,*) 'Matrix V'
  call printMatrix(matV)
  write(6,*) 'before'
  allocate(ritzVector2(ndimA,ndimV))
!  allocate(ritzVector(ndimA,ndimV))
  !call dgemv('n', ndimA, ndimV, 1.0d0, matV, ndimA, eigenvals, 1, 0.0d0, ritzVector, 1)
  !write(*,*) 'Ritzvector'
  !call printMatrix(ritzVector)
  write(6,*) 'matrix', ritzVector 

  contains
    subroutine printMatrix(mat) 
        integer :: i
        real(wp), allocatable :: mat(:,:)
        do i = 1, size(mat(:,i))
          print *, mat(i, :)
        end do
    end subroutine printMatrix
      

  end program davidson
      
