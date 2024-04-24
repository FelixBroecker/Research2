program davidson
  implicit none
 
  intrinsic :: selected_real_kind
  integer, parameter :: wp = selected_real_kind(15)
  real(wp), allocatable :: matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), work(:)
  integer :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(2), lwork, info
  
  ndimA = 5
  ndimV = 20
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
  lwork = max(1, 3 * ndimV - 1)
  allocate(work(lwork))
  call dsyev('V', 'U', ndimV, matP, ndimV, eigenvals, work, lwork, info)
  write(*,*) 'Eigenvalues:', eigenvals
  
!output
  write(*,*) 'Input Matrix A:'
  call printMatrix(matA)



contains
  subroutine printMatrix(mat) 
      integer :: i
      real(wp), allocatable :: mat(:,:)
      do i = 1, size(mat(:,i))
        print *, mat(i, :)
      end do
  end subroutine printMatrix
    

end program davidson
    
