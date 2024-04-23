program davidson
    implicit none
    
    real, allocatable :: matA(:,:), matV(:,:) 
    integer :: i, j, it, ndimA, ndimV, maxiter, idxMaxVal(2)
    
    ndimA = 5
    ndimV = 20
    maxiter = 10

!get matrix
    allocate(matA( ndimA, ndimA ))

    do i = 1, ndimA 
        do j = 1, ndimA
            if (j == i) then
                matA(i,j) = i + 1
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

    do it = 1, maxiter
        
    end do
    
!output
    write(*,*) 'Input Matrix A:'
    call printMatrix(matA)



contains
    subroutine printMatrix(mat) 
        integer :: i
        real, allocatable :: mat(:,:)
        do i = 1, size(mat(:,i))
            print *, mat(i, :)
        end do
    end subroutine printMatrix
    

end program davidson
    
