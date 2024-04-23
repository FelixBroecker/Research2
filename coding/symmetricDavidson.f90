program davidson
    implicit none
    
    real, allocatable :: matA(:,:), matV(:,:) 
    integer :: i, j, ndimA, ndimV
    
!get matrix
    ndimA = 5
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

    allocate(matV( ndimV, ndimA ))
    
    write(*,*) 'Input Matrix A:'
    call printMatrix(matA)



contains
    subroutine printMatrix(mat) 
        integer :: i
        real, allocatable :: mat(:,:)
        do i = 1, size(mat(1,:))
            print *, mat(i, :)
        end do
    end subroutine printMatrix
    

end program davidson
    
