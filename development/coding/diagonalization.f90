module diagonalization
  implicit none
  
contains
  subroutine cholesky(m, n, A)
  intrinsic                 ::  selected_real_kind, sqrt
  integer,  parameter       ::  wp = selected_real_kind(15)
!
!  takes matrix A with m rows and n columns
!
    integer,    intent(in)    :: m, n
    real(wp),   intent(inout) :: A(:,:)
    real(wp),   allocatable   :: L(:,:)
    integer                   :: i, k, j 
    real(wp)                  :: sum1, sum2

    allocate(L(m,n))
    L = 0.d0
  
    do i=1, n
      do k=1, i
        sum1 = 0.d0
        sum2 = 0.d0
        do j=1, k-1 
          if (i .eq. k) then
            sum1 = sum1 + L(k,j) * L(k,j)
          else if (i .gt. k) then
            sum2 = sum2 + L(i, j) * L(k, j)
          else 
            L(i,k) = 0.d0
          end if
        end do
!
        if (k .eq. i) then
          L(k,k) = sqrt(A(k,k) - sum1)
        else if (i .gt. k) then
          L(i,k) = (1.d0 / L(k,k)) * (A(i,k) - sum2)
        else 
          L(i,k) = 0.d0
        end if
      end do
    end do

  do i=1, m
    print "(3(1X,F6.2))", L(i,:)
  end do

  end subroutine cholesky

end module diagonalization
