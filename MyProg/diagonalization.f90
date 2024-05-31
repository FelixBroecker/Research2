module diagonalization
  implicit none
  
contains
  subroutine symmetricDavidson(mat_in, dim_mat_in, return_eigenvecs, return_eigenvals, eigen_in, verbose, &
      threshold_residual, maxiter)

    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)

    real(wp), intent(in)      ::  mat_in(:,:), threshold_residual
    real(wp), intent(out)     ::  return_eigenvecs(:,:), return_eigenvals(:)
    integer, intent(in)       ::  dim_mat_in, eigen_in, verbose, maxiter
    real(wp)                  ::  lw(1), zero, check_GS, thresh_GS
    real(wp), allocatable     ::  matA(:,:), matV(:,:), matW(:,:), matP(:,:),  diagonalA(:), eigenvals(:), eigenvecs(:,:), work(:)
    real(wp), allocatable     ::  ritzVector(:,:), temp_mat(:,:), tau(:)
    real(wp), allocatable     ::  residual(:,:), temp_mat_prime(:,:),  temp(:,:)
    integer                   ::  i, j, it, ndimA, ndimV, idxMaxVal(1), lwork, info,  n_grow,  max_orth
    real(wp)                  ::  dnrm2, diff
    logical, allocatable      ::  mask(:), converged(:)
    logical                   ::  matrix_not_full, GS_in_loop, not_orthogonal, size_not_exceeded



    ndimA                 = dim_mat_in
    ndimV                 = 20 * eigen_in
    max_orth              = 10
    thresh_GS             = 1.d-12
    GS_in_loop            = .false.

!   allocate space to create matrix A
    allocate(matA( ndimA, ndimA), mask(ndimA), diagonalA(ndimA))

!   allocate space to obtain reduced space
    allocate(matV( ndimA, ndimV ), matW(ndimA, ndimV), matP(ndimV, ndimV))

!   allocate space for diagonalization
    allocate(eigenvals(ndimV), eigenvecs(ndimV, ndimV))

!   allocate space for ritz vector and residual
    allocate(ritzVector(ndimA,ndimV), residual(ndimA, eigen_in))

!   allocate space for diagonalization
    allocate(temp_mat(ndimV, eigen_in), temp_mat_prime(ndimA, eigen_in), temp(ndimA, eigen_in))

!   allocate space for else
    allocate(converged(eigen_in))



    zero            = 0.0d0
    
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


    if (verbose .ge. 1) then
      print *
      print *, '-------------------- start of Davidson routine -----------------'
      print *
    end if

!   get matrix diagonal
    do i = 1, dim_mat_in
      do j = i, dim_mat_in
        if (j .eq. i ) then
          diagonalA(i) = mat_in(i,j)
        end if
      end do
    end do

    matA = mat_in

!   get initial vector
!   search for column with maximum value and use that column for initial vector

    do i = 1, eigen_in
      idxMaxVal = minloc(diagonalA, mask=mask)
      matV(idxMaxVal, i) = 1.0d0
      mask(idxMaxVal) = .false.
    end do 

    n_grow = eigen_in

!   start loop

    outer: do it = 1, maxiter
      if (verbose .ge. 2) then
        print *
        print *, '----------------------------------------------------------------'
        print *, 'Iteration', it
        print *, '----------------------------------------------------------------'
      end if

      eigenvecs = 0.0d0
      eigenvals = 0.0d0

!     get projection matrix P = (W)^T * V;  W = A * V

      call dgemm('n', 'n', ndimA, n_grow, ndimA, 1.0d0, matA, ndimA, matV, ndimA, 0.0d0, matW, ndimA)
      call dgemm('t', 'n', ndimV, n_grow, ndimA, 1.0d0, matW, ndimA, matV, ndimA, 0.0d0, matP, ndimV)

!
!     diagonalize and obtain eigenvalues and eigenvectors 
!
      eigenvecs = matP
      call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, lw, -1, info)
      lwork = int(lw(1))
      allocate(work(lwork))
      call dsyev('V', 'u', n_grow, eigenvecs, ndimV, eigenvals, work, lwork, info)
      deallocate(work)
      call checkInfo(info, 'Diagonalization of projection Space P')

      if (verbose .ge. 2) then
        print *
        print *, 'Eigenvalues:'
        call printVector(eigenvals, eigen_in)
        print*
      end if

!     get ritz vector 

      call dgemm('n', 'n', ndimA, n_grow, ndimV, 1.0d0, matV, ndimA, eigenvecs, ndimV, 0.0d0, ritzVector, ndimA)


!     get residual

      call dgemm('n', 'n', ndimA, eigen_in, ndimV, 1.0d0, matW, ndimA, eigenvecs, ndimV, 0.0d0, residual, ndimA)
      do i = 1, eigen_in
        call daxpy(ndimA, -eigenvals(i), ritzVector(:,i), 1, residual(:,i), 1)
      end do


!     compute norm and check convergency

      do i = 1, eigen_in
        if (verbose .ge. 2) then
          write(*,'(A, I2, A)', advance='no') 'residual norm of eigenvector', i, ':' 
          write(*,'(E13.3)') dnrm2(ndimA, residual(:,i), 1)
        end if
        if(dnrm2(ndimA, residual(:,i), 1) .le. threshold_residual) then
          if (verbose .ge. 2) then
            print *, 'Converged for:', i
          end if
          converged(i) = .true.
        end if
      end do

      if (all(converged)) then
        if (verbose .ge. 1) then
          print *, 'Converged for all sought eigenpairs after', it, 'iterations.'
        end if


!       copy eigenpairs in output

        do i = 1, eigen_in 
          return_eigenvals(i) = eigenvals(i)
          return_eigenvecs(:,i) = ritzVector(:,i)
        end do
       
        exit outer
      end if
      

!     check if matrix is not full

      matrix_not_full = .true.
      if (n_grow .gt. ndimV - eigen_in) then
        matrix_not_full = .false.
        if (verbose .gt. 2) then
          print *, 'Subspace V is full'
        end if
      end if


!     check if reduced space would exceeds size of input Matrix
      size_not_exceeded = .true.
      if (n_grow + eigen_in .gt. ndimA) then
        size_not_exceeded = .false.  
      end if


!     precondition y_i = r_i / (D - lambda_i)

      if (matrix_not_full .and. size_not_exceeded) then
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
        

!       Gram Schmidt orthogonalization

!       loop implementation

        if (GS_in_loop) then
          do i = 1, eigen_in
            temp = 0.0d0
            do j = 1, n_grow
              temp(:,i) = temp(:,i) + dot_product(residual(:,i), matV(:,j)) / dot_product(matV(:,j), matV(:,j)) * matV(:,j)
            end do
            residual(:,i) = residual(:,i) - temp(:,i)  
            residual(:,i) = residual(:,i) / dnrm2(ndimA, residual(:,i), 1)
          end do 

!       check if orthogonal

          not_orthogonal = .false.
          call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS, not_orthogonal, .true.)
        
            
!       matrix implementation

        else
          do j = 1, max_orth

            temp_mat        = zero
            temp_mat_prime  = zero

!           matrix product V^T * y  
            call dgemm('t', 'n',  n_grow, eigen_in ,ndimA, 1.0d0, matV, ndimA, residual, ndimA, 0.0d0, temp_mat, ndimV)

!           matrix product V * (V^T * y)
            call dgemm('n', 'n', ndimA, eigen_in, n_grow, 1.0d0, matV, ndimA, temp_mat, ndimV, 0.0d0, temp_mat_prime, ndimA)

!           y_prime = y - V * (V^T * y) 
            do i = 1, eigen_in
              residual(:, i) = residual(:, i) - temp_mat_prime(:, i)
            end do

!           print a warning if matrix is not orthogonal 
            not_orthogonal = .false.
            call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS, not_orthogonal, .true.)


!
!         orthogonalize residual with itself (QR)
!
            if (eigen_in .gt. 1) then
!
              allocate(tau(eigen_in))
              tau = zero
              call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, lw, -1, info)
              lwork = int(lw(1))
              allocate(work(lwork))
              work = zero
              call dgeqrf(ndimA, eigen_in, residual, ndimA, tau, work, lwork, info)
              call checkInfo(info, 'Orthogonalization of residual step 1')
!
!
              call dorgqr(ndimA, eigen_in, eigen_in, residual, ndimA, tau, work, lwork, info)
              call checkInfo(info, 'Orthogonalization of residual step 2')
              deallocate(work)
              deallocate(tau)


!             check if orthogonal
              
!             print a warning if matrix is not orthogonal 
              not_orthogonal = .false.
              call checkOrth2mat(residual, eigen_in, 'Residual', matV, n_grow, 'Matrix V', thresh_GS, not_orthogonal, .false.)
              call checkOrth1mat(residual, eigen_in, 'Precondition', thresh_GS, not_orthogonal, .false.)
            end if
            
            if (.not. not_orthogonal) then
              if (verbose .gt. 2) then
                print *, 'orthogonal after', j, 'iterations'
              end if
                exit
            end if
              
            if (j .eq. max_orth) then
              print *
              print *, '*** GS orthogonalization did not result in an orthogonormal matrix after ', max_orth, &
                                                                        'orthogonalizations.***'
              print *
              exit outer
            end if

          end do
            

!         orthonormalization
          do i = 1, eigen_in
            residual(:,i) = residual(:,i) /  dnrm2(ndimA, residual(:,i), 1)
          end do

        end if  




!       add orthonormal residual to subspace

        do i = 1, eigen_in
          matV(:, n_grow  + i) = residual(:,i)
        end do


!       handle case if not converged after maxiter steps

        if (it .eq. maxiter) then
          if (verbose .ge. 1) then
            print *
            print *, 'Not converged after', it, 'iterations'
          end if
          return_eigenvecs = zero
          return_eigenvals = zero
        end if


!       restart
        n_grow = n_grow + eigen_in 


!     treat case if matrix full
      else

        matV            = zero
!
        do i = 1, eigen_in
          matV(:, i) = ritzVector(:, i)
        end do
!
!       restart
!
        matW            = zero
        matP            = zero
        eigenvals       = zero
        eigenvecs       = zero
        ritzVector      = zero
        temp_mat        = zero
        temp_mat_prime  = zero
        residual        = zero
        temp            = zero

        n_grow          =  eigen_in
        
      end if
      
      if (it .eq. maxiter .and. verbose .ge. 1) then
        print *, '** Davidson routine not converged after ', it,' iterations. **' 
      end if

    end do outer

    deallocate(matA, mask, diagonalA, matV, matW, matP, eigenvals, eigenvecs, &
               ritzVector, residual, temp_mat, temp_mat_prime, temp, converged)


    if (verbose .ge. 1) then
      print *
      print *, '--------------------- end of Davidson routine ------------------'
      print *
    end if

  end subroutine symmetricDavidson
!
!
!
  subroutine cholesky(m, n, A, verbose)
  intrinsic                 ::  selected_real_kind, sqrt
  integer,  parameter       ::  wp = selected_real_kind(15)
!
!  takes matrix A with m rows and n columns
!
    integer,    intent(in)    :: m, n, verbose
    real(wp),   intent(inout) :: A(:,:)
    real(wp),   allocatable   :: L(:,:)
    integer                   :: i, k, j 
    real(wp)                  :: sum1, sum2
!
    allocate(L(m,n))
    L = 0.d0
! 
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
!
  if (verbose .gt. 1) then
    call printMatrix(L, m, n)
  end if
  A = L
!
  end subroutine cholesky
!
!
!
  subroutine lu(m, n, A, p, verbose)
  intrinsic                 ::  selected_real_kind, sqrt
  integer,  parameter       ::  wp = selected_real_kind(15)
!
!  takes matrix A with m rows and n columns
!
    integer,    intent(in)    :: m, n, verbose
    real(wp),   intent(inout) :: a(:,:)
    integer,    intent(out)   :: p(:)
!
    real(wp),   allocatable   :: l(:,:)
    integer                   :: i, k, j 
    real(wp)                  :: sum1, sum2
!
    allocate(L(m,n))
    p = (/ (i, i=1,3) /)
    L = 0.d0
    print *, 'Challo'
! 
!
  end subroutine lu
!
!
!
  subroutine lapackDiag(mat, eigenvals, ndimMat)
!
! get diagonalized matrix with lapack function
!
    
    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    real(wp),   intent(inout)     :: mat(:,:)
    real(wp),   intent(out)       :: eigenvals(:)
    integer,    intent(in)        :: ndimMat
    real(wp)                      :: lw(1)
    real(wp), allocatable         :: work(:)
    integer                       :: lwork, info

    call dsyev('V', 'u', ndimMat, mat, ndimMat, eigenvals, lw, -1, info)
    lwork = int(lw(1))
    allocate(work(lwork))
    call dsyev('V', 'u', ndimMat, mat, ndimMat, eigenvals, work, lwork, info)
    call checkInfo(info, 'diagonalize whole mat A')
    deallocate(work)

  end subroutine lapackDiag
!
!
!
  subroutine printMatrix(mat, nrows, ncols) 
!   
! print formatted matrix
!
    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    real(wp), intent(in)  :: mat(:,:)
    integer , intent(in)  :: nrows, ncols
    integer :: i,j 

    do i = 1, ncols
      do j = 1, nrows
        write(*,'(F6.3)', advance='no') mat(i,j)
        if (j .lt. nrows) then
          write(*, '(A)', advance='no') ' '
        end if
      end do
      print *
    end do

  end subroutine printMatrix
!
!
! 
  subroutine printVector(vec, lenRow)
!    
! print formatted vector
!
    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    real(wp), intent(in)  :: vec(:)
    integer,  intent(in)  :: lenRow
    integer               :: i

    do i = 1, lenRow
      !print *, vec(i)
      write(*,'(F6.3)') vec(i)
    end do

  end subroutine printVector
!
!
!
  subroutine checkInfo(info, occasion)
!    
! check if info is zero and print error message if not
!
    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    integer               :: info, zero
    character(len=*)     :: occasion

    zero = 0
    if (info .NE. zero) then
      print *
      print *, '--- WARNING ---'
      print *, occasion
      print *, 'Process terminated with info not equal to 0'
      print*
    end if
!
  end subroutine checkInfo
!
!
  subroutine checkOrth2mat(mat1, ncols1, mat1_name, mat2, ncols2, mat2_name, thresh, not_orthogonal, verbose)
!    
! check if vectors of two matrices are orthogonal
!
    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    integer,              intent(in)      ::  ncols1, ncols2
    real(wp),             intent(in)      ::  thresh, mat1(:,:), mat2(:,:)
    character(len=*),    intent(in)       ::  mat1_name, mat2_name
    real(wp)                              ::  dot_prod
    integer                               :: i,j 
    logical,              intent(inout)   :: not_orthogonal
    logical,              intent(in)      :: verbose

    do i = 1, ncols1
      do j = 1, ncols2
        dot_prod = abs(dot_product( mat2(:,j), mat1(:,i) ))
        if ( dot_prod .GT. thresh ) then
          not_orthogonal = .true.
          if (verbose) then
            print *
            print *, '--- WARNING ---'
            print *, 'vector ', i, ' of matrix ', mat1_name, ' is not orthogonal to vector ', j, ' of matrix', mat2_name
            print *, 'Result of dot product:', dot_prod
            print *
          end if 
        end if
      end do
    end do

  end subroutine checkOrth2mat


  subroutine checkOrth1mat(mat1, ncols1, mat1_name, thresh, not_orthogonal, verbose)
! check if the vectors within the matrix are orthogonal

    implicit none
    intrinsic                 ::  selected_real_kind, abs
    integer,  parameter       ::  wp = selected_real_kind(15)
    integer,              intent(in)      ::  ncols1
    real(wp),             intent(in)      ::  thresh, mat1(:,:)
    character(len=*),     intent(in)      ::  mat1_name
    logical,              intent(inout)   ::  not_orthogonal
    logical,              intent(in)      :: verbose
    real(wp)                              ::  dot_prod
    integer :: i,j 

    do i = 1, ncols1
      do j = 1, ncols1
        if (j .NE. i) then
          dot_prod  = abs(dot_product( mat1(:,j), mat1(:,i) ))
          if ( dot_prod .GT. thresh ) then
            not_orthogonal = .true.
            if (verbose) then 
              print *
              print *, '--- WARNING ---'
              print *, 'vector ', i, ' of matrix ', mat1_name, ' is not orthogonal to vector ', j, ' of matrix', mat1_name
              print *, 'Result of dot product: ', dot_prod
              print *
            end if
          end if
        end if
      end do
    end do

  end subroutine checkOrth1mat

end module diagonalization
