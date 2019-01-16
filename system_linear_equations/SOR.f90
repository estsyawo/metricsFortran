! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 15, 2019
! Successive over-relaxation (SOR) algorithm for system of linear equations AX = b
! Reference: https://en.wikipedia.org/wiki/Successive_over-relaxation


subroutine SOR(A,b,X,n)
! solve the system AX=b
    implicit none
    integer :: n
    real, dimension(n,n), intent(in) :: A
    real, dimension(n,1), intent(in) :: b
    real, dimension(n,1), intent(out) :: X
    real :: sig, tol, w, xi, dev, mxdev
    integer :: i, j, k, maxiter

! initialise parameters
    w=1.2
    maxiter = 1000
    k = 0
    X(1:n,1) = 0.0
    dev = 0.0
    mxdev = 0.0
    tol = 1e-07
    
! main while loop
    do 
    k = k+1
      do i=1,n
      sig = 0.0
        do j=1,n
          if(i.ne.j) then
            sig = sig + A(i,j)*X(j,1)
          end if
        end do ! end do j
        xi = X(i,1)
        X(i,1) = X(i,1) + w*((b(i,1)-sig)/A(i,i) - X(i,1))
        dev = abs(X(i,1)-xi)
        if(dev>mxdev) then
        mxdev = dev
        end if
      end do ! end do i
    
    !check convergence
    if(mxdev<=tol) then
    exit
    end if
    if(k>=maxiter) then
    print*, "Maximum number of iterations reached."
    exit
    end if
    end do !end main while loop



end subroutine SOR

! main programme to execute subroutine
program main
    ! main program for execution
    implicit none
    integer, parameter :: n = 10
    real, dimension(n,n):: A
    real, dimension(100,n) :: M
    real, dimension(n,1) :: b
    real, dimension(n,1) :: X
    real, dimension(n,1) :: Z
    integer :: i
   
    
    call random_number(M) ! randomly fill matrix A    
    A = matmul(transpose(M),M)
    
    call random_number(Z)
    b = matmul(A,Z) !compute b
    
    print *, "Matrix A is "
    print "(10f5.1)", (A(i,:),i=1,n)
    print*
    print*, "The vector b is"
    print "(1f5.1)", b(:,1)
    print*
    print*, "The true unknown vector to be solved for is Z = "
    print "(1f5.3)", Z(:,1)
    print*
    call SOR(A,b,X,n)
    
    print*, "The solution is X="
    print "(1f5.3)", X(:,1)
    
end program main
