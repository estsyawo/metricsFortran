! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! December 24, 2018
! Conjugate Gradient algorithm for system of linear equations AX = b
! Reference: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

! subroutine gauss_seidel()
subroutine gauss_seidel(A,b,X,n)
! solve the system AX=b
    implicit none
    integer :: n
    real, dimension(n,n), intent(in) :: A
    real, dimension(n,1), intent(in) :: b
    real, dimension(n,1), intent(out) :: X
    real :: sig, xi, dev, mxdev, tol
    integer :: i, j, mxcnt, cnt
    
! initialise elements of X and parameters
    X(1:n,1) = 1e-4
    dev = 0.0
    mxdev = 0.0
    tol = 1e-6
    mxcnt = 1000
    cnt = 0
    
    do !begin while loop
    cnt = cnt + 1
      do i=1,n
      sig = 0.0
      if(abs(A(i,i))<=tol) then
      print*, i, "diagonal element of A too small, floating point error likely"
      exit
      end if
        do j=1,n
          if(j/=i) then
          sig = sig + A(i,j)*X(j,1)
          end if
        end do ! end j-loop
        xi = X(i,1) ! hold X_i
        X(i,1) = 1/A(i,i)*(b(i,1)-sig)
        dev = abs(X(i,1)-xi)
        if(dev>mxdev) then
        mxdev = dev
        end if
      end do ! end i-loop
    if(mxdev<=tol) exit 
    if(cnt>=mxcnt) then
      print*,"Maximum number of iterations reached."
      exit
    end if
    end do
    print*, "Number of iterations = ", cnt, ", deviation=",dev
end subroutine gauss_seidel


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
    call gauss_seidel(A,b,X,n)
    
    print*, "The solution is X="
    print "(1f5.3)", X(:,1)
    
end program main
