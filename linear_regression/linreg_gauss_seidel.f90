! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 28, 2019
! Linear regression based on the Guass-Seidel algorithm
! Reference: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

! main programme to execute subroutine
program main
    ! main program for execution
    implicit none
    integer, parameter :: k = 5
    integer, parameter :: n = 1000
    real, dimension(n,k) :: X
    real, dimension(n,1) :: Y,e
    real, dimension(k,1) :: Z,beta !true parameter values to be solved for
    integer :: i
   
    
    call random_number(X) ! randomly fill matrix X    
    X(:,1) = 1.0 !1's in first column for intercept
    
    call random_number(e) ! generate error
    e = e/100.0 !scale down error
    call random_number(Z) ! generate true parameters
    Y = matmul(X,Z) + e ! generate outcome variable Y
    
    print*
    print*, "The true unknown vector to be solved for is Z = "
    print "(1f5.3)", Z(:,1)
    print*
    call linreg_gauss_seidel(X,Y,beta,n,k)
    
    print*, "The solution is beta ="
    print "(1f5.3)", beta(:,1)

    
contains

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

! linear regression based on Gauss-Seidel algorithm
! X - nxk matrix of covariates
! Y - nx1 matrix of outcome 
! beta - kx1 matrix of parameters to be estimated
! n - number of observations
! k - number of parameters (including the intercept)

subroutine linreg_gauss_seidel(X,Y,beta,n,k)
    implicit none
    integer :: n,k
    real, dimension(n,k), intent(in) :: X
    real, dimension(n,1), intent(in) :: Y
    real, dimension(k,1), intent(out) :: beta
    real, dimension(k,k) :: A
    real, dimension(k,1) :: b
    
    ! construct kxk Gram matrix
    A = matmul(transpose(X),X)
    b = matmul(transpose(X),Y)
    
    call gauss_seidel(A,b,beta,k)
     
end subroutine linreg_gauss_seidel
    
end program main
