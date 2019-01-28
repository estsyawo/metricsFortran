! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 28, 2019
! Linear regression based on the SOR algorithm
! Reference: https://en.wikipedia.org/wiki/Successive_over-relaxation

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
    call linreg_sor(X,Y,beta,n,k)
    
    print*, "The solution is beta ="
    print "(1f5.3)", beta(:,1)

contains

! system of equations based on SOR
! A - nxn matrix of coefficients
! b - nx1 matrix
! X - nx1 matrix of unknowns to be solved for
! n - the number of unknowns in X

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

! linear regression based on SOR
! X - nxk matrix of covariates
! Y - nx1 matrix of outcome 
! beta - kx1 matrix of parameters to be estimated
! n - number of observations
! k - number of parameters (including the intercept)

subroutine linreg_sor(X,Y,beta,n,k)
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
    
    call SOR(A,b,beta,k)
     
end subroutine linreg_sor
    
end program main
