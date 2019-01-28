! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 28, 2019
! Linear regression based on the Conjugate gradient algorithm
! Reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method

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
    call linreg_conjgrad(X,Y,beta,n,k)
    
    print*, "The solution is beta ="
    print "(1f5.3)", beta(:,1)

    
contains

! subroutine conjgrad()
subroutine conjgrad(A,b,X,n)
! solve the system AX=b
    implicit none
    integer :: n
    real, dimension(n,n), intent(in) :: A
    real, dimension(n,1), intent(in) :: b
    real, dimension(n,1), intent(out) :: X
    real, dimension(n,1) :: r, p, vec
    real :: dev, alf, beta, tol, dpr1, dpr0
    integer :: i, k, maxiter

! initialise elements of X and parameters
    X(1:n,1) = 0.0
    dev = 0.0
    tol = 1e-7
    maxiter = 1000
    k = 0
    ! fill in initial values
    r = b ! initialise with X=0 vector
    p = r ! initialise p with r
    dpr0 = dot_product(r(:,1),r(:,1)) ! take dot product
    
    do ! begin (while) do loop
    vec = matmul(A,p); ! assign vec = A*p
    dpr1 = dot_product(p(:,1),vec(:,1)) ! store dot product temporarily in dpr1
    alf = dpr0/dpr1 ! compute alpha
    X = X + alf*p ! update X
    r = r - alf*vec ! update r
    vec = abs(r) ! pass absolute values of r to vec
    !dev = vec(maxloc(vec(:,1)),1) !extract infinity norm of r 
    dev = maxval(vec)
    
    ! checking for convergence
    if(dev<tol) then
    exit
    end if 
    
    if(k>=maxiter) then
    print*,"Maximum number of iterations reached. CG algorithm fails to", &
            "converge."
            exit
    end if
    
    dpr1 = dot_product(r(:,1),r(:,1))
    beta = dpr1/dpr0
    p = r + beta*p !update p
    k = k +1
    dpr0 = dpr1 
    
    end do ! end (while) do loop
    print*, "Number of iterations = ",k
end subroutine conjgrad

! linear regression based on conjugate gradient algorithm
! X - nxk matrix of covariates
! Y - nx1 matrix of outcome 
! beta - kx1 matrix of parameters to be estimated
! n - number of observations
! k - number of parameters (including the intercept)

subroutine linreg_conjgrad(X,Y,beta,n,k)
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
    
    call conjgrad(A,b,beta,k)
     
end subroutine linreg_conjgrad

    
end program main
