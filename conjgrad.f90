! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 8, 2019
! Conjugate Gradient algorithm for system of linear equations AX = b
! Reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method


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
end subroutine conjgrad


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
    call conjgrad(A,b,X,n)
    
    print*, "The solution is X="
    print "(1f5.3)", X(:,1)
    
end program main
