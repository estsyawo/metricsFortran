! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! May 1, 2019
! Cholesky decomposition of a symmetric positive definite matrix




subroutine choldecom(A,n)
    implicit none
    integer::n,j
    real, dimension(n,n), intent(inout)::A

    do j=1:n
    
    if(j>1) then
    A(j:n,j) = A(j:n,j) - A(j:n,1:(j-1))*transpose(A(j,1:(j-1)))
    
    end if
    A(j:n,j) = A(j:n,j)/sqrt(A(j,j))
    end do


end subroutine
