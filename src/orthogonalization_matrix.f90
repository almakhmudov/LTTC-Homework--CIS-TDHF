subroutine orthogonalization_matrix(nBas,S,X)
    implicit none

    ! Input variables
    integer,intent(in)            :: nBas
    double precision,intent(in)   :: S(nBas,nBas)

    ! Local variables
    integer                       :: i
    double precision,allocatable  :: e_vectors(:,:)
    double precision,allocatable  :: e_values(:)
    double precision,allocatable  :: temp_matrix(:,:)
    double precision,allocatable  :: Identity(:,:)

    ! Output variables
    double precision,intent(out)   :: X(nBas,nBas)

    ! Allocate the variables
    allocate(e_vectors(nBas,nBas),e_values(nBas),temp_matrix(nBas,nBas),Identity(nBas,nBas))

    ! The Lowdin orthogonalisation procedure
    e_vectors = S
    call diagonalize_matrix(nBas,e_vectors,e_values)

    ! Construct the diagonal matrix
    temp_matrix = 0.0d0
    do i=1,nBas
        temp_matrix(i,i) = 1.0d0/sqrt(e_values(i))
    end do

    ! X = e_vectors * diag_matrix with e_values^(-1/2) * e_vectors^T
    X = matmul(e_vectors,matmul(temp_matrix,transpose(e_vectors)))

    ! For debugging purposes
    Identity = matmul(X,matmul(S,transpose(X)))
    print*,''
    print*,'Identity matrix I:'
    print*,''
    call matout(nBas,nBas,Identity)
    print*,''

    ! Deallocate the variables
    deallocate(e_vectors,e_values,temp_matrix,Identity)

    ! For debugging purposes
    print*,''
    print*,'Orthogonalization matrix X:'
    print*,''
    call matout(nBas,nBas,X)
    print*,''

end subroutine orthogonalization_matrix