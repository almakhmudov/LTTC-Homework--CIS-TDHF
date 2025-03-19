subroutine TDHF(nBas, nO, e, ERI_MO, nTDHFstates, omega)
    implicit none

    ! Input/Output variables
    integer, intent(in) :: nBas, nO, nTDHFstates
    double precision, intent(in) :: e(nBas), ERI_MO(nBas,nBas,nBas,nBas)
    double precision, intent(out) :: omega(nTDHFstates)
    
    ! Local variables
    integer :: i, j, a, b, ia, jb, nV, nStates
    double precision, allocatable :: A_matrix(:,:), B_matrix(:,:), C_matrix(:,:)
    double precision, allocatable :: AminB(:,:), AplusB(:,:), sqrt_AminB_temp(:,:), AminB_temp(:,:)
    
    ! Hello world
    write(*,*)
    write(*,*) '**************************************'
    write(*,*) '|          TDHF calculation          |'
    write(*,*) '**************************************'
    write(*,*)
    
    ! Number of virtual orbitals
    nV = nBas - nO
    
    ! Dimension of the matrices
    nStates = nO * nV
    
    ! Allocate memory
    allocate(A_matrix(nStates, nStates), B_matrix(nStates, nStates), C_matrix(nStates, nStates))
    allocate(AminB(nStates, nStates), AplusB(nStates, nStates), sqrt_AminB_temp(nStates, nStates), AminB_temp(nStates, nStates))

    ! Initialize matrices
    A_matrix = 0.0d0
    B_matrix = 0.0d0
    
    ! Construct A and B matrices
    do ia = 1, nStates
        i = (ia - 1) / nV + 1
        a = mod(ia - 1, nV) + nO + 1
        
        do jb = 1, nStates
            j = (jb - 1) / nV + 1
            b = mod(jb - 1, nV) + nO + 1
            
            ! A matrix
            ! Diagonal term
            if (i == j .and. a == b) then
                A_matrix(ia,jb) = e(a) - e(i)
            end if
            A_matrix(ia,jb) = A_matrix(ia,jb) + 2.0d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
            
            ! B matrix
            ! Exchange terms
            B_matrix(ia,jb) = B_matrix(ia,jb) + 2.0d0 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a)
        end do
    end do
    
    ! Compute (A - B) and (A + B)
    AminB = A_matrix - B_matrix
    AplusB = A_matrix + B_matrix
    AminB_temp = AminB

    ! Compute (A - B)^(-1/2) using matrix diagonalisation
    call diagonalize_matrix(nStates, AminB_temp, omega)
    do i = 1, nStates
        sqrt_AminB_temp(i,i) = sqrt(omega(i))  ! Square root of eigenvalue
    end do
    AminB = matmul(AminB_temp, matmul(sqrt_AminB_temp, transpose(AminB_temp)))  ! (A-B)^(-1/2)

    ! Compute C matrix
    C_matrix = matmul(AminB, matmul(AplusB, AminB))

    ! Diagonalize C to get excitation energies
    ! Reuse omega array
    omega = 0.0d0
    call diagonalize_matrix(nStates, C_matrix, omega)
    omega = sqrt(omega)  ! Excitation energies
    
    ! Print excitation energies
    write(*,*) '----------------------------------'
    write(*,*) '|    Excitation energies (eV):   |'
    write(*,*) '----------------------------------'
    do i = 1, nTDHFstates
        write(*,'(1X,A1,1X,I5,3X,A1,1X,F20.3,1X,A1)') '|', i, '|', omega(i) * 27.2114d0, '|'
    end do
    write(*,*) '----------------------------------'
    
    ! Deallocate memory
    deallocate(A_matrix, B_matrix, C_matrix, AminB, AplusB, sqrt_AminB_temp, AminB_temp)

end subroutine TDHF