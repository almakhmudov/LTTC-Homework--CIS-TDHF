subroutine CIS(nBas, nO, e, ERI_MO, nCISstates, omega)
    implicit none

    ! Input/Output variables
    integer, intent(in) :: nBas              ! Number of basis functions
    integer, intent(in) :: nO                ! Number of occupied orbitals
    double precision, intent(in) :: e(nBas)  ! Orbital energies
    double precision, intent(in) :: ERI_MO(nBas,nBas,nBas,nBas) ! Two-electron integrals in MO basis
    integer, intent(in) :: nCISstates        ! Number of CIS states to compute
    double precision, intent(out) :: omega(nCISstates) ! Excitation energies
    
    ! Local variables
    integer :: i, j, a, b, ia, jb
    integer :: nV, nStates
    double precision, allocatable :: A_matrix(:,:)

    ! Hello World
    write(*,*)
    write(*,*) '**************************************'
    write(*,*) '|          CIS calculation           |'
    write(*,*) '**************************************'
    write(*,*)
    
    ! Number of virtual orbitals
    nV = nBas - nO
    
    ! Dimension of the CIS matrix
    nStates = nO * nV
    
    ! Allocate memory for the CIS matrix
    allocate(A_matrix(nStates, nStates))
    
    ! Construct the CIS matrix
    A_matrix = 0.0d0
    do ia = 1, nStates
        i = (ia - 1) / nV + 1
        a = mod(ia - 1, nV) + nO + 1
        
        do jb = 1, nStates
            j = (jb - 1) / nV + 1
            b = mod(jb - 1, nV) + nO + 1
            
            ! Diagonal term
            if (i == j .and. a == b) then
                A_matrix(ia,jb) = e(a) - e(i)
            end if
            
            ! Exchange terms
            A_matrix(ia,jb) = A_matrix(ia,jb) + 2.0d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a)
        end do
    end do
    
    ! Call the diagonalisation subroutine
    call diagonalize_matrix(nStates, A_matrix, omega)
    
    ! Print the excitation energies in eV
    write(*,*) '----------------------------------'
    write(*,*) '|    Excitation energies (eV):   |'
    write(*,*) '----------------------------------'
    
    do i = 1, nCISstates
        write(*,'(1X,A1,1X,I5,3X,A1,1X,F20.3,1X,A1)') &
              '|', i, '|', omega(i)*27.2114d0, '|'
    end do
    write(*,*) '----------------------------------'
    
    ! Deallocate memory
    deallocate(A_matrix)
    
end subroutine CIS