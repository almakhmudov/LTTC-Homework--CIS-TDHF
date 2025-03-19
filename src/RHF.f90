subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)
    ! Perform a restricted Hartree-Fock calculation
    implicit none
    include 'parameters.h'
    
    ! Input variables
    integer,intent(in) :: nBas
    integer,intent(in) :: nO
    double precision,intent(in) :: S(nBas,nBas)
    double precision,intent(in) :: T(nBas,nBas)
    double precision,intent(in) :: V(nBas,nBas)
    double precision,intent(in) :: Hc(nBas,nBas)
    double precision,intent(in) :: X(nBas,nBas)
    double precision,intent(in) :: ERI(nBas,nBas,nBas,nBas)
    double precision,intent(in) :: ENuc
    
    ! Local variables
    integer :: i,j,k,l, mu,nu,si,la
    integer,parameter :: maxSCF = 64
    double precision,parameter :: thresh = 1d-5
    integer :: nSCF
    double precision :: Conv
    double precision :: Gap
    double precision :: ET,EV,EJ,EK
    double precision,allocatable :: cp(:,:)
    double precision,allocatable :: P(:,:)
    double precision,allocatable :: J_en(:,:)
    double precision,allocatable :: K_en(:,:)
    double precision,allocatable :: F(:,:),Fp(:,:)
    double precision,allocatable :: error(:,:)
    double precision,allocatable :: tmp(:,:)
    double precision,external :: trace_matrix
    
    ! Output variables
    double precision,intent(out) :: EHF
    double precision,intent(out) :: e(nBas)
    double precision,intent(out) :: c(nBas,nBas)
    
    ! Hello world
    write(*,*)
    write(*,*)'************************************************'
    write(*,*)'|      Restricted Hartree-Fock calculation     |'
    write(*,*)'************************************************'
    write(*,*)
    
    ! Memory allocation
    allocate(cp(nBas,nBas),P(nBas,nBas), &
             J_en(nBas,nBas),K_en(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
             error(nBas,nBas),tmp(nBas,nBas))
    
    ! Initialize variables
    F(:,:) = Hc(:,:)
    nSCF = 0
    Conv = 1.0d0
    
    ! Main SCF loop header
    write(*,*)
    write(*,*)'----------------------------------------------------'
    write(*,*)'|                 RHF calculation                  |'
    write(*,*)'----------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
    '|','#','|','HF energy','|','Conv','|','HL Gap','|'
    write(*,*)'----------------------------------------------------'
    
    do while(Conv > thresh .and. nSCF < maxSCF)
  
        ! Increment iteration counter
        nSCF = nSCF + 1
  
        ! Transform Fock matrix
        Fp = matmul(transpose(X), matmul(F, X))
        
        ! Diagonalize F' to get C' and energies e
        cp = Fp
        call diagonalize_matrix(nBas, cp, e)
        
        ! Transform C' back to original basis
        c = matmul(X, cp)
        
        ! Form new density matrix
        P = 2.0d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
        
        ! Build J and K matrices
        J_en = 0d0
        K_en = 0d0
        do mu = 1, nBas
          do nu = 1, nBas
              do la = 1, nBas
                  do si = 1, nBas
                      J_en(mu,nu) = J_en(mu,nu) + P(la,si) * ERI(mu,la,nu,si)
                      K_en(mu,nu) = K_en(mu,nu) + P(la,si) * ERI(mu,la,si,nu)
                  end do
              end do
          end do
        end do
      
        ! Form new Fock matrix F = Hc + J - 1/2·K
        F = Hc + J_en - 0.5d0 * K_en
  
        ! Calculate energy components
        ET = trace_matrix(nBas, matmul(P, T))
        EV = trace_matrix(nBas, matmul(P, V))
        EJ = 0.5d0 * trace_matrix(nBas, matmul(P, J_en))
        EK = 0.5d0 * trace_matrix(nBas, matmul(P, K_en))
        
        ! Calculate current HF energy
        EHF = 0.5 * trace_matrix(nBas, matmul(P, F + Hc))
        
        ! Convergence criteria
        error = matmul(matmul(F, P), S) - matmul(matmul(S, P), F)
        Conv = maxval(abs(error))  ! Added abs() for proper convergence check
        
        ! Calculate HOMO-LUMO gap
        Gap = e(nO+1) - e(nO)
        
        ! Print iteration results
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
        '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
        
    end do
    
    write(*,*)'----------------------------------------------------'
    
    ! Check convergence
    if(nSCF == maxSCF) then
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'             Convergence failed                     '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        stop
    endif
    
    ! Compute final HF energy
    call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)
    
    ! Deallocate arrays
    deallocate(cp,P,J_en,K_en,F,Fp,error,tmp)
    
  end subroutine RHF