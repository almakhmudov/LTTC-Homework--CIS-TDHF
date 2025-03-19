program EST

  ! Electronic Structure Theory program
    implicit none
    include 'parameters.h'
  
    integer                       :: nAt,nBas,nEl,nO,nV
    double precision              :: ENuc,EHF
  
    double precision,allocatable  :: ZNuc(:)
    double precision,allocatable  :: rAt(:,:)
  
    integer                       :: nShell
    integer,allocatable           :: TotAngMomShell(:)
    integer,allocatable           :: KShell(:)
    double precision,allocatable  :: CenterShell(:,:)
    double precision,allocatable  :: DShell(:,:)
    double precision,allocatable  :: ExpShell(:,:)
  
    double precision,allocatable  :: S(:,:)
    double precision,allocatable  :: T(:,:)
    double precision,allocatable  :: V(:,:)
    double precision,allocatable  :: Hc(:,:)
    double precision,allocatable  :: X(:,:)
    double precision,allocatable  :: ERI(:,:,:,:)
    double precision,allocatable  :: ERI_MO(:,:,:,:)
    double precision,allocatable  :: e(:)
    double precision,allocatable  :: c(:,:)
    double precision,allocatable  :: omega(:)
  
    double precision              :: start_HF,end_HF,t_HF
    double precision              :: start_AOtoMO,end_AOtoMO,t_AOtoMO
    double precision              :: start_CIS,end_CIS,t_CIS
    double precision              :: start_TDHF,end_TDHF,t_TDHF
    integer                       :: nCISstates,nTDHFstates
  
  ! Hello World
  
    write(*,*)
    write(*,*) '***************************'
    write(*,*) '* TCCM winter school 2025 *'
    write(*,*) '***************************'
    write(*,*)
  
  !------------------------------------------------------------------------
  ! Read input information
  !------------------------------------------------------------------------
  
    call read_molecule(nAt,nEl,nO)
    allocate(ZNuc(nAt),rAt(nAt,3))
  
  ! Read geometry
  
    call read_geometry(nAt,ZNuc,rAt,ENuc)
  
    allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
             DShell(maxShell,maxK),ExpShell(maxShell,maxK))
  
  !------------------------------------------------------------------------
  ! Read basis set information
  !------------------------------------------------------------------------
  
    call read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)
  
  !------------------------------------------------------------------------
  ! Read one- and two-electron integrals
  !------------------------------------------------------------------------
  
    allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), &
             ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas))
  
    call read_integrals(nBas,S,T,V,Hc,ERI)  
  
  !------------------------------------------------------------------------
  ! Orthogonalization X = S^(-1/2)
  !------------------------------------------------------------------------
  
    call orthogonalization_matrix(nBas,S,X)
  
  !------------------------------------------------------------------------
  ! Compute restricted HF energy
  !------------------------------------------------------------------------
  
    call cpu_time(start_HF)
    call RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)
    call cpu_time(end_HF)
  
    t_HF = end_HF - start_HF
  
  !------------------------------------------------------------------------
  ! AO to MO transformation
  !------------------------------------------------------------------------
  
    allocate(ERI_MO(nBas,nBas,nBas,nBas))

    call cpu_time(start_AOtoMO)
    call AO_to_MO(nBas,c,ERI,ERI_MO)
    call cpu_time(end_AOtoMO)

    t_AOtoMO = end_AOtoMO - start_AOtoMO
  
  !------------------------------------------------------------------------
  ! CIS Calculation
  !------------------------------------------------------------------------
  
    ! Determine number of CIS states to calculate
    nCISstates = nO * nV
    allocate(omega(nCISstates))
    
    call cpu_time(start_CIS)
    call CIS(nBas, nO, e, ERI_MO, nCISstates, omega)
    call cpu_time(end_CIS)
    
    t_CIS = end_CIS - start_CIS
    deallocate(omega)
  !------------------------------------------------------------------------
  ! TDHF Calculation
  !------------------------------------------------------------------------
    
    ! Determine number of TDHF states to calculate
    nTDHFstates = nO * nV
    ! Reusing omega array
    allocate(omega(nTDHFstates))

    call cpu_time(start_TDHF)
    call TDHF(nBas, nO, e, ERI_MO, nTDHFstates, omega)
    call cpu_time(end_TDHF)
    
    t_TDHF = end_TDHF - start_TDHF

  !------------------------------------------------------------------------
  ! Timing information
  !------------------------------------------------------------------------

    write(*,*)
    write(*,*) '**************************************'
    write(*,*) '|         Timing Information         |'
    write(*,*) '**************************************'
    write(*,*)
    write(*,'(A45,1X,F9.5,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,'(A45,1X,F9.5,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,'(A45,1X,F9.5,A8)') 'Total CPU time for CIS = ',t_CIS,' seconds'
    write(*,'(A45,1X,F9.5,A8)') 'Total CPU time for TDHF = ',t_TDHF,' seconds'
    write(*,*) '--------------------------------------------------------------------------------'
    write(*, '(A45,1X,F9.5,A8)') 'Total CPU time for the program',t_HF + t_AOtoMO + t_CIS + t_TDHF,' seconds'
    write(*,*)
    write(*,*) 'Your calculation is done. Have a great day!'

  !------------------------------------------------------------------------
  ! End of EST
  !------------------------------------------------------------------------
  
  ! Clean up memory
  deallocate(ZNuc, rAt)
  deallocate(CenterShell, TotAngMomShell, KShell, DShell, ExpShell)
  deallocate(S, T, V, Hc, X, ERI, ERI_MO, e, c, omega)

end program EST