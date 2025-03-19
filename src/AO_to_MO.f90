subroutine AO_to_MO(nBas, c, ERI_AO, ERI_MO)
  implicit none

  ! Input/Output variables
  integer, intent(in) :: nBas
  double precision, intent(in) :: c(nBas,nBas)
  double precision, intent(in) :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision, intent(out) :: ERI_MO(nBas,nBas,nBas,nBas)
  
  ! Local variables
  integer :: mu, nu, la, si
  integer :: p, q, r, s
  double precision, allocatable :: tmp1(:,:,:,:), tmp2(:,:,:,:)
  
  ! Hello World
  write(*,*)
  write(*,*) '**************************************'
  write(*,*) '|       AO to MO Transformation      |'
  write(*,*) '**************************************'
  write(*,*)
  
  ! Memory allocation
  allocate(tmp1(nBas, nBas, nBas, nBas))
  allocate(tmp2(nBas, nBas, nBas, nBas))
  
  ! Initialisation
  tmp1 = 0.0d0
  tmp2 = 0.0d0
  
  ! First quarter-transformation
  do p = 1, nBas
    do mu = 1, nBas
      tmp1(p,:,:,:) = tmp1(p,:,:,:) + c(mu,p)*ERI_AO(mu,:,:,:)
    end do
  end do
  
  ! Second quarter-transformation
  do p = 1, nBas
    do q = 1, nBas
      do nu = 1, nBas
        tmp2(p,q,:,:) = tmp2(p,q,:,:) + c(nu,q)*tmp1(p,nu,:,:)
      end do
    end do
  end do
  
  tmp1 = 0.0d0

  ! Third quarter-transformation
  do p = 1, nBas
    do q = 1, nBas
      do r = 1, nBas
        do la = 1, nBas
          tmp1(p,q,r,:) = tmp1(p,q,r,:) + c(la,r)*tmp2(p,q,la,:)
        end do
      end do
    end do
  end do
  
  tmp2 = 0.0d0

  ! Fourth quarter-transformation
  do p = 1, nBas
    do q = 1, nBas
      do r = 1, nBas
        do s = 1, nBas
          do si = 1, nBas
            tmp2(p,q,r,s) = tmp2(p,q,r,s) + c(si,s)*tmp1(p,q,r,si)
          end do
        end do
      end do
    end do
  end do
  
  ERI_MO = tmp2
  
  ! Deallocate temporary arrays
  deallocate(tmp1, tmp2)
  
  write(*,*) 'AO to MO transformation completed successfully!'
  write(*,*)
  
end subroutine AO_to_MO