!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: decay_mod
!
!\\
!\\
! !INTERFACE: 
!
      MODULE CHEM_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: do_emiss 
      PUBLIC  :: do_decay 
!
! !PRIVATE MEMBER FUNCTIONS:
! 
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
! !DEFINED PARAMETERS:
!
      ! To convert kg to atoms
      REAL*8, PARAMETER :: XNUMOL_Rn = ( 6.0225d23 / 222.0d-3 )    
      REAL*8, PARAMETER :: XNUMOL_Pb = ( 6.0225d23 / 210.0d-3 )    
      REAL*8, PARAMETER :: XNUMOL_Be = ( 6.0225d23 /   7.0d-3 )

      CONTAINS
!EOC
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE do_emiss ( DT, STT, PS )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE ADVECTION_MOD,      ONLY : GET_AIR_MASS
!
! !INPUT PARAMETERS: 
!
      REAL*8,         INTENT(IN)    :: DT     
      REAL*8,         INTENT(IN)    :: PS(:,:) 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: STT(:,:,:)   
!
! !OUTPUT PARAMETERS:
! 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: I,          J
      REAL*8        :: A_CM2,      Rn_emiss
      REAL*8        :: AD

!$OMP PARALLEL DO           &
!$OMP DEFAULT( SHARED )     &
!$OMP PRIVATE( I, J, A_CM2, Rn_emiss, AD   )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [cm2]
         A_CM2      = GET_AREA_CM2( I, J, 1 )

         ! Compute mass of air in grid box
         AD = GET_AIR_MASS(I,J,1,PS(I,J))

         ! Baseline 222Rn emissions over land [kg]
         ! Rn_LAND [atoms] = [1 atom 222Rn/cm2/s] * [s] * [cm2]
         Rn_emiss    = 1d0 * DT * A_CM2

        ! Convert STT to [atoms Rn] first 
        ! [atoms Rn] = [atoms Rn/atoms air] * [kg air] * [mole air/kg
        ! air] * [atoms air/mole air]
        STT(I,J,1) = STT(I,J,1) * AD * 6.0225d23 / 0.02897

        ! Add to STT array
         STT(I,J,1)     = STT(I,J,1) + Rn_emiss    ! For TURBDAY

        ! Convert STT back to [v/v]
        STT(I,J,1) = STT(I,J,1) * 0.02897 / (AD * 6.0225d23)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

     print*, 'end chem mod', SUM(STT)
      END SUBROUTINE do_emiss
!EOC
!BOC 
      SUBROUTINE do_decay( DT, STT )
!
! !USES:
!
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS:
!
      REAL*8,    INTENT(IN)    :: DT
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,     INTENT(INOUT) :: STT(:,:,:)
!
! !OUTPUT PARAMETERS:
!
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, J, L, N
      REAL*8            :: Tracer_lost(IIPAR, JJPAR, LLPAR)
      REAL*8            :: EXP_Rn
!
! START OF CODE
!
      ! Fraction of Rn left after decay
      EXP_Rn = EXP( -DT * 2.097d-6)

!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L) 
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        Tracer_lost(I,J,L) = STT(I,J,L) * (1d0 - EXP_Rn)
        STT(I,J,L) = STT(I,J,L) - Tracer_lost(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE do_decay
!EOC
      END MODULE CHEM_MOD

