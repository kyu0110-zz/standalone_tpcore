!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cleanup
!
! !DESCRIPTION: Subroutine CLEANUP deallocates the memory assigned to 
!  dynamically allocatable arrays just before exiting a GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP( am_I_Root, RC )
!
! !USES:
!
      USE ERROR_MOD,               ONLY : DEBUG_MSG
      USE TPCORE_FVDAS_MOD,        ONLY : EXIT_TPCORE
      USE CMN_SIZE_Mod,            ONLY : Cleanup_CMN_SIZE

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC           ! Success or failure
!
! !REVISION HISTORY: 
!  29 Nov 1999 - R. Yantosca - Initial version
!  (1 ) CLEANUP is written in Fixed-Format F90.
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) 
 100  FORMAT( '     - CLEANUP: deallocating arrays now...' )

      ! Call cleanup routines from individual F90 modules
      CALL EXIT_TPCORE

      ! Call extra cleanup routines, mostly from modules in Headers/
      CALL Cleanup_CMN_SIZE   ( am_I_Root, RC )

      END SUBROUTINE CLEANUP
!EOC
