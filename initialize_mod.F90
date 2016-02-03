!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE: 
!
MODULE Initialize_Mod
!
! !USES:
!
  ! GEOS-Chem modules
  USE CMN_SIZE_MOD
  USE CMN_GCTM_MOD                        ! Physical constants
  USE ERROR_MOD,        ONLY : ERROR_STOP       
  USE GRID_MOD                           
 
  IMPLICIT NONE
  PRIVATE

!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: Setup_Grid
!
! !REMARKS:
!
! !REVISION HISTORY:
!  30 Jul 2015 - K. Yu       - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_resolution_string
!
! !DESCRIPTION: Function Get\_Resolution\_String returns the proper filename 
!  extension for the GEOS-Chem horizontal grid resolution.  This is used to
!  construct the various file names.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Resolution_String() RESULT( resString )
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255) :: resString
! 
! !REVISION HISTORY:
!  10 Feb 2012 - R. Yantosca - Initial version
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocssor switch
!  14 Jan 2014 - R. Yantosca - Now add NESTED_SE option
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GRID4x5 )
    resString = '4x5.nc'
     
#elif defined( GRID2x25 ) 
    resString = '2x25.nc'

#elif defined( GRID1x125 )
    resString = '1x125.nc'

#elif defined( GRID1x1 ) 
    resString = '1x1.nc'

#elif defined( GRID05x0666 )
    resString = '05x0666.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_CH )
    resString = '025x03125.CH.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_EU )
    resString = '025x03125.EU.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_NA )
    resString = '025x03125.NA.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_SE )
    resString = '025x03125.SE.nc'

#elif defined( GRID025x03125 )
    resString = '025x03125.nc'

#endif

  END FUNCTION Get_Resolution_String
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS checks to see if dimensions read 
!  from the netCDF file match the defined GEOS-Chem dimensions.  If not, then 
!  it will stop the GEOS-Chem simulation with an error message.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Setup_Grid( )
!
! !INPUT PARAMETERS:
!
! 
! !REMARKS:
!
! !REVISION HISTORY:
!  30 Jul 2015 - K. Yu       - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER     :: RC

    ! Start
    CALL INIT_GRID(IM=IIPAR, JM=JJPAR, LM=LLPAR, RC=RC)

#if defined(GRID4x5)
    DLON = 5.0d0
    DLAT = 4.0d0
    DLAT(:,1,:) = 2.0d0
    DLAT(:,JM_WORLD,:) = 2.0d0
#elif defined(GRID2x25)
    DLON = 2.5d0
    DLAT = 2.0d0
    DLAT(:,1,:) = 1.0d0
    DLAT(:,JM_WORLD,:) = 1.0d0
#elif defined(GRID025x03125)
    DLON = 0.3125d0
    DLAT = 0.25d0  
    DLAT(:,1,:) = 0.125d0
    DLAT(:,JM_WORLD,:) = 0.125d0
#endif

    CALL COMPUTE_GRID(I1=1, I2=IIPAR, J1=1, J2=JJPAR, JSP=1, JNP=JM_WORLD, &
                      L1=1, L2=LLPAR, DLON=DLON, DLAT=DLAT, I_LO=I_LO, &
                      J_LO=J_LO, RC=RC)
  END SUBROUTINE Setup_Grid
!EOC
END MODULE Initialize_Mod
