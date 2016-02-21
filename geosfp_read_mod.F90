!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: read_input_mod
!
! !DESCRIPTION: Module read_input_mod contains subroutines for reading the 
!  GEOS-FP data and restart files from disk (in netCDF format). Modifed from
!  GEOSFP_READ_MOD of the GEOS-Chem code.
!\\
!\\
! !INTERFACE: 
!
MODULE read_input_mod
!
! !USES:
!
  ! NcdfUtil modules for netCDF I/O
  USE m_netcdf_io_open                    ! netCDF open
  USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
  USE m_netcdf_io_read                    ! netCDF data reads
  USE m_netcdf_io_close                   ! netCDF close

  ! GEOS-Chem modules
  USE CMN_SIZE_MOD
  USE CMN_GCTM_MOD                        ! Physical constants
  USE ERROR_MOD,        ONLY : ERROR_STOP, ALLOC_ERR, GEOS_CHEM_STOP       
  USE TIME_MOD                            ! Date and time routines

  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"                    ! Include file for netCDF library
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Check_Dimensions
  PRIVATE :: Get_Resolution_String
  PRIVATE :: LUMP_2
  PRIVATE :: TRANSFER_3D
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: GeosFp_Read_A3dyn
  PUBLIC  :: geosfp_read_omega
  PUBLIC  :: read_wc
  PUBLIC  :: read_Kz
  PUBLIC  :: GeosFp_Read_I3
  PUBLIC  :: restart_read
!
! !REMARKS:
!  Assumes that you have a netCDF library (either v3 or v4) installed on 
!  your system. 
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  03 Feb 2012 - R. Yantosca - Add Geos57_Read_A3 wrapper function
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Add function Get_Resolution_String
!  05 Apr 2012 - R. Yantosca - Convert units for specific humidity properly
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields via Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to geosfp_read_mod.F90
!  14 Jan 2014 - R. Yantosca - Remove "define GEOS572_FILES #ifdef blocks
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
  SUBROUTINE Check_Dimensions( lon, lat, lev, time, time_expected, caller )
!
! !INPUT PARAMETERS:
!
    INTEGER,          OPTIONAL, INTENT(IN)  :: lon            ! Lon dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: lat            ! Lat dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: lev            ! Alt dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: time           ! Time dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: time_expected  ! Expected # of 
                                                              !  time slots
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: caller         ! Name of caller
                                                              !  routine
! 
! !REMARKS:
!  Call this routine with keyword arguments, e.g
!     CALL CHECK_DIMENSION( lon=X,  lat=Y,           lev=Z,         &
!                           time=T, time_expected=8, caller=caller )
!
! !REVISION HISTORY:
!  02 Feb 2012 - R. Yantosca - Initial version
!  03 Feb 2012 - R. Yantosca - Now pass the caller routine name as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Error message string
    CHARACTER(LEN=255) :: errMsg                  

    ! Error check longitude dimension 
    IF ( PRESENT( lon ) ) THEN
       IF ( lon /= IIPAR ) THEN
          WRITE( errMsg, 100 ) lon, IIPAR
 100      FORMAT( 'Longitude dimension (', i5, &
                  ' ) does not match IIPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lat ) ) THEN
       IF ( lat /= JJPAR ) THEN
          WRITE( errMsg, 110 ) lat, JJPAR
 110      FORMAT( 'Latitude dimension (', i5, &
                  ' ) does not match JJPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lev ) ) THEN
       IF ( lev /= LGLOB ) THEN
          WRITE( errMsg, 120 ) lev, LGLOB
 120      FORMAT( 'Levels dimension (', i5, &
                  ' ) does not match LGLOB ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF

    ! Error check time dimension 
    IF ( PRESENT( time ) .and. PRESENT( time_expected ) ) THEN
       IF ( time /= time_expected ) THEN
          WRITE( errMsg, 130 ) time, time_expected
 130      FORMAT( 'Time dimension (', i5, &
                  ' ) does not match expected # of times ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF

  END SUBROUTINE Check_Dimensions
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3dyn
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr time-averaged (A3) data (dynamics fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A3dyn( YYYYMMDD, HHMMSS, wz, U, V )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,         INTENT(INOUT)   :: wz(:,:,:)         ! omega
    REAL*8,         INTENT(INOUT)   :: U(:,:,:)         ! U
    REAL*8,         INTENT(INOUT)   :: V(:,:,:)         ! V
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directories with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3dyn
!  15 Nov 2013 - R. Yantosca - Now convert RH from [1] to [%], in order
!                              to be consistent with GEOS-Chem convention
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                             
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q (IIPAR,JJPAR,LGLOB  )  ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

        print*, 'Reading in met data'
    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A3dyn (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( 'GEOS_FP/YYYY/MM/' )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
#if defined ( GRID025x03125 )
    nc_file = 'GEOS.fp.asm.tavg3_3d_asm_Nv.YYYYMMDD_hhmm.V01.nc4'
#else
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A3dyn.' // TRIM( nc_file )
#endif
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
#if defined ( GRID4x5 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/4x5/' ) // TRIM( dir ) // TRIM( nc_file )
#elif defined ( GRID2x25 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/2x25/' ) // TRIM( dir ) // TRIM( nc_file )
#elif defined ( GRID025x03125 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/Native/' ) // TRIM(dir) // TRIM(nc_file)
#endif    

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
#if defined ( GRID025x03125 )
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=1, caller=caller )
#else
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )
#endif

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
#if defined ( GRID025x03125 ) 
    time_index = 1
#else
    time_index = ( HHMMSS / 030000 ) + 1
#endif
        print*, 'time index', time_index

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, 'GEOS57_READ_A1 (geos57_read_mod.F90)' )
    ENDIF

    !--------------------------------
    ! Read data on level centers
    !--------------------------------

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read OMEGA  from file
    v_name = "OMEGA"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
#if defined ( GRID025x03125 ) 
    call TRANSFER_3D(Q(:,:,LGLOB:1:-1), wz)
#else
    call TRANSFER_3D(Q, wz)
#endif

    ! Read U
    v_name = "U"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
#if defined ( GRID025x03125 )
    call TRANSFER_3D(Q(:,:,LGLOB:1:-1), U)
#else
    call TRANSFER_3D(Q, U)
#endif

    ! Read V
    v_name = "V"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
#if defined ( GRID025x03125 ) 
    call TRANSFER_3D(Q(:,:,LGLOB:1:-1), V)
#else
    call TRANSFER_3D(Q, V)
#endif

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  GEOS-FP A3dyn  met fields for ', a )

    !======================================================================
    ! Unit conversions, diagnostics, cleanup, and quit
    !======================================================================
    
    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_Read_A3dyn
!EOC
!BOP
  SUBROUTINE geosfp_read_omega( YYYYMMDD, HHMMSS, negomega )
!
! !INPUT PARAMTERS:
! 
  INTEGER,      INTENT(IN)      :: YYYYMMDD
  INTEGER,      INTENT(IN)      :: HHMMSS
! 
! Output paramters:
!
  REAL*8,       INTENT(OUT)     :: negomega(:,:,:)
!EOP
!BOC
! 
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q(IIPAR,JJPAR,LGLOB  )    ! 2D temporary data arrray

    !==========================================================================
    ! Opn the netCDF file
    !==========================================================================

    caller = 'GEOSFP_READ_OMEGA'

    ! Replace time & date tokens in the file name
    dir     = TRIM( 'GEOS_FP/YYYY/MM/' )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.OMEGA.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
#if defined( GRID4x5 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/4x5/' ) // TRIM( nc_file )
#elif defined( GRID2x25 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/2x25/' ) // TRIM(nc_file)
#elif defined(  GRID025x03125 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/025/' ) // TRIM(nc_file)
#endif

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st4d   = (/ 1,     1,     1,       time_index /)
    ct4d   = (/ IIPAR, JJPAR, LGLOB,   1          /)

    ! Read OMEGA
    v_name = "negOmega"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d(Q, negomega)

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found negomega GEOS-FP I3     met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE geosfp_read_omega
!EOC
!BOP
  SUBROUTINE read_wc( YYYYMMDD, HHMMSS, wc )
!
! !INPUT PARAMTERS:
! 
  INTEGER,      INTENT(IN)      :: YYYYMMDD
  INTEGER,      INTENT(IN)      :: HHMMSS
! 
! Output paramters:
!
  REAL*8,       INTENT(OUT)     :: wc(:,:,:)
!EOP
!BOC
! 
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q(IIPAR,JJPAR,LLPAR  )    ! 2D temporary data arrray

    !==========================================================================
    ! Opn the netCDF file
    !==========================================================================

    caller = 'READ_WC'

    ! Replace time & date tokens in the file name
    !dir     = TRIM( '' )
    !CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    !nc_file = Get_Resolution_String()
    !nc_file = 'GEOSFP.YYYYMMDD.OMEGA.' // TRIM( nc_file )
    !CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    nc_file = 'wc_YYYYMMDD.nc'
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
#if defined( GRID4x5 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/wc/4x5/' ) // TRIM( nc_file )
#elif defined( GRID2x25 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/2x25/' ) // TRIM('wc.nc')
#elif defined(  GRID025x03125 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/025/' ) // TRIM('wc.nc')
#endif

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = 1

    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st4d   = (/ 1,     1,     1,       time_index /)
    ct4d   = (/ IIPAR, JJPAR, LLPAR,   1          /)

    ! Read OMEGA
    v_name = "wc"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    wc = Q

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found wc met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE read_wc
!EOC
!-----------------------------------------------------------------
!BOP
  SUBROUTINE read_Kz( YYYYMMDD, HHMMSS, Kz )
!
! !INPUT PARAMTERS:
! 
  INTEGER,      INTENT(IN)      :: YYYYMMDD
  INTEGER,      INTENT(IN)      :: HHMMSS
! 
! Output paramters:
!
  REAL*8,       INTENT(OUT)     :: Kz(:,:,:)
!EOP
!BOC
! 
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q(IIPAR,JJPAR,LLPAR  )    ! 2D temporary data arrray

    !==========================================================================
    ! Opn the netCDF file
    !==========================================================================

    caller = 'READ_Kz'

    ! Replace time & date tokens in the file name
    !dir     = TRIM( '' )
    !CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    !nc_file = Get_Resolution_String()
    !nc_file = 'GEOSFP.YYYYMMDD.OMEGA.' // TRIM( nc_file )
    !CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    nc_file = 'YYYYMMDDhhmmss.nc'
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
#if defined( GRID4x5 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/Kz/' ) // TRIM( nc_file )
#elif defined( GRID2x25 )
    nc_file = TRIM( '/n/regal/jacob_lab/kyu/MAX_OMEGA/2x25/' ) // TRIM('wc.nc')
#endif

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = 1

    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st4d   = (/ 1,     1,     1,       time_index /)
    ct4d   = (/ IIPAR, JJPAR, LLPAR,   1          /)

    ! Read OMEGA
    v_name = "Kz"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    Kz = Q

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found wc met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE read_Kz
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_I3_1
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_I3( YYYYMMDD, HHMMSS, PS )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,         INTENT(OUT)   :: PS(:,:)
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Now convert QV1 from [kg/kg] to [g/kg]
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  06 Sep 2013 - R. Yantosca - Bug fix: we need to initialize State_Met%T
!                              with State_Met%TMPU1 to avoid errors.  The
!                              State_Met%T field will be set again in INTERP.
!  29 Oct 2013 - R. Yantosca - Now read T_FULLGRID_1 for offline simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR      )    ! 2D temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

        print*, 'Reading I3 met fields'
    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_I3 (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( 'GEOS_FP/YYYY/MM/' )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
#if defined( GRID025x03125 )
    nc_file = 'GEOS.fp.asm.inst3_3d_asm_Nv.YYYYMMDD_hhmm.V01.nc4'
#else
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.I3.' // TRIM( nc_file )
#endif
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
#if defined( GRID4x5 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/4x5/' ) // TRIM( dir ) // TRIM( nc_file )
#elif defined( GRID2x25 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/2x25/' ) // TRIM(dir) // TRIM(nc_file)
#elif defined(  GRID025x03125 )
    nc_file = TRIM( '/n/home11/kyu/globaltransport/Met/Native/' ) // TRIM(dir) // TRIM(nc_file)
#endif

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
#if defined( GRID025x03125 )
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=1, caller=caller )
#else
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )
#endif

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
#if defined( GRID025x03125 )
    time_index = 1
#else
    time_index = ( HHMMSS / 030000 ) + 1
#endif

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     time_index /)
    ct3d   = (/ IIPAR, JJPAR, 1          /)

    ! Read PS
    v_name = "PS"
    CALL NcRd( Q2, fId, TRIM(v_name), st3d, ct3d )
#if defined( GRID025x03125 )
    PS = Q2 / 100d0
#else
    PS = Q2 
#endif

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found PS  GEOS-FP I3     met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

  END SUBROUTINE GeosFp_Read_I3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_read
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE restart_read( fname, t_inst )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    CHARACTER(LEN=*), INTENT(IN) :: fname  ! name of file to read 
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,         INTENT(INOUT)   :: t_inst(:,:,:)
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR, LLPAR)    ! 2D temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

        print*, 'Reading restart'
    ! Name of this routine (for error printout)
    caller  = "restart_read (input_read_mod.F90)"

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( fname ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=LGLOB,         &
                           time=T, time_expected=1, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     1 /)
    ct3d   = (/ IIPAR, JJPAR, LLPAR /)

    ! Read bottom tracer
    v_name = "tracer"
    CALL NcRd( Q2, fId, TRIM(v_name), st3d, ct3d )
    t_inst(:,:,:) = Q2 

    ! Echo info
 10 FORMAT( '     - Found restart ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE restart_read 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_3d
!
! !DESCRIPTION: Subroutine TRANSFER\_3D transfers 3-dimensional data from a 
!  REAL*4  array to a REAL*8 array.  Vertical layers are collapsed (from LGLOB
!  to LLPAR) if necessary.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0, J0 are local variables. (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  (6 ) Rewritten for clarity (bmy, 2/8/07)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      INTEGER :: L_COPY = 36
      REAL*4  :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_3D begins here!
      !================================================================

      ! Copy the first L_COPY levels
      OUT(:,:,1:L_COPY) = IN( 1:IIPAR, 1:JJPAR, 1:L_COPY )

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Collapse levels in the stratosphere
      !================================================================

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Copy a vertical column into INCOL
         INCOL = IN( I, J, 1:LGLOB )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(I,J,37) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(I,J,38) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(I,J,39) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(I,J,40) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(I,J,41) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(I,J,42) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(I,J,43) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(I,J,44) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(I,J,45) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(I,J,46) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(I,J,47) = LUMP_4( INCOL, LGLOB, 69 ) 

      ENDDO
      ENDDO
      END SUBROUTINE TRANSFER_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_2
!
! !DESCRIPTION: Function LUMP\_2\_R8 lumps 2 sigma levels into one thick 
!  level.  Input arguments must be REAL*8. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_2( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*8              :: OUT        ! Data on output grid: 2 lumped levels
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      REAL*8            :: EDGE_IN(LGLOB+1)

      !=================================================================
      ! LUMP_2_R8 begins here!
      !=================================================================      

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+2 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR: L < 1 or L > L_IN or L+2 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_2 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      EDGE_IN = (/ & 
              0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
              1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
              4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
              7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
              1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
              1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
              2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
              2.243630d+02, 2.168650d+02, 2.011920d+02, & 
!------- EDGES OF GEOS-5 FIXED PRESSURE LEVELS OCCUR BELOW THIS LINE ------
                                                        1.769300d+02, &
              1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
              7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01, &
              4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01, &
              1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01, &
              9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00, &
              4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00, &
              1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01, &
              6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01, &
              2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02, &
              6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02, &
              1.000000d-02 /)

      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) +
      !               ( Q(L+2) * ( PEDGE(L+2) - PEDGE(L+3) ) ) +
      !               ( Q(L+3) * ( PEDGE(L+3) - PEDGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+4)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by EQUATION 2:
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) ]
      !  Q(L') = ----------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+2)
      !=================================================================     

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Othewise,   EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) + &
             ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) )

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+2) )
       
      END FUNCTION LUMP_2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_4_r4
!
! !DESCRIPTION: Function LUMP\_4\_R4 lumps 4 sigma levels into one thick 
!  level.  Input arguments must be REAL*4.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_4( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*4              :: OUT        ! Data on output grid: 4 lumped levels
!
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      REAL*8            :: EDGE_IN(LGLOB+1)
      !=================================================================
      ! LUMP_4_R4 begins here!
      !=================================================================      

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+4 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         WRITE( 6, '(a)' ) 'ERROR: L < 1 or L > L_IN or L+4 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_4 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         CALL GEOS_CHEM_STOP
      ENDIF

      EDGE_IN = (/ &
              0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
              1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
              4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
              7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
              1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
              1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
              2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
              2.243630d+02, 2.168650d+02, 2.011920d+02, &
!------- EDGES OF GEOS-5 FIXED PRESSURE LEVELS OCCUR BELOW THIS LINE ------
                                                        1.769300d+02, &
              1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
              7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01, &
              4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01, &
              1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01, &
              9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00, &
              4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00, &
              1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01, &
              6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01, &
              2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02, &
              6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02, &
              1.000000d-02 /)

      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) +
      !               ( Q(L+2) * ( PEDGE(L+2) - PEDGE(L+3) ) ) +
      !               ( Q(L+3) * ( PEDGE(L+3) - PEDGE(L+4) ) ) ]
      !  Q(L') = --------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+4)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by EQUATION 2:
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) +
      !             ( Q(L+2) * ( SIG_EDGE(L+2) - SIG_EDGE(L+3) ) ) +
      !             ( Q(L+3) * ( SIG_EDGE(L+3) - SIG_EDGE(L+4) ) ) ]
      !  Q(L') = ----------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+4)
      !=================================================================     

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Otherwise,  EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) + &
             ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) ) + &
             ( IN(L+2) * ( EDGE_IN(L+2) - EDGE_IN(L+3) ) ) + &
             ( IN(L+3) * ( EDGE_IN(L+3) - EDGE_IN(L+4) ) ) 

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+4) )
       
      END FUNCTION LUMP_4
!EOC

END MODULE read_input_mod
