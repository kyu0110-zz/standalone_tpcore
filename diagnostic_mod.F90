! ----------------------------------- 
! Write out diagnostics
! Karen Yu 
! 2 August 2015
! ------------------------------------
  MODULE DIAGNOSTIC_MOD
 
! !USES: 
! 
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: WRITE_DIAGNOSTICS
!
! !REMARKS:
!
! !REVISION HISTORY:
!  02 Aug 2015 - K. Yu - Initial version
!EOP
!--------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
  INTEGER, PARAMETER    :: dp = kind(1.0d0)
  INTEGER, PARAMETER    :: sp = kind(1.0) 
!
CONTAINS
!EOC
!--------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: WRITE_DIAGNOSTICS
! 
! !DESCRIPTION:
!
! !INTERFACE:
!
  SUBROUTINE WRITE_DIAGNOSTICS( ncFile, YYYY, MM, DD, h, m, s, wz, wz2, STT )
!
! !USES:
!
  USE CMN_SIZE_MOD
  USE m_netCDF_io_define
  USE Ncdf_Mod,                 ONLY : NC_Create, NC_Close
  USE Ncdf_Mod,                 ONLY : NC_Var_Def, NC_Var_Write
  USE JulDay_Mod,               ONLY : JulDay
  USE Grid_Mod,                 ONLY : get_XMID, get_YMID, get_area_m2
! 
! !INPUT PARAMETERS:
!
  CHARACTER(LEN=*), INTENT(IN) :: ncFile
  INTEGER,          INTENT(IN) :: YYYY
  INTEGER,          INTENT(IN) :: MM
  INTEGER,          INTENT(IN) :: DD
  INTEGER,          INTENT(IN) :: h, m, s
  REAL(dp),         INTENT(IN) :: wz(:,:,:), wz2(:,:,:)
  REAL(dp),         INTENT(IN) :: STT(:,:,:)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! 02 Aug 2015 - K. Yu   - Initial version
!EOP
!--------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
  INTEGER               :: I, J, K
  REAL(dp)              :: GMT, JD1, JD1985, JD_DELTA, THISDAY
  REAL(sp)              :: TMP, JD_DELTA_RND
  REAL(sp), POINTER     :: Arr1D(:)
  INTEGER,  POINTER     :: Int1D(:)
  REAL(sp), POINTER     :: Arr3D(:,:,:)
  REAL(sp), POINTER     :: Arr4D(:,:,:,:)
  INTEGER               :: fId, lonId, latId, levId, TimeId
  INTEGER               :: VarCt
  INTEGER               :: nTime
  INTEGER               :: Prc
  CHARACTER(LEN=31)     :: timeunit, myName, myUnit
  REAL(sp), POINTER     :: nctime(:)
  REAL(dp)              :: XMID(IIPAR)
  REAL(dp)              :: YMID(JJPAR)
  REAL(dp)              :: A_M2(IIPAR,JJPAR)
  REAL(sp), POINTER     :: AREA_M2(:,:) 

  print*, 'ncfile = ', ncfile
  nTime = 1

  ! Get grid properties 
  DO I = 1, IIPAR
    DO J = 1, JJPAR
        A_M2(I,J) = Get_Area_m2(I,J,1)
        YMID(J) = Get_yMid(1,J,1)
    ENDDO
    XMID(I) = Get_xMid(I,1,1)
  ENDDO

  ! Create output file
  CALL NC_CREATE( ncFile, IIPAR,  JJPAR,  LLPAR,  nTime, &
                  fId,    lonId, latId, levId, timeId, VarCt )

  print*, 'timeid', timeid

  ! Add longitude
  CALL NC_VAR_DEF( fId, lonId, -1, -1, -1, &
                   'lon', 'Longitutde', 'degrees_east', 4, VarCt )
  ALLOCATE( Arr1d( IIPAR ) ) 
  Arr1D = XMID
  CALL NC_VAR_WRITE( fId, 'lon', Arr1D=Arr1D )
  DEALLOCATE( Arr1D) 

  ! Add latitude
  CALL NC_VAR_DEF( fId, -1, latId, -1, -1, & 
                   'lat', 'Latitude', 'degrees_north', 4, VarCt )
  ALLOCATE( Arr1D( JJPAR ) ) 
  Arr1D = YMID
  CALL NC_VAR_WRITE( fId, 'lat', Arr1D=Arr1D ) 
  DEALLOCATE( Arr1D )

  ! Add level
  CALL NC_VAR_DEF( fId, -1, levId, -1, -1, &
                   'lev', 'GEOS-Chem level', 'unitless', 1, VarCt ) 
  ALLOCATE( Int1D( LLPAR ) )
  DO I = 1, LLPAR
    Int1D(I) = I
  ENDDO
  CALL NC_VAR_WRITE( fId, 'lev', Arr1D=Int1D ) 
  DEALLOCATe( Int1D ) 

  ! Add time 
  timeunit = 'hours since 1985-01-01 00:00:00 GMT'
  GMT = REAL(h,dp) + (REAL(m,dp)/60.0_dp) + (REAL(s,dp)/3600.0_dp)
  THISDAY  = DD + ( GMT / 24.0_dp )
  JD1      = JULDAY ( YYYY, MM, THISDAY ) 
  JD1985   = JULDAY ( 1985, 1, 0.0_dp ) + 1.0_dp
  JD_DELTA = (JD1 - JD1985) * 24.0_dp

  ! Roun to 2 digits after comma
  JD_DELTA_RND = REAL(JD_DELTA,sp) * 100.0_sp
  TMP          = ANINT( JD_DELTA_RND )
  JD_DELTA_RND = TMP / 100.0_sp
  
  ALLOCATE( nctime(1) ) 
  nctime(1) = JD_DELTA_RND
  CALL NC_VAR_DEF( fId, -1, -1, -1, timeId, & 
                   'time', 'Time', TRIM(timeunit), 4, VarCt ) 
  CALL NC_VAR_WRITE( fId, 'time', Arr1D=nctime )
  DEALLOCATE( nctime ) 

  ! Write out grid box areas
  myName = 'area' 
  myUnit = 'm2' 
  ALLOCATE( AREA_M2(IIPAR, JJPAR) ) 
  CALL NC_VAR_DEF( fId, lonId, latId, -1, -1, & 
                   TRIM(myName), 'Grid box area', TRIM(myUnit), 4, VarCt )
  AREA_M2 = A_M2
  CALL NC_VAR_WRITE( fId, TRIM(myName), Arr2D=AREA_M2 ) 

  ! Write out GEOS omega
  ALLOCATE ( Arr3D(IIPAR, JJPAR, LLPAR) )
  print*, size(wz)
  print*, size(Arr3D)
  Arr3D = wz
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, & 
                   'wzg', 'GEOS-FP omega', 'Pa/s', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'wzg', Arr3D=Arr3D)  
  DEALLOCATE( Arr3D )

  ! Write out tpcore omega
  ALLOCATE( Arr3D(IIPAR, JJPAR, LLPAR) ) 
  Arr3D = wz2
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'wzt', 'TPCORE omega', 'Pa/s', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'wzt', Arr3D=Arr3D ) 
  DEALLOCATE( Arr3D )

  ! Write out time averaged tracers
  ALLOCATE( Arr3D(IIPAR, JJPAR, LLPAR) )
  Arr3D = STT
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'tracer', 'tracer conc', 'v/v', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'tracer', Arr3D=Arr3D )
  DEALLOCATE( Arr3D )

  ! Close file
  CALL NC_CLOSE( fId )

  END SUBROUTINE WRITE_DIAGNOSTICS
END MODULE DIAGNOSTIC_MOD 
