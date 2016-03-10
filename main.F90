!--------------------------------------------------------------------------------
! GEOS-Chem advection-only model
! Runs tpcore indenpendently from GEOS-Chem
!--------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_advection
!
! !DESCRIPTION: Main level driver program for the advection program that runs
! the GEOS-Chem tpcore advection independently from GEOS-Chem 
!
! INTERFACE:
!
    PROGRAM gc_advection
!
! !USES:
!
    USE CMN_SIZE_MOD            ! borrowed from GEOS-Chem 
    USE CMN_GCTM_MOD            ! borrowed from GEOS-Chem
    USE ERROR_MOD               ! borrowed from GEOS-Chem
    USE read_input_mod,        ONLY : GeosFp_Read_A3dyn, GeosFp_Read_I3
    USE read_input_mod,        ONLY : restart_read, geosfp_read_omega 
    USE read_input_mod,        ONLY : read_wc, read_Kz
    USE TPCORE_FVDAS_MOD,       ONLY : Tpcore_FvDas 
    USE PRESSURE_MOD
    USE PJC_PFIX_MOD
    USE GRID_MOD
    USE INITIALIZE_MOD
    USE ADVECTION_MOD
    USE CHEM_MOD,               ONLY: DO_EMISS, DO_DECAY
    USE DIAGNOSTIC_MOD
    USE TIME_MOD,               ONLY : SET_BEGIN_TIME, SET_END_TIME
    USE TIME_MOD,               ONLY : SET_TIMESTEPS, PRINT_CURRENT_TIME
    USE TIME_MOD,               ONLY : SET_CURRENT_TIME, SET_CT_DYN
    USE TIME_MOD,               ONLY : SET_ELAPSED_MIN, ITS_TIME_FOR_EXIT
    USE TIME_MOD,               ONLY : ITS_TIME_FOR_MET, GET_NYMD
    USE TIME_MOD,               ONLY : GET_NHMS, GET_TIME_AHEAD, GET_CT_DYN
    USE TIME_MOD,               ONLY : GET_HOUR, GET_MINUTE, GET_SECOND

    IMPLICIT NONE

! !REVISION HISTORY:
! 24 Jul 2015 - K. Yu     - Initial version
!
!EOP
!--------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER        :: IM, JM, KM          ! counters for loops            
    INTEGER        :: RC                  ! return int for program success 
    INTEGER        :: I, J, L, T
    INTEGER        :: DDATE                 ! diagnostics variables
    INTEGER        :: NUM_ARGS
    CHARACTER(LEN=255), ALLOCATABLE    :: ARGS(:)
    INTEGER        :: N_ARGS
    INTEGER        :: start_time, end_time, start_date, end_date
 
    REAL*8         :: DT                    ! transport timestep [s]

    LOGICAL        :: LFILL = .TRUE.
    LOGICAL        :: FIRST = .TRUE.      ! used for initializations

    CHARACTER(LEN=255)      :: RESTART = 'F'   ! read restart file?
    CHARACTER(LEN=255)    :: fname
    CHARACTER(LEN=255)    :: ffname

    ! Arrays
    INTEGER               :: A3_TIME(2)           ! array for u,v,w time
    INTEGER               :: I3_TIME(2)           ! array for ps time
    REAL*8, ALLOCATABLE   :: wzg(:,:,:)           ! geosfp omega
    REAL*8, ALLOCATABLE   :: wzt(:,:,:)           ! tpcore omega
    REAL*8, ALLOCATABLE   :: wzgneg(:,:,:)         ! negative omega
    REAL*8, ALLOCATABLE, TARGET   :: U(:,:,:)     ! u wind
    REAL*8, ALLOCATABLE, TARGET   :: V(:,:,:)     ! v wind
    REAL*8, ALLOCATABLE, TARGET   :: t_inst(:,:,:)   ! tracer concentrations
    REAL*8, ALLOCATABLE   :: PS(:,:)              !surface pressure
    REAL*8, ALLOCATABLE   :: PS2(:,:)              !surface pressure
    REAL*8, ALLOCATABLE   :: PFLT(:,:)            ! floating surface pressure
    REAL*8                :: TCVV                 ! tracer mw
    REAL*8, ALLOCATABLE   :: t_avg(:,:,:)    ! time avg tracer conc
    REAL*8, ALLOCATABLE   :: wzg_AVG(:,:,:)          ! time avg omega
    REAL*8, ALLOCATABLE   :: wzt_AVG(:,:,:)         ! time avg tpcore omega
    REAL*8, ALLOCATABLE   :: TRACERFLUX(:,:,:)      ! mass flux of tracer
    REAL*8, ALLOCATABLE   :: Kz(:,:,:)          ! Kz
    REAL*8, ALLOCATABLE   :: c_bar(:,:,:)          ! Kz
    REAL*8, ALLOCATABLE   :: wc_prime(:,:,:)
    REAL*8, ALLOCATABLE   :: dcdp(:,:,:)

    !============================================================================
    ! Code begins here
    !============================================================================

    print*, 'Starting advection model'
 
    !============================================================================
    ! Set up grid and time variables
    !============================================================================

    ! Initialize grid 
    CALL INIT_CMN_SIZE( RC )

    ! Allocate arrays
    ALLOCATE( wzg(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( wzt(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( U(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( V(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( t_avg(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( t_inst(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( PS(IIPAR, JJPAR) )
    ALLOCATE( PS2(IIPAR, JJPAR) )
    ALLOCATE( PFLT(IIPAR, JJPAR) )
    ALLOCATE( wzt_AVG(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( wzg_AVG(IIPAR, JJPAR, LLPAR ) )
#if defined( GRID4x5 ) || defined( GRID2x25 )
    ALLOCATE( wzgneg(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( tracerflux(IIPAR, JJPAR, LLPAR))
    ALLOCATE( Kz(IIPAR, JJPAR, LLPAR) )
    ALLOCATE( wc_prime(IIPAR, JJPAR, LLPAR) )
#else 
    ALLOCATE( Kz(72, 46, LLPAR))
    ALLOCATE( c_bar(72, 46, LLPAR))
    ALLOCATE( dcdp(72, 46, LLPAR) )
    ALLOCATE( wc_prime(72, 46, LLPAR) )
#endif

    ! Select timestep based on grid
#if defined( GRID4x5 )  
    DT = 1800d0
    !DT = 300d0
#elif defined( GRID2x25 ) 
    DT = 300d0
    !DT = 900d0
#elif defined( GRID025x03125 ) 
    DT = 300d0
#endif

    NUM_ARGS = COMMAND_ARGUMENT_COUNT()
    ALLOCATE(ARGS(NUM_ARGS))
    DO N_ARGS=1, NUM_ARGS
        CALL GET_COMMAND_ARGUMENT(N_ARGS, ARGS(N_ARGS))
        print*, args(n_args)
    ENDDO

    ! Set up time
    read(ARGS(1), '(i10)') start_date
    read(ARGS(2), '(i10)') start_time
    read(ARGS(3), '(i10)') end_date
    read(ARGS(4), '(i10)') end_time
    CALL SET_BEGIN_TIME( start_date, start_time )
    CALL SET_END_TIME( end_date, end_time )
    CALL SET_TIMESTEPS( INT(DT/60) )
    CALL SET_CURRENT_TIME
    CALL SET_CT_DYN

    ! TCVV is MW air / MW tracer. Set to 1 for tracer that has the same MW as air
    !TCVV = 1.0
    TCVV = 28.97 / 222.0d0

    ! Set up grid
    CALL SETUP_GRID

    ! Read restart file or initialize tracer array with 1 ppt mixing 
    ! ratio tracer on top and bottom
    ! rows
    restart = TRIM(ARGS(5))
    IF (restart == 'F') THEN
      t_inst = 0        ! start with zero concentrations
      !t_inst(:,:,1) = 0.001      ! start with zero concentrations
      wzt = 0
    ELSE 
      t_inst = 0 
      call restart_read(restart, t_inst)
    ENDIF
      t_avg = t_inst

    Kz = 0

    ! Read in surface pressure before starting loop
    CALL GeosFp_Read_I3( GET_NYMD(), GET_NHMS(), PS2 )

    !============================================================================
    ! Loop over timesteps
    !============================================================================
    DO 

      IF (.NOT. FIRST) CALL SET_ELAPSED_MIN
      CALL SET_CURRENT_TIME
      CALL PRINT_CURRENT_TIME

      !==========================================================================
      ! Read in met if appropriate
      !==========================================================================
      IF ( ITS_TIME_FOR_MET() ) THEN
        ! Set time variables
        A3_TIME = GET_TIME_AHEAD(90)
        I3_TIME= GET_TIME_AHEAD(180)

        ! Set previous next pressure to current presure
        PS = PS2

        ! Set floating pressure to current pressure
        PFLT = PS

        ! Read in next surface pressure
        CALL GeosFp_Read_I3( I3_TIME(1), I3_TIME(2), PS2 )

        ! Read in omega, u, and v
        CALL GeosFp_Read_A3dyn( A3_TIME(1), A3_TIME(2), wzg, U, V )
     
#if defined(GRID4x5) || defined(GRID2x25)
        ! Read in omegas 
        CALL geosfp_read_omega(A3_TIME(1), A3_TIME(2), wzgneg)
        
#endif

     ENDIF

      !============================================================================
      ! Initialize transport if haven't yet
      !============================================================================
      IF (FIRST) CALL INIT_ADVECTION(DT, RC)
      IF (FIRST) wzt = wzg
      FIRST = .FALSE.

      !============================================================================
      ! Do emissions
      !============================================================================
      CALL DO_EMISS(DT, t_inst, PFLT)
    print*, 'after emission', sum(t_inst)

#if defined(GRID4x5) || defined(GRID2x25)
      ! compute flux before doign transport
      !CALL COMPUTE_FLUX(DT, wzg, wzgneg, t_inst, PS2, TCVV, TRACERFLUX)
#endif

      !============================================================================
      ! Do transport
      !============================================================================

#if defined(GRID025x03125)
        ! write out Kz
        CALL COMPUTE_KZ(wzt, t_inst, Kz, c_bar, PFLT, 100.0d0, dcdp, wc_prime)

        IF (ITS_TIME_FOR_MET()) THEN

        DO L = 1, LLPAR
        DO J = 1, 46
        DO I = 1, 72
        Kz(I,J,L) = Kz(I,J,L) / 36.0d0
        ENDDO
        ENDDO
        ENDDO

        Kz = 0
        ENDIF

        DDATE = GET_NYMD()
        WRITE(fname, '(i8)') DDATE
        DDATE = GET_NHMS()
        WRITE(ffname, '(i6.6)') DDATE
        fname = TRIM(fname) // TRIM(ffname)
        fname = TRIM(fname) //'.nc'
        fname = '/n/regal/jacob_lab/kyu/wc_prime_withdcdp/' // TRIM(fname)
        CALL WRITE_WCprime(TRIM(fname), DDATE/10000, &
                MOD(DDATE, 10000)/100, MOD(DDATE, 100), GET_HOUR(), &
                GET_MINUTE(), GET_SECOND(), wc_prime, dcdp)
#else
        IF (ITS_TIME_FOR_MET()) THEN 
        ! read in Kz
        print*, GET_NHMS()
        CALL READ_KZ(GET_NYMD(), GET_NHMS(), Kz)
        ENDIF
      !CALL COMPUTE_FLUX(DT, Kz, wzt, t_inst, PFLT, TCVV, TRACERFLUX)
#endif 

      CALL DO_ADVECTION(LFILL, DT, PFLT, PS2, U, V, t_inst, TCVV, wzt, RC)
    print*, 'after advection', sum(t_inst)


      ! ==========================================================================
      ! Do flux echange
      ! ==========================================================================
#if defined(GRID4x5) || defined(GRID2x25)
        print*, 'sum wc 4x5', sum((wzt / 9.81 ) * (t_inst / TCVV) )

      CALL COMPUTE_FLUX(DT, Kz, wzt, t_inst, PFLT, TCVV, TRACERFLUX, wc_prime)
      CALL DO_WC_TRANSPORT(DT, t_inst, PFLT, TCVV, TRACERFLUX )

        DDATE = GET_NYMD()
        WRITE(fname, '(i8)') DDATE
        DDATE = GET_NHMS()
        WRITE(ffname, '(i6.6)') DDATE
        fname = TRIM(fname) // TRIM(ffname)
        fname = TRIM(fname) //'.nc'
        fname = '/n/regal/jacob_lab/kyu/wc_prime/4x5/' // TRIM(fname)

        CALL WRITE_WC(TRIM(fname), DDATE/10000, &
                MOD(DDATE, 10000)/100, MOD(DDATE, 100), GET_HOUR(), &
                GET_MINUTE(), GET_SECOND(), wc_prime)
#endif

      !===========================================================================
      ! Do radioactive decay of species
      !============================================================================

      CALL DO_DECAY(DT, t_inst)
    print*, 'after decay', sum(t_inst)
 
      !============================================================================
      ! Write diagnostics if end of the day
      !============================================================================
      ! Save to diagnostics array
      t_avg = t_avg + t_inst
      wzg_AVG = wzg_AVG + wzg
      wzt_AVG = wzt_AVG + wzt

      IF (GET_HOUR() == 23 .AND. GET_MINUTE() == 0) THEN
        DDATE = GET_NYMD()
        WRITE(fname, '(i8)') DDATE
        fname = TRIM(fname)//'.nc'
        CALL WRITE_DIAGNOSTICS(TRIM(fname), DDATE/10000, &
                MOD(DDATE, 10000)/100, MOD(DDATE, 100), GET_HOUR(), & 
                GET_MINUTE(), GET_SECOND(), wzg, wzt, t_inst )
      ENDIF

      !============================================================================
      ! Increment counters
      !============================================================================
      CALL SET_CT_DYN( INCREMENT=.TRUE. )

      IF (ITS_TIME_FOR_EXIT() ) EXIT
    ENDDO

    !==============================================================================
    ! Write out diagnostics
    !==============================================================================
    t_avg = t_avg / GET_CT_DYN()
    wzg_AVG = wzg_AVG / GET_CT_DYN()
    wzt_AVG = wzt_AVG / GET_CT_DYN()
    DDATE = GET_NYMD()
    CALL WRITE_DIAGNOSTICS('average.nc', DDATE/10000, &
                MOD(DDATE, 10000)/100, MOD(DDATE, 100), GET_HOUR(), &
                GET_MINUTE(), GET_SECOND(), wzg_AVG, wzt_AVG, t_avg)

    write( 6, '(a)' ) repeat( '=', 79 ) 
    write( 6, '(a)' ) 'Simulation completed'
    write( 6, '(a)' ) repeat( '=', 79 )
 
    END PROGRAM GC_ADVECTION 
