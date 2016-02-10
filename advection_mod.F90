!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transport_mod
!
! !DESCRIPTION: Module TRANSPORT\_MOD is used to call the proper version of 
!  the TPCORE advection scheme for GCAP, GEOS-4, GEOS-5, or GEOS-5.7
!  nested-grid or global simulations.
!\\
!\\
! !INTERFACE: 
!
      MODULE ADVECTION_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: INIT_ADVECTION
      PUBLIC  :: DO_ADVECTION
      PUBLIC  :: COMPUTE_FLUX
      PUBLIC  :: DO_FLUX_EXCHANGE
      PUBLIC  :: DO_WC_TRANSPORT
      PUBLIC  :: HOMOGENIZE
      PUBLIC  :: GET_AIR_MASS
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! (1 ) Ap     (REAL*8 ) : Vertical coordinate array for TPCORE
      ! (2 ) A_M2   (REAL*8 ) : Grid box surface areas [m2]
      ! (3 ) Bp     (REAL*8 ) : Vertical coordinate array for TPCORE
      ! (4 ) IORD   (REAL*8 ) : TPCORE E/W option flag
      ! (5 ) JORD   (REAL*8 ) : TPCORE N/S option flag
      ! (6 ) KORD   (REAL*8 ) : TPCORE vertical option flag
      ! (7 ) JLAST  (INTEGER) : For fvDAS TPCORE
      ! (8 ) MG     (INTEGER) : For fvDAS TPCORE
      ! (9 ) NG     (INTEGER) : For fvDAS TPCORE
      ! (10) N_ADJ  (INTEGER) : For fvDAS TPCORE
      !=================================================================
      INTEGER             :: IORD,  JORD, KORD, JFIRST 
      INTEGER             :: JLAST, NG,   MG,   N_ADJ
      REAL*8, ALLOCATABLE, TARGET :: Ap(:)
      REAL*8, ALLOCATABLE :: A_M2(:)
      REAL*8, ALLOCATABLE, TARGET :: Bp(:)
      REAL*8, ALLOCATABLE :: Ap_fullgrid(:)
      REAL*8, ALLOCATABLE :: Bp_fullgrid(:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_transport
!
! !DESCRIPTION: Subroutine INIT\_TRANSPORT initializes all module variables 
!  and arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ADVECTION( DT, RC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_YMID_R
      USE PRESSURE_MOD
      USE TPCORE_FVDAS_MOD,   ONLY : INIT_TPCORE

      USE CMN_SIZE_MOD         ! Size parameters
      USE CMN_GCTM_MOD         ! Re
!
! !INPUT PARAMETERS:
!
      REAL*8,         INTENT(IN) :: DT
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J, K, L, N_DYN
      REAL*8  :: YMID_R(JJPAR)

      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates 
      !
      ! For TPCORE v7.1.m (for GEOS-3 met fields):
      ! 
      !    P(I,J,L) = ( Ap(L) * PTOP ) + ( Bp(L) * ( Psurf(I,J)-PTOP ) )
      !
      ! For fvDAS TPCORE with for GEOS-4 or GEOS-5 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )
 
      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ALLOCATE( Ap_fullgrid(LGLOB+1), STAT=RC)
      ALLOCATE( Bp_fullgrid(LGLOB+1), STAT=RC)

      CALL INIT_PRESSURE(Ap, Bp, Ap_fullgrid, Bp_fullgrid)

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the GEOS-4/fvDAS version of TPCORE
      !=================================================================

      ! Initialize
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! YMID_R is latitude of grid box center [radians]
      DO J = 1,JJPAR
         YMID_R(J) = GET_YMID_R( 1, J, 1 )
      ENDDO

      ! Call INIT routine from "tpcore_fvdas_mod.f" 
      CALL INIT_TPCORE( IIPAR,  JJPAR, LLPAR,  JFIRST, JLAST, & 
                     NG, MG,  DT, Re,     YMID_R )

      ! Set IORD, JORD, KORD
      IORD = 3
      JORD = 3
      KORD = 7

      END SUBROUTINE INIT_ADVECTION
!EOC
!BOC
! !INTERFACE:
      SUBROUTINE DO_ADVECTION(LFILL, DT, PS, PS2, U, V, STT, &
                              TCVV, wz, RC)
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
        USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
        USE ERROR_MOD
! !INPUT PRAMETERS: 
! 
        LOGICAL,        INTENT(IN)      :: LFILL
        REAL*8,         INTENT(IN)      :: DT
        REAL*8,         INTENT(INOUT)      :: PS(IIPAR,JJPAR)
        REAL*8,         INTENT(IN)      :: PS2(IIPAR,JJPAR)
        REAL*8, TARGET, INTENT(IN)      :: U(IIPAR,JJPAR,LLPAR)
        REAL*8, TARGET, INTENT(IN)      :: V(IIPAR,JJPAR,LLPAR)
        REAL*8, TARGET, INTENT(INOUT)   :: STT(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)      :: TCVV
! !OUTPUT PARAMETERS:
!
        INTEGER,        INTENT(OUT)     :: RC
        REAL*8,         INTENT(OUT)     :: wz(IIPAR,JJPAR,LLPAR)
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        REAL*8          :: P_TP1(IIPAR,JJPAR)
        REAL*8          :: P_TP2(IIPAR,JJPAR)
        REAL*8          :: P_TEMP(IIPAR,JJPAR)
        REAL*8, TARGET  :: XMASS(IIPAR, JJPAR, LLPAR)
        REAL*8, TARGET  :: YMASS(IIPAR, JJPAR, LLPAR)

        REAL*8, POINTER :: p_U(:,:,:)
        REAL*8, POINTER :: p_V(:,:,:)
        REAL*8, POINTER :: p_XMASS(:,:,:)
        REAL*8, POINTER :: p_YMASS(:,:,:)
        REAL*8, POINTER :: p_Ap(:)
        REAL*8, POINTER :: p_Bp(:)
        REAL*8, POINTER :: p_STT(:,:,:)

        ! Variable for mass conservation
        REAL*8          :: SUMADA
        REAL*8          :: TR_DIFF, A_DIFF
      
        REAL*8          :: AD_A
        REAL*8          :: AD_B 
        REAL*8          :: TR_A(IIPAR, JJPAR, LLPAR)
        REAL*8          :: TR_B(IIPAR, JJPAR, LLPAR)
         

        print*, 'at start of do_advection', sum(STT)
        ! Get surface pressures
        DO J = 1, JJPAR
        DO I = 1, IIPAR
          ! True surface pressure at midpoint of dynamic timestep [hPa]
          P_TP1(I,J) = PS(I,J)
          ! True surface pressure at end of dynamic timestep [hPa]
          P_TP2(I,J) = PS2(I,J)
        ENDDO
        ENDDO

        !==========================================================
        ! Get air and tracer mass before advection
        !==========================================================

        ! Tracer mass before transport
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE(I,J,L,AD_B) 
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
          AD_B = GET_AIR_MASS(I, J, L, P_TP1(I,J))
          TR_B(I,J,L) = STT(I,J,L) * AD_B / TCVV
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        print*, 'trb', SUM(TR_B(:,:,:))
        print*, 'p_tp1', 'p_tp2', sum(p_tp1)/(IIPAR*JJPAR*LLPAR), sum(p_tp2)/(IIPAR*JJPAR*LLPAR)
        !===========================================================
        ! Pressure fixer to compute xmass and ymass
        !===========================================================
        CALL DO_PJC_PFIX( dt, P_TP1+PTOP, P_TP2+PTOP, U, V, XMASS, &
                          YMASS, Ap, Bp )
        print*, 'xmass, ymass', SUM(XMASS), SUM(YMASS)

        !===========================================================
        ! TPCORE 
        !===========================================================

        ! Flip arrays for tpcore
        p_U => U(:,:,LLPAR:1:-1)
        p_V => V(:,:,LLPAR:1:-1)
        p_Ap => Ap(LLPAR+1:1:-1)
        p_Bp => Bp(LLPAR+1:1:-1)
        p_XMASS => XMASS(:,:,LLPAR:1:-1)
        p_YMASS => YMASS(:,:,LLPAR:1:-1)
        p_STT => STT(:,:,LLPAR:1:-1)

        print*, 'before tpcore', sum(STT)
        ! Tpcore for transport
        CALL TPCORE_FVDAS( dt, Re, IIPAR, JJPAR, LLPAR, JFIRST, JLAST, &
                           0, 0, 1, p_Ap, p_Bp, p_U, p_V, P_TP1, &
                           P_TP2,  P_TEMP, p_STT, IORD, JORD, KORD, 0, &
                           p_XMASS, p_YMASS, LFILL, A_M2, TCVV, wz )

        print*, 'after tpcore', sum(STT)
        ! Flip wz
        wz = wz(:,:,LLPAR:1:-1)
      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      PS = P_TP2

      ! Adjust tracer to correct residual non-conservation of mass
      ! Zero summing variable
       SUMADA = 0.d0

!$OMP PARALLEL DO           &
!$OMP DEFAULT( SHARED )     &
!$OMP PRIVATE( I, J, L, AD_A, SUMADA )    
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR

          ! Air mass [kg] after transport
          AD_A = GET_AIR_MASS( I, J, L, P_TP2(I,J))
         
          ! Tracer mass [kg] after transport
          TR_A(I,J,L) = STT(I,J,L) * AD_A / TCVV

          IF ( STT(I,J,L) > 0.d0   .or. STT(I,J,L) < 0.d0 ) THEN
               SUMADA = SUMADA + AD_A
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

       ! Residual mass difference [kg]: before - after
       TR_DIFF = SUM( TR_B ) - SUM( TR_A )

       ! Convert from [kg] to [v/v]
       TR_DIFF = SAFE_DIV(TR_DIFF, SUMADA, 0.d0) * TCVV

!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L ) 
       ! Add mass difference [v/v] back to STT
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR

          IF ( STT(I,J,L) > 0.d0   .or. STT(I,J,L) < 0.d0 ) THEN
               STT(I,J,L) = STT(I,J,L) + TR_DIFF
          ENDIF

          STT(I,J,L) = MAX( STT(I,J,L), 0d0 )

       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

        print*, 'after residual correction', sum(STT)         
      END SUBROUTINE DO_ADVECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOC
! !INTERFACE:
      SUBROUTINE COMPUTE_FLUX( DT, OMEGA, NEGOMEGA, &
                              STT, PS, TCVV, TRACERFLUX)
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE ERROR_MOD
        USE GRID_MOD
! !INPUT PRAMETERS: 
! 
        REAL*8,         INTENT(IN)      :: DT
        REAL*8,         INTENT(IN)   :: OMEGA(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)   :: NEGOMEGA(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(INOUT)   :: STT(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)      :: TCVV
        REAL*8,         INTENT(IN)      :: PS(:,:)
        REAL*8,         INTENT(OUT)   :: TRACERFLUX(IIPAR, JJPAR, LLPAR)
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        REAL*8          :: ZMASS(IIPAR, JJPAR, LLPAR)
        REAL*8          :: AREA_M2(IIPAR,JJPAR)
        REAL*8          :: AD
        REAL*8          :: BOXMASS

! START EXECUTION
        ! Print info
        WRITE( 6, '(a)' ) 'Compute fluxes'

        ! Compute surface area of grid box
!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J )   
        DO I = 1, IIPAR
          DO J = 1, JJPAR
            AREA_M2(I,J) = GET_AREA_M2(I,J,1)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO      

        ! up/down mass flux [kg air / m^2]
        ! If the net omega is in the same direction as negomega, the
        ! additional flux ist he difference between the two. If the next
        ! flux is in the opposite direction, the additional flux is just
        ! negomega.
        WHERE ( OMEGA .lt. 0d0 ) & 
         ZMASS = ((NEGOMEGA / -9.81) * DT) - ((OMEGA / -9.81) * DT)
        WHERE ( OMEGA .ge. 0d0 ) & 
         ZMASS = ( NEGOMEGA / -9.81) * DT

        ! UP/DOWN mass flux [kg air]
        DO L = 1, LLPAR
            ZMASS(:,:,L) = ZMASS(:,:,L) * AREA_M2(:,:)
        ENDDO

        print*, 'sum ZMASS', sum(ZMASS)
            ! Additional mass flux of each tracer is the mass fluxo f
            ! air multiplied by the mixign ratio of the tracer divided
            ! by the ratio of molec w eight of tracer to the molec
            ! weight of air
            ! [kg tracer] = [kg air] * [mol tracer/mol air] * 
            !               [kg tracer/mol tracer * mol air/kg air]
            TRACERFLUX(:,:,:) = ZMASS(:,:,:) * STT(:,:,:) / TCVV

!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L, AD, BOXMASS)
        ! tracerflux cannot be more than half the mass that's currently
        ! in the box
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
            AD = GET_AIR_MASS(I, J, L, PS(I,J))
            BOXMASS = 0.5 * (STT(I,J,L) * AD / TCVV)
            TRACERFLUX(I,J,L) = min(TRACERFLUX(I,J,L), BOXMASS)
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        print*, 'sum tracerflux', sum(tracerflux)
        print*, 'max tracerflux', maxval(TRACERFLUX)
       END SUBROUTINE COMPUTE_FLUX 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOC
! !INTERFACE:
      SUBROUTINE DO_FLUX_EXCHANGE( DT, STT, PS, TCVV, TRACERFLUX)
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE ERROR_MOD
        USE GRID_MOD
! !INPUT PRAMETERS: 
! 
        REAL*8,         INTENT(IN)      :: DT
        REAL*8,         INTENT(INOUT)   :: STT(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)      :: TCVV
        REAL*8,         INTENT(IN)      :: PS(:,:)
        REAL*8,         INTENT(IN)      :: TRACERFLUX(IIPAR,JJPAR,LLPAR)
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        REAL*8          :: Q_NEW(IIPAR, JJPAR, LLPAR)
        REAL*8          :: AD(IIPAR, JJPAR, LLPAR)

! START EXECUTION
        ! Print info
        WRITE( 6, '(a)' ) 'Performing vertical flux exchange'

        print*, 'at start of flux ex', sum(STT)

!$OMP PARALLEL DO           &
!$OMP DEFAULT( SHARED )     &
!$OMP PRIVATE( I, J, L ) 
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
            AD(I,J,L) = GET_AIR_MASS(I, J, L, PS(I,J))
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        ! Bottom grid box first: new mass = mass currently in cell
        ! minus mass out + mass in from box above
        Q_NEW(:,:,1) = (STT(:,:,1) * AD(:,:,1) / TCVV) - TRACERFLUX(:,:,1) + &
                           TRACERFLUX(:,:,2)

        print*, 'maxval qnew', MAXval(Q_NEW(:,:,1))
        DO L = 2, LLPAR-1

            ! New amass = mass currently in cell - 2 * mass out (one
            ! for up and one for down) + mass in from cell above +
            ! mass in from cell below
            Q_NEW(:,:,L) = (STT(:,:,L) * AD(:,:,L)) / &
                            TCVV - 2*TRACERFLUX(:,:,L) + &
                            TRACERFLUX(:,:,L-1) + TRACERFLUX(:,:,L+1)
        ENDDO

        print*, 'maxval qnew', MAXval(Q_NEW)
        ! Top layer
        ! New mass = mass currently in cell - mass out (down) + mass
        ! in from cell below
         Q_NEW(:,:,47) = STT(:,:,47) * AD(:,:,47) / &
                        TCVV - TRACERFLUX(:,:,47) + &
                        TRACERFLUX(:,:,46)

        print*, 'maxval qnew', MAXval(Q_NEW)

!$OMP PARALLEL DO           &
!$OMP DEFAULT ( SHARED )    &
!$OMP PRIVATE( I, J, L )    
        ! Reset original array with the new values, converting from
        ! [kg] to [v/v]
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
        STT(I,J,L) = Q_NEW(I,J,L) * TCVV / AD(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        print*, 'maxval stt', maxval(stt)

        print*, 'at end of flux ex stt', sum(STT(:,:,:))
            
       END SUBROUTINE DO_FLUX_EXCHANGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOC
! !INTERFACE:
      SUBROUTINE DO_WC_TRANSPORT( DT, STT, PS, TCVV, TRACERFLUX)
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE ERROR_MOD
        USE GRID_MOD
! !INPUT PRAMETERS: 
! 
        REAL*8,         INTENT(IN)      :: DT
        REAL*8,         INTENT(INOUT)   :: STT(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)      :: TCVV
        REAL*8,         INTENT(IN)      :: PS(:,:)
        REAL*8,         INTENT(IN)      :: TRACERFLUX(IIPAR,JJPAR,LLPAR)
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        REAL*8          :: Q_NEW(IIPAR, JJPAR, LLPAR)
        REAL*8          :: AD(IIPAR, JJPAR, LLPAR)
        REAL*8          :: T

! START EXECUTION
        ! Print info
        WRITE( 6, '(a)' ) 'Performing vertical flux exchange'

        print*, 'at start of flux ex', sum(STT)

!$OMP PARALLEL DO           &
!$OMP DEFAULT( SHARED )     &
!$OMP PRIVATE( I, J, L ) 
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
            AD(I,J,L) = GET_AIR_MASS(I, J, L, PS(I,J))
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        ! Do bottom and top first
        Q_NEW(:,:,:) = (STT(:,:,:) * AD(:,:,:) / TCVV)

        print*, 'maxval qnew', MAXval(Q_NEW(:,:,1))

        DO L = 2, LLPAR-1
        DO I = 1, IIPAR
        DO J = 1, JJPAR

           !Q_NEW(I,J,L) = (STT(I,J,L) * AD(I,J,L)) / &
                            !TCVV !- TRACERFLUX(I,J,L)

            !T = Q_NEW(I,J,L) * 0.00001
            IF (TRACERFLUX(I,J,L) > 0) THEN
            Q_NEW(I,J,L) = Q_NEW(I,J,L) - TRACERFLUX(I,J,L)
              Q_NEW(I,J,L-1) = Q_NEW(I,J,L-1) + &
                            TRACERFLUX(I,J,L)
            ELSE 
                Q_NEW(I,J,L) = Q_NEW(I,J,L) + TRACERFLUX(I,J,L)
                Q_NEW(I,J,L+1) = Q_NEW(I,J,L+1) - TRACERFLUX(I,J,L)
            ENDIF
        ENDDO
        ENDDO
        ENDDO

        print*, 'maxval qnew', MAXval(Q_NEW)

!$OMP PARALLEL DO           &
!$OMP DEFAULT ( SHARED )    &
!$OMP PRIVATE( I, J, L )    
        ! Reset original array with the new values, converting from
        ! [kg] to [v/v]
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
        STT(I,J,L) = Q_NEW(I,J,L) * TCVV / AD(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

        print*, 'maxval stt', maxval(stt)

        print*, 'at end of flux ex stt', sum(STT(:,:,:))
            
       END SUBROUTINE DO_WC_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOC
! !INTERFACE:
      SUBROUTINE HOMOGENIZE( STT )
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE ERROR_MOD
        USE GRID_MOD
        USE REGRID_A2A_MOD, ONLY : MAP_A2A
! !INPUT PRAMETERS: 
! 
        REAL*8,TARGET, INTENT(INOUT)   :: STT(IIPAR,JJPAR,LLPAR)
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        REAL*8          :: LON_025(IIPAR+1)
        REAL*8          :: LON_4x5(73)
        REAL*8          :: LAT_025(JJPAR+1)
        REAL*8          :: LAT_4x5(47)
        REAL*8, POINTER :: Q_025(:, :)
        REAL*8          :: Q_4x5(72, 46)

! START EXECUTION
        ! Print info
        WRITE( 6, '(a)' ) 'Horizontal homogenization'

        ! Compute lat and lons
        DO I = 1, IIPAR+1
            LON_025(I) = GET_XEDGE(I, 1, 1)
        ENDDO
 
        DO J = 1, JJPAR+1
            LAT_025(J) = GET_YSIN(1, J, 1)
        ENDDO

        ! Read lats and lons from file
        LON_4x5 = (/ -182.5, -177.5, -172.5, -167.5, -162.5, -157.5, -152.5, -147.5, -142.5, -137.5, -132.5, -127.5, -122.5, -117.5, -112.5, -107.5, -102.5, -97.5, -92.5, -87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5, 102.5, 107.5, 112.5, 117.5, 122.5, 127.5, 132.5, 137.5, 142.5, 147.5, 152.5, 157.5, 162.5, 167.5, 172.5, 177.5 /)

        LAT_4x5 = (/ -1.0, -0.999390827, -0.9945218954, -0.984807753, -0.9702957263, -0.9510565163, -0.9271838546, -0.8987940463, -0.8660254038, -0.8290375726, -0.7880107536, -0.7431448255, -0.6946583705, -0.6427876097, -0.5877852523, -0.5299192642, -0.4694715628, -0.4067366431, -0.3420201433, -0.2756373558, -0.2079116908, -0.139173101, -0.0697564737, 0.0, 0.0697564737, 0.139173101, 0.2079116908, 0.2756373558, 0.3420201433, 0.4067366431, 0.4694715628, 0.5299192642, 0.5877852523, 0.6427876097, 0.6946583705, 0.7431448255, 0.7880107536, 0.8290375726, 0.8660254038, 0.8987940463, 0.9271838546, 0.9510565163, 0.9702957263, 0.984807753, 0.9945218954, 0.999390827, 1.0 /)
 
      DO L=1, LLPAR
        Q_025 => STT(:, :, L)

        ! Regrid input field from 0.25x0.3125 to 4x5 grid
        CALL MAP_A2A(IIPAR, JJPAR, LON_025, LAT_025, Q_025, & 
               72, 46, LON_4x5, LAT_4x5, Q_4x5, 0, 0 )  

        ! Regrid from 4x5 back to 0.25x0.3125
        CALL MAP_A2A(72, 46, LON_4x5, LAT_4x5, Q_4x5,   &
               IIPAR, JJPAR, LON_025, LAT_025, Q_025, 0, 0 )

      ENDDO

        print*, 'maxval stt', maxval(stt)

        print*, 'at end of homogenization', sum(STT)
            
       END SUBROUTINE HOMOGENIZE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_air_mass
!
! !DESCRIPTION: Function GET\_AIR\_MASS returns the air mass based on the 
!  pressures returned before and after the call to the GEOS-4/fvDAS TPCORE 
!  code. (bmy, 6/24/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_AIR_MASS( I, J, L, P_SURF ) RESULT( AIR_MASS )
!
! !USES:
!
      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_GCTM_MOD               ! g0_100
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I, J, L   ! GEOS-Chem lon, lat, level indices
      REAL*8,  INTENT(IN) :: P_SURF    ! Surface pressure [hPa] at (I,J,L=1)

! 
! !REVISION HISTORY: 
!  24 Jun 2003 - R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L2
      REAL*8  :: P_BOT, P_TOP, AIR_MASS
      
      !=================================================================
      ! GET_AIR_MASS begins here!
      !=================================================================
      
      ! Index for Ap, Bp from atmosphere top down to surface
      ! since the Ap's and Bp's have been flipped for TPCORE
!      L2       = ( LLPAR + 1 ) - L + 1
              
      ! Hybrid-grid formulation for air mass
      P_BOT    = Ap(L)   + ( Bp(L)   * P_SURF )
      !P_BOT    = Ap(L2)   + ( Bp(L2)   * P_SURF )
      P_TOP    = Ap(L+1) + ( Bp(L+1) * P_SURF )
      !P_TOP    = Ap(L2-1) + ( Bp(L2-1) * P_SURF )
      AIR_MASS = ( P_BOT - P_TOP ) * G0_100 * A_M2(J)

      END FUNCTION GET_AIR_MASS
!EOC

      END MODULE ADVECTION_MOD
