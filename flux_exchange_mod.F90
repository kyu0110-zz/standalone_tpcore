!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transport_mod
!
! !DESCRIPTION: Module for routines related to flux change method for
! vertical transport
!\\
!\\
! !INTERFACE: 
!
      MODULE FLUX_EXCHANGE_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_FLUX_EXCHANGE
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOC
! !INTERFACE:
      SUBROUTINE DO_FLUX_EXCHANGE( N_TRC, DT, OMEGA, NEGOMEGA, &
                              STT, AD, TCVV)
! !USES:
        USE CMN_GCTM_MOD
        USE CMN_SIZE_MOD
        USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
        USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
        USE ERROR_MOD
        USE GRID_MOD
! !INPUT PRAMETERS: 
! 
        INTEGER,        INTENT(IN)      :: N_TRC
        REAL*8,         INTENT(IN)      :: DT
        REAL*8,         INTENT(IN)      :: OMEGA(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(IN)      :: NEGOMEGA(IIPAR,JJPAR,LLPAR)
        REAL*8,         INTENT(INOUT)   :: STT(:,:,:,:)
        REAL*8,         INTENT(IN)      :: TCVV(:)
        REAL*8,         INTENT(IN)      :: AD(:,:,:)
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!
        INTEGER         :: I, J, K, N, L
        INTEGER         :: internal_timestep
        REAL*8          :: ZMASS(IIPAR, JJPAR, LLPAR)
        REAL*8          :: Q_NEW(IIPAR, JJPAR, LLPAR)
        REAL*8          :: AREA_M2(IIPAR,JJPAR)
        REAL*8          :: TRACERFLUX(IIPAR, JJPAR, LLPAR)

! START EXECUTION
        ! Print info
        WRITE( 6, '(a)' ) 'Performing vertical flux exchange'

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

        ! Loop through each tracer
        DO N = 1, N_TRC

            ! Additional mass flux of each tracer is the mass fluxo f
            ! air multiplied by the mixign ratio of the tracer divided
            ! by the ratio of molec w eight of tracer to the molec
            ! weight of air
            ! [kg tracer] = [kg air] * [mol tracer/mol air] * 
            !               [kg tracer/mol tracer * mol air/kg air]
            TRACERFLUX(:,:,:) = ZMASS(:,:,:) * STT(:,:,:,N) &
                                / TCVV(N)

        print*, 't', tracerflux(10,10,1), tracerflux(10,10,2)
        print*, 'stt a', STT(10,10,1,N)*AD(10,10,1) / TCVV(N)
        print*, 'stt, ad, tcvv', STT(10,10,1,N), AD(10,10,1), TCVV(N)
            ! Bottom grid box first: new mass = mass currently in cell
            ! minus mass out + mass in from box above
            Q_NEW(:,:,1) = (STT(:,:,1,N) * AD(:,:,1)) / &
                           TCVV(N) - TRACERFLUX(:,:,1) + &
                           TRACERFLUX(:,:,2)

            DO L = 2, LLPAR-1
!$OMP PARALLEL DO           &
!$OMP DEFAULT( SHARED )     &
!$OMP PRIVATE( I, J )
            DO J = 1, JJPAR
            DO I = 1, IIPAR 
                ! New amass = mass currently in cell - 2 * mass out (one
                ! for up and one for down) + mass in from cell above +
                ! mass in from cell below
                Q_NEW(I,J,L) = (STT(I,J,L,N) * AD(I,J,L)) / &
                               TCVV(N) - 2*TRACERFLUX(I,J,L) + &
                            TRACERFLUX(I,J,L-1) + TRACERFLUX(I,J,L+1)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
            ENDDO

            ! Top layer
            ! New mass = mass currently in cell - mass out (down) + mass
            ! in from cell below
            Q_NEW(:,:,47) = STT(:,:,47,N) * AD(:,:,L) / &
                            TCVV(N) - TRACERFLUX(:,:,47) + &
                            TRACERFLUX(:,:,46)

            ! Reset original array with the new values, converting from
            ! [kg] to [v/v]
            STT(:,:,:,N) = Q_NEW(:,:,:) * TCVV(N)/AD(:,:,:)
            print*, 'stt', sum(STT(:,:,:,N))
            
        ENDDO

       END SUBROUTINE DO_FLUX_EXCHANGE
!EOC
      END MODULE FLUX_EXCHANGE_MOD
