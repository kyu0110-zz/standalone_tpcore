!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: convection_mod
!
! !DESCRIPTION: Module CONVECTION\_MOD contains routines which select the 
!  proper convection code for GEOS-3, GEOS-4, GEOS-5, MERRA, or GCAP met 
!  field data sets. 
!\\
!\\
! !INTERFACE: 
!
      MODULE CONVECTION_MOD
! 
! !USES:
!
      USE GIGC_Errcode_Mod
      USE CMN_GCTM_Mod

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_CONVECTION
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DO_MERRA_CONVECTION
!
! !REVISION HISTORY:
!  27 Jan 2004 - R. Yantosca - Initial version
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also now references
!        "error_mod.f".  Now make F in routine NFCLDMX a 4-D array to avoid
!        memory problems on the Altix. (bmy, 1/27/04)
!  (2 ) Bug fix: Now pass NTRACE elements of TCVV to FVDAS_CONVECT in routine 
!        DO_CONVECTION (bmy, 2/23/04)  
!  (3 ) Now references "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now also references "ocean_mercury_mod.f" and "tracerid_mod.f" 
!        (sas, bmy, 1/19/05)
!  (5 ) Now added routines DO_GEOS4_CONVECT and DO_GCAP_CONVECT by breaking 
!        off code from DO_CONVECTION, in order to implement GCAP convection
!        in a much cleaner way. (swu, bmy, 5/25/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Shut off scavenging in shallow convection for GCAP (swu, bmy, 11/1/05)
!  (8 ) Modified for tagged Hg simulation (cdh, bmy, 1/6/06)
!  (9 ) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (10) Fix for GEOS-5 met fields in routine NFCLDMX (swu, 8/15/07)
!  (11) Resize DTCSUM array in NFCLDMX to save memory (bmy, 1/31/08)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Added modifications for MERRA
!  05 Oct 2010 - R. Yantosca - Added ND14 and ND38 diagnostics to 
!                              DO_MERRA_CONVECTION routine
!  16 Aug 2011 - J. Fisher   - Minor bug fixes in DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  02 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  23 Jun 2015 - E. Lundgren - Convert tracer units from v/v dry air to
!                              kg/kg total air for convection
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
! !IROUTINE: do_convection
!
! !DESCRIPTION: Subroutine DO\_CONVECTION calls the appropriate convection 
!  driver program for different met field data sets.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_CONVECTION( DT, CMFMC, tracer_conc
     &                          )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE TIME_MOD,           ONLY : GET_TS_DYN

!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)    :: CMFMC(IIPAR, JJPAR, LLPAR+1) 
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT) :: tracer_conc   
!
! !OUTPUT PARAMETERS:
!
! 
! !REVISION HISTORY: 
!  08 Feb 2007 - R. Yantosca - Now reference "CMN_SIZE".  Now references 
!                              CLDMAS, CMFMC, DTRAIN from "dao_mod.f" so that 
!                              we can pass either GEOS-5 or GEOS-3 meteorology
!                              to NFCLDMX. 
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Now call DO_MERRA_CONVECTION for MERRA met
!  05 Oct 2010 - R. Yantosca - Now attach diagnostics to MERRA conv routine
!  06 Oct 2010 - R. Yantosca - Parallelized call to DO_MERRA_CONVECTION
!  15 Oct 2010 - H. Amos     - Now get BXHEIGHT, T from dao_mod.f
!  15 Oct 2010 - R. Yantosca - Now get LDYNOCEAN from logical_mod.f
!  15 Oct 2010 - R. Yantosca - Now get ITS_A_MERCURY_SIM from tracer_mod.f
!  15 Oct 2010 - R. Yantosca - Now get IDTHg2, IDTHgP from tracerid_mod.f
!  15 Oct 2010 - R. Yantosca - Now get H2O2s, SO2s from wetscav_mod.f
!  15 Oct 2010 - H. Amos     - Now pass BXHEIGHT, T, to DO_MERRA_CONVECTION
!  15 Oct 2010 - R. Yantosca - Now pass H2O2s, SO2s to DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  29 Aug 2011 - R. Yantosca - Bug fix: reposition #if defined(APM) statement
!  09 Feb 2012 - R. Yantosca - For GEOS-5.7, PFICU and PFLCU fields are
!                              defined on level edges.  Pass the top edge
!                              of each level to DO_MERRA_CONVECTION
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  21 Jun 2012 - R. Yantosca - Now use poiners to pass array slices to routines
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  17 Apr 2014 - R. Yantosca - Speed up MERRA and GEOS-FP convection by
!                              adding !$OMP+SCHEDULE( DYNAMIC )
!  17 Apr 2014 - R. Yantosca - Test if (ND14>0) and (ND38>0) outside OMP loop
!  18 Apr 2014 - R. Yantosca - Now use a pointer to pass a slice of the
!                              F array to subroutine COMPUTE_F
!  18 Apr 2014 - R. Yantosca - Clean up call to DO_MERRA_CONVECTION, remove
!                              stuff that was leftover from the column code
!  18 Apr 2014 - R. Yantosca - Now use proper # of tracers for APM in the
!                              call to DO_MERRA_CONVECTION (MERRA & GEOS-FP)
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to COMPUTE_F
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  11 Aug 2015 - R. Yantosca - Added support for MERRA2 data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars

      ! Pointers
      ! Define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. Units are 
      ! initially v/v dry air but will be converted to
      ! kg/kg total air for convection (ewl, 23/15)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! MERRA or GEOS-FP met fields definitions
      !=================================================================
      
      ! Scalars          
      INTEGER           :: I, J, L, N, NN, TS_DYN
      REAL*8          :: AREA_M2, DT
      
      ! Arrays        
      REAL*8          :: PEDGE  (LLPAR+1)

      !=================================================================
      ! DO_CONVECT begins here!
      !=================================================================

      ! Convert State_Chm%Tracers from [v/v dry air] to 
      ! [kg/kg total air] prior to convection (ewl, 6/23/15)
      !CALL Convert_DryVV_to_DryKgKg( am_I_Root, N_TRACERS, 
      !CALL Convert_DryKgKg_to_MoistKgKg( am_I_Root, N_TRACERS,  

      !=================================================================
      ! MERRA or GEOS-FP met fields
      !=================================================================

      ! Initialize
      TS_DYN    = GET_TS_DYN()                    ! Dynamic timestep [min]
      DT        = DBLE( TS_DYN )                  ! Dynamic timestep [min]

      !-----------------------------------------------------------------
      ! Do convection column by column
      !-----------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,  I, AREA_M2, L, PEDGE, N, NN )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         AREA_M2     =  GET_AREA_M2( I, J, 1 ) 

         ! Pressure edges
         DO L = 1, LLPAR+1
            PEDGE(L) =  GET_PEDGE( I, J, L, Ap, Bp, PS )
         ENDDO
         
         !--------------------------
         ! Do the cloud convection
         !--------------------------
         CALL DO_MERRA_CONVECTION( CMFMC      = CMFMC(I, J, 2:LLPAR+1),
     &                             Q          = STT(I, J, :, :),
     &                             I          = I,
     &                             J          = J,
     &                             N_TRACERS  = 1,
     &                             AREA_M2    = AREA_M2,
     &                             PEDGE      = PEDGE,
     &                             TS_DYN     = DT          )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Convert tracer concentration from kg/kg total air to 
      ! mol/mol dry air after convection (ewl, 6/23/15)
      !CALL Convert_MoistKgKg_to_DryKgKg( am_I_Root, N_TRACERS,  
      !CALL Convert_DryKgKg_to_DryVV( am_I_Root, N_TRACERS, 

      END SUBROUTINE DO_CONVECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_merra_convection
!
! !DESCRIPTION: Subroutine DO\_MERRA\_CONVECTION (formerly called NFCLDMX)
!  is S-J Lin's cumulus transport module for 3D GSFC-CTM, modified for the 
!  GEOS-Chem model.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_MERRA_CONVECTION( CMFMC,  
     &                                Q,
     &                                I,          
     &                                J,
     &                                N_TRACERS,  
     &                                AREA_M2,    
     &                                PEDGE,      
     &                                TS_DYN,
     &                                BXHEIGHT,
     &                                DQRCU,
     &                                DTRAIN,
     &                                DELP,
     &                                T )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP  ! hma Nov 3, debug
!
! !INPUT PARAMETERS:
!  
      REAL*8,       INTENT(IN)    :: CMFMC(:)
      INTEGER,        INTENT(IN)    :: I, J       ! Lon & lat indices
      INTEGER,        INTENT(IN)    :: N_TRACERS  ! # of convected tracers
      REAL*8,       INTENT(IN)    :: AREA_M2    ! Surface area [m2]
      REAL*8,       INTENT(IN)    :: PEDGE(:)   ! P @ level box edges [hPa]
      REAL*8,       INTENT(IN)    :: TS_DYN     ! Dynamic timestep [min]
      REAL*8,  INTENT(IN) :: BXHEIGHT (:        )
      REAL*8,  INTENT(IN)  :: DQRCU    (:        )
      REAL*8,  INTENT(IN)  :: DTRAIN   (:        )
      REAL*8,  INTENT(IN)  :: DELP     (:        )
      REAL*8,  INTENT(IN)  :: T        (:        )
!                                                  
! !INPUT/OUTPUT PARAMETERS:                        
!        
      REAL*8,       INTENT(INOUT) :: Q(:)  ! Tracer array
!                                                  
! !OUTPUT PARAMETERS:                              
!                     
!
! !REMARKS:
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!                                                                             .
!  Unit conversion for BMASS:
!
!      Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!     -------------+---------+--------+-------+----+--------  =  -----
!                  | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!                                                                             .
!  NOTE: We are passing I & J down to this routine so that it can call the
!  proper code from "mercury_mod.f".  Normally, we wouldn't pass I & J as
!  arguments to columnized code.  This prevents rewriting the mercury_mod.f
!  routines ADD_Hg2_
!
! !REVISION HISTORY:
!  15 Jul 2009 - R. Yantosca - Columnized and cleaned up.
!                            - CLDMAS renamed to CMFMC and DTRN renamed
!                              to DTRAIN for consistency w/ GEOS-5.
!  17 Jul 2009 - R. Yantosca - Now do unit conversion of Q array from
!                              [kg] --> [v/v] and vice versa internally
!  14 Dec 2009 - R. Yantosca - Now remove internal unit conversion, since
!                              Q now comes in as [mol/mol] (=[v/v]) from the
!                              calling routine.
!  14 Dec 2009 - R. Yantosca - Remove COEF from the argument list
!  06 May 2010 - R. Yantosca - Now add IDENT via the argument list
!  29 Sep 2010 - R. Yantosca - Modified for MERRA met fields
!  05 Oct 2010 - R. Yantosca - Now pass COEF via the argument list
!  05 Oct 2010 - R. Yantosca - Attach ND14 and ND38 diagnostics
!  15 Oct 2010 - H. Amos     - Added BXHEIGHT and T as arguments
!  15 Oct 2010 - R. Yantosca - Added I, J, H2O2s and SO2s as arguments
!  15 Oct 2010 - H. Amos     - Added scavenging below cloud base
!  06 Apr 2011 - M.Fu, H.Amos- Bug fix: make sure washout adheres to the same
!                              algorithm as in the wet deposition code.
!  27 Jul 2011 - R. Yantosca - Declare CLDBASE as INTEGER to avoid PGI errors
!  16 Aug 2011 - J. Fisher   - Bug fix: use IS_Hg2() and IS_HgP to test if 
!                              a tracer is Hg2 or HgP (for tagged species)
!  16 Aug 2011 - J. Fisher   - Now use WETLOSS instead of T0_SUM in the ND38 
!                              diagnostic below the cloud.  Using T0_SUM leads 
!                              us to over-count the tracer scavenged out of 
!                              the column.
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  31 May 2013 - R. Yantosca - Now pass State_Chm to WASHOUT
!  05 Sep 2013 - R. Yantosca - Bug fix: DT is apparently undefined, but still
!                              passed to WASHOUT.  Use SDT instead.  This
!                              avoids a floating-point error.
!  18 Apr 2014 - R. Yantosca - Now point to 3-D arrays internally
!  18 Apr 2014 - R. Yantosca - Now also pass N_TRACERS (to facilitate APM)
!  18 Apr 2014 - R. Yantosca - Remove code that we don't need anymore
!  04 Feb 2015 - M. Sulprizio- Fix calculation of WETLOSS for non-aerosol 
!                              tracers below the cloud base (C. Friedman)
!  20 Apr 2015 - E. Lundgren - Use DELP*100/g instead of AD/area as grid box
!                              moist mass per area and remove AD from routine
!  20 May 2015 - M. Sulprizio- Apply bug fixes provided by Viral Shah:
!                              -- Remove F(K,IC) > 0 condition that prevents
!                                 archiving of deposited mass in DIAG38
!                              -- Add statement that subtracts the wet deposited
!                                 amount from the atmospheric mass
!                              -- Fix inconsistency in units when T0_SUM is used
!  04 Jun 2015 - E. Lundgren - Adapt Viral Shah bug fixes to moist mixing ratio
!  09 Jun 2015 - R. Yantosca - Now deposit Hg2, HgP to snowpack regardless of
!                              whether the dynamic ocean is used
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  22 Jun 2015 - E. Lundgren - Move QB_NUM calculation to within timestep loop
!  12 Aug 2015 - R. Yantosca - Treat MERRA2 in same way as we do for GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: TINYNUM = 1d-14
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER             :: IC,          ISTEP,    K,       KTOP
      INTEGER             :: NC,          NDT,      NLAY,    NS    
      INTEGER             :: CLDBASE
      REAL*8            :: CMFMC_BELOW, ALPHA,    ALPHA2
      REAL*8            :: CMOUT,       DELQ,     DQ,      DNS
      REAL*8            :: ENTRN,       QC,       QC_PRES, QC_SCAV 
      REAL*8            :: SDT,         T0,       T0_SUM,  T1    
      REAL*8            :: T2,          T3,       T4 
      REAL*8            :: TSUM
      REAL*8            :: LOST,        GAINED
      REAL*8            :: MASS_NOWASH, QDOWN,    DT
      REAL*8            :: MB,          QB,       QB_NUM

      ! Arrays
      REAL*8            :: BMASS       (LLPAR    )    
      REAL*8            :: PDOWN       (LLPAR    )

      !========================================================================
      ! (0)  I n i t i a l i z a t i o n
      !========================================================================
      
      ! # of levels and # of tracers
      NLAY     = LLPAR

      ! Top level for convection
      KTOP     = NLAY - 1

      ! Convection timestep [s]
      NDT      = TS_DYN * 60d0

      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      NS       = NDT / 300                ! Num internal timesteps (int)
      NS       = MAX( NS, 1 )             ! Set lower bound to 1
      DNS      = DBLE( NS )               ! Num internal timesteps (real)
      SDT      = DBLE( NDT ) / DBLE( NS ) ! seconds in internal timestep
      
      !-----------------------------------------------------------------
      ! Determine location of the cloud base, which is the level where
      ! we start to have non-zero convective precipitation formation
      !-----------------------------------------------------------------

      ! Minimum value of cloud base is the surface level
      CLDBASE = 1

      ! Find the cloud base
      DO K = 1, NLAY
         IF ( DQRCU(K) > 0d0 ) THEN
            CLDBASE = K
            EXIT
         ENDIF
      ENDDO

      !-----------------------------------------------------------------
      ! Compute PDOWN and BMASS
      !-----------------------------------------------------------------
      DO K = 1, NLAY

         ! BMASS is the total air mass per unit area for the grid box
         ! bounded by level K and K+1 [kg/m2] 
         ! BMASS is equivalent to deltaP * 100 / g
         ! This is done to keep BMASS in the same units as CMFMC * SDT
         BMASS(K) = DELP(K) * G0_100

      ENDDO

      !-----------------------------------------------------------------
      ! Compute MB, the mass per unit area of total air below the cloud
      ! base [kg/m2]. Calculate MB by looping over levels below the 
      ! cloud base.
      !-----------------------------------------------------------------
      MB = 0d0
      DO K = 1, CLDBASE-1
         MB = MB + BMASS(K)
      ENDDO


      !========================================================================
      ! (1)  T r a c e r   L o o p 
      !========================================================================
      DO IC = 1, NC

         !=====================================================================
         ! (2)  I n t e r n a l   T i m e   S t e p   L o o p
         !=====================================================================
         DO ISTEP = 1, NS

            ! Initialize
            QC     = 0d0    ! [kg tracer/kg total air]                  
            T0_SUM = 0d0    ! [kg tracer/m2/timestep]

            !----------------------------------------------------------
            ! B e l o w   C l o u d   B a s e   (K < CLDBASE) 
            !
            ! QB is the "weighted avg" mixing ratio below the cloud
            ! base [kg/kg total air].
            ! QC is the mixing ratio of the air that moved in cumulus 
            ! transport up to the next level [kg/kg total air]. 
            ! MB is the total mass of air below the cloud base per 
            ! unit area [kg/m2] (see calculation before loop).
            !-----------------------------------------------------------

            ! We need to make this a nested IF statement so that we don't
            ! get an out-of-bounds error when CLDBASE=1 (bmy, 11/18/10)
            IF ( CLDBASE > 1 ) THEN
               
               IF ( CMFMC(CLDBASE-1) > TINYNUM ) THEN
 
                  !-----------------------------------------------------
                  ! %%% Non-negligible Cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! Calculate QB_NUM, the numerator for QB. QB is the 
                  ! weighted average mixing ratio below the cloud base.
                  ! QB_NUM is equal to the grid box tracer concentrations 
                  ! [kg/kg total air] weighted by the adjacent level pressure 
                  ! differences and summed over all levels up to just 
                  ! below the cloud base (ewl, 6/22/15)           
                  QB_NUM  = 0d0  
                  DO K  = 1, CLDBASE-1
                     QB_NUM = QB_NUM + Q(K,IC) * DELP(K)
                  ENDDO

                  ! Compute QB, the weighted avg mixing ratio below
                  ! the cloud base [kg/kg total air]
                  QB = QB_NUM / ( PEDGE(1) - PEDGE(CLDBASE) )

                  ! Compute QC, the mixing ratio of the air that moved
                  ! in cumulus transport up to the next level [kg/kg]
                  !
                  !        Total mass of tracer below cloud base  +   
                  !        Subsidence into cloud base from above  
                  ! QC =  --------------------------------------------
                  !            Total air mass below cloud base 
                  !
                  QC = ( MB*QB + CMFMC(CLDBASE-1) * 
     &                           Q(CLDBASE,IC)    * SDT  ) /
     &                 ( MB    + CMFMC(CLDBASE-1) * SDT  ) 

                  ! Copy QC to all levels of the tracer array Q
                  ! that are below the cloud base level [kg/kg]
                  Q(1:CLDBASE-1,IC) = QC

               ELSE                  

                  !-----------------------------------------------------
                  ! %%% Negligible cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! When CMFMC is negligible, then set QC to the tracer
                  ! concentration at the cloud base level [kg/kg]
                  QC = Q(CLDBASE,IC)
   
               ENDIF

            ELSE
               
               !-----------------------------------------------------
               ! If the cloud base happens at level 1, then just
               ! set QC to the tracer concentration at the surface 
               ! level [kg/kg]
               !-----------------------------------------------------
               QC = Q(CLDBASE,IC) 

            ENDIF

            !==================================================================
            ! (3)  A b o v e   C l o u d   B a s e
            !==================================================================
            DO K = CLDBASE, KTOP
               
               ! Initialize
               ALPHA   = 0d0
               ALPHA2  = 0d0
               CMOUT   = 0d0
               ENTRN   = 0d0
               QC_PRES = 0d0
               
               ! CMFMC_BELOW is the air mass [kg/m2/s] coming into the
               ! grid box (K) from the box immediately below (K-1).
               IF ( K == 1 ) THEN
                  CMFMC_BELOW = 0d0
               ELSE
                  CMFMC_BELOW = CMFMC(K-1)
               ENDIF

               ! If we have a nonzero air mass flux coming from 
               ! grid box (K-1) into (K) ...
               IF ( CMFMC_BELOW > TINYNUM ) THEN

                  !------------------------------------------------------------
                  ! (3.1)  M a s s   B a l a n c e   i n   C l o u d
                  !
                  ! F(K,IC) = fraction of tracer IC in level K that is 
                  !           available for wet-scavenging by cloud updrafts.  
                  !
                  ! If ENTRN > 0 then compute the new value of QC:
                  !
                  !      tracer mass from below      (i.e. level K-1) + 
                  !      tracer mass from this level (i.e. level K)
                  !  = -----------------------------------------------------
                  !             total mass coming into cloud
                  !
                  ! Otherwise, preserve the previous value of QC.  This will 
                  ! ensure that TERM1 - TERM2 is not a negative quantity (see 
                  ! below).
                  !  
                  ! Entrainment must be >= 0 (since we cannot have a negative 
                  ! flux of air into the cloud).  This condition is strong 
                  ! enough to ensure that CMOUT > 0 and will prevent floating-
                  ! point exception.
                  !------------------------------------------------------------

                  ! Air mass flowing out of cloud at grid box (K) [kg/m2/s]
                  CMOUT   = CMFMC(K) + DTRAIN(K)

                  ! Air mass flowing into cloud at grid box (K) [kg/m2/s]
                  ENTRN   = CMOUT - CMFMC_BELOW

                  ! Amount of QC preserved against scavenging [kg/kg]
                  QC_PRES = QC * ( 1d0 - F(K,IC) )

                  ! Update QC taking entrainment into account [kg/kg]
                  ! Prevent div by zero condition
                  IF ( ENTRN >= 0e+0_fp .and. CMOUT > 0e+0_fp ) THEN
                     QC   = ( CMFMC_BELOW * QC_PRES   + 
     &                       ENTRN       * Q(K,IC) ) / CMOUT
                  ENDIF

                  !------------------------------------------------------------
                  ! (3.2)  M a s s   B a l a n c e   i n   L e v e l  ==> Q
                  !
                  ! Terminology:
                  !
                  !  C_k-1   = cloud air mass flux from level k-1 to level k
                  !  C_k     = cloud air mass flux from level k   to level k+1
                  !  QC_k-1  = mixing ratio of tracer INSIDE CLOUD at level k-1
                  !  QC_k    = mixing ratio of tracer INSIDE CLOUD at level k
                  !  Q_k     = mixing ratio of tracer in level k
                  !  Q_k+1   = mixing ratio of tracer in level k+1
                  ! 
                  ! For convenience we denote:
                  !
                  !  QC_SCAV = Amount of tracer wet-scavenged in updrafts
                  !          = QC_k-1 * F(k,IC)    [kg/kg]
                  !
                  !  QC_PRES = Amount of tracer preserved against
                  !            wet-scavenging in updrafts [kg/kg]
                  !          = QC_k-1 * ( 1 - F(k,IC) )   
                  !
                  ! Where F(k,IC) is the fraction of tracer IC in level k
                  ! that is available for wet-scavenging by cloud updrafts.
                  ! F(k,IC) is computed by routine COMPUTE_UPDRAFT_FSOL
                  ! and passed to this routine as an argument.
                  !
                  ! The cumulus transport above the cloud base is done as 
                  ! follows:
                  !  
                  !                 ||///////////////////||
                  !                 ||//// C L O U D ////||
                  !                 ||                   ||
                  !   k+1     ^     ||         ^         ||3)   C_k * Q_k+1
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k      C_k    ||2)   C_k * QC_k    ||         V
                  !                 ||                   ||
                  !                 ||                   ||
                  !           ^     ||         ^         ||4)   C_k-1 * Q_k 
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k-1   C_k-1   ||1) C_k-1 * QC_k-1  ||         V
                  !                 ||         * (1 - F) ||
                  !                 ||                   ||
                  !                 ||//// C L O U D ////||
                  !                 ||///////////////////||
                  !
                  ! There are 4 terms that contribute to mass flow in 
                  ! and out of level k:
                  !
                  ! 1) C_k-1 * QC_PRES = tracer convected from k-1 to k 
                  ! 2) C_k   * QC_k    = tracer convected from k   to k+1 
                  ! 3) C_k   * Q_k+1   = tracer subsiding from k+1 to k 
                  ! 4) C_k-1 * Q_k     = tracer subsiding from k   to k-1 
                  !
                  ! Therefore the change in tracer concentration is given by
                  !
                  !    DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
                  !
                  ! and Q(K,IC) = Q(K,IC) + DELQ.  
                  !
                  ! The term T0 is the amount of tracer that is scavenged 
                  ! out of the box.
                  !
                  ! Units of T0, T1, T2, T3, T4, and TSUM are
                  ! [kg/m2/s * kg tracer / kg total air]
                  !------------------------------------------------------------
                  T0      =  CMFMC_BELOW * QC_SCAV      
                  T1      =  CMFMC_BELOW * QC_PRES       
                  T2      = -CMFMC(K  )  * QC
                  T3      =  CMFMC(K  )  * Q(K+1,IC)
                  T4      = -CMFMC_BELOW * Q(K,  IC)
                 
                  TSUM    = T1 + T2 + T3 + T4          
                  
                  DELQ    = ( SDT / BMASS(K) ) * TSUM    ! change in [kg/kg]

                  ! If DELQ > Q then do not make Q negative!!!
                  IF ( Q(K,IC) + DELQ < 0 ) THEN
                     DELQ = -Q(K,IC)
                  ENDIF

                  ! Increment the tracer array [kg/kg]
                  Q(K,IC) = Q(K,IC) + DELQ

                  ! check for infinity
                  IF ( .not. IT_IS_FINITE( Q(K,IC) ) ) THEN
                     CALL GEOS_CHEM_STOP
                  ENDIF

                  ! Return if we encounter NaN
                  IF ( IT_IS_NAN( Q(K,IC) ) ) THEN 
                     WRITE( 6, 250 )
                     WRITE( 6, 255 ) K, IC, Q(K,IC) 
 250                 FORMAT( 'NaN encountered in DO_MERRA_CONVECTION!' )
 255                 FORMAT( 'K, IC, Q(K,IC): ', 2i4, 1x, es13.6 )
                     RC = GIGC_FAILURE
                     RETURN
                  ENDIF

                  ! Pass T0_SUM in units of [kg tracer/m2/timestep].
                  ! Converting kg total air to kg tracer requires use
                  ! of the molecular weight of air including moisture
                  ! (ewl, 6/5/15)
                  T0_SUM = T0_SUM + T0 * SDT

               ELSE

                  !------------------------------------------------------------
                  ! (3.5)  N o   C l o u d   M a s s   F l u x   B e l o w 
                  !------------------------------------------------------------

                  ! If there is no cloud mass flux coming from below, set
                  ! QC to the tracer concentration at this level [kg/kg]
                  QC = Q(K,IC)
                  
                  ! Bug fix for the cloud base layer, which is not necessarily
                  ! in the boundary layer, and for MERRA, there could be 
                  ! "secondary convection" plumes - one in the PBL and another 
                  ! one not.  NOTE: T2 and T3 are the same terms as described 
                  ! in the above section.  (swu, 08/13/2007)
                  IF ( CMFMC(K) > TINYNUM ) THEN 

                     ! Tracer convected from K -> K+1 
                     ! [kg/m2/s * kg tracer/kg total air]
                     T2   = -CMFMC(K) * QC

                     ! Tracer subsiding from K+1 -> K [kg/m2/s]
                     ! [kg/m2/s * kg tracer/kg total air]
                     T3   =  CMFMC(K) * Q(K+1,IC)

                     ! Change in tracer concentration [kg/kg]
                     DELQ = ( SDT / BMASS(K) ) * (T2 + T3)

                     ! If DELQ > Q then do not make Q negative!!!
                     IF ( Q(K,IC) + DELQ < 0.0e+0_fp ) THEN 
                        DELQ = -Q(K,IC)
                     ENDIF
  
                     ! Add change in tracer to Q array [kg/kg]
                     Q(K,IC) = Q(K,IC) + DELQ

                  ENDIF
               ENDIF
            ENDDO     ! End of loop over levels above cloud base
            
         ENDDO      ! End internal timestep loop
      ENDDO      ! End loop over tracers

      !================================================================
      ! Succesful return!
      !================================================================
      
      ! Nullify pointers
      NULLIFY( BXHEIGHT )
      NULLIFY( CMFMC    )
      NULLIFY( DQRCU    )
      NULLIFY( DTRAIN   )
      NULLIFY( PFICU    )
      NULLIFY( PFLCU    )
      NULLIFY( REEVAPCN )
      NULLIFY( DELP     )
      NULLIFY( T        )
      NULLIFY( Q        )

      ! Set error code to success
      RC                      = GIGC_SUCCESS

      END SUBROUTINE DO_MERRA_CONVECTION
!EOC
      END MODULE CONVECTION_MOD
