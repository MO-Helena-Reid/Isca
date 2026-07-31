module bechtold_alg_mod
  use mo_cumastrn, only: cumastrn
  
  USE mo_kind   ,ONLY: JPRB=>wp, vp   , &
    &                  jpim=>i4
  USE constants_mod, ONLY: rdgas
  implicit none
  subroutine bechtold_alg( is,     js,      Time,      temp0,   qvap0,     &
                  uwnd0,  vwnd0, omega0, pres0,   pres0_int,      &
                  flux_t0, flux_q0, zfull0, zhalf0,  coldT0,      &
                  dtime,  dtemp0,  dqvap0,    duwnd0,  dvwnd0,    &
                  rain0,  snow0,   do_strat,                      &
                  klzbs,  klcls,                                  &
                  !OPTIONAL IN
                  mask,    kbot,                                  &
                  !OPTIONAL OUT
                  mc0,    ql0, qi0, qa0,  dl0, di0, da0,          &
                  ras_tracers, qtrras)
  !---------------------------------------------------------------------
! Arguments (Intent in)
!     is, js    - Starting indices for window
!     Time      - Time used for diagnostics [time_type]
!     dtime     - Size of time step in seconds
!     omega0    - Vertical velocity in Pa/s
!     pres0     - Pressure
!     pres0_int - Pressure at layer interface
!     flux_t0   - Sensible heat flux
!     flux_q0   - Moisture flux
!     zfull0    - Height at layer centre
!     zhalf0    - Height at layer interface
!     temp0     - Temperature
!     qvap0     - Water vapor 
!     uwnd0     - U component of wind
!     vwnd0     - V component of wind
!     coldT0    - should the precipitation assume to be frozen?
!Optional:
!     R0        - OPTIONAL;prognostic tracers to move around
!                 note that R0 is assumed to be dimensioned
!                 (nx,ny,nz,nt), where nt is the number of tracers
!     ql0       - OPTIONAL;cloud liquid
!     qi0       - OPTIONAL;cloud ice
!     qa0       - OPTIONAL;cloud/saturated volume fraction
!     mask      - OPTIONAL; mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!     kbot      - OPTIONAL;index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!---------------------------------------------------------------------
  type(time_type), intent(in)                   :: Time
  integer,         intent(in)                   :: is, js
  real,            intent(in)                   :: dtime
  real,            intent(in), dimension(:,:,:) :: pres0, pres0_int, zhalf0, omega0
  real,            intent(in), dimension(:,:,:) :: temp0, qvap0, uwnd0, vwnd0
  logical,         intent(in), dimension(:,:)   :: coldT0
  logical,         intent(in)                   :: do_strat
  real,  intent(inout), dimension(:,:,:) ,optional :: ql0, qi0, qa0
  real,  intent(in), dimension(:,:,:,:),optional :: ras_tracers
  !optional
  real, intent(in) , dimension(:,:,:), OPTIONAL :: mask
  integer, intent(in), OPTIONAL, dimension(:,:) :: kbot
!---------------------------------------------------------------------
! Arguments (Intent out)
!       rain0  - surface rain
!       snow0  - surface snow
!       dtemp0 - Temperature change 
!       dqvap0 - Water vapor change 
!       duwnd0 - U wind      change 
!       dvwnd0 - V wind      change 
!       mc0    - OPTIONAL; cumulus mass flux
!       Dl0    - OPTIONAL; cloud liquid change
!       Di0    - OPTIONAL; cloud ice change
!       Da0    - OPTIONAL; cloud fraction change
!       DR0    - OPTIONAL; increment to prognostic tracers
!---------------------------------------------------------------------
  integer, intent(out), dimension(:,:)  :: klzbs, klcls
  real, intent(out), dimension(:,:,:) :: dtemp0, dqvap0, duwnd0, dvwnd0
  real, intent(out), dimension(:,:)   :: rain0,  snow0
  
  real, intent(out), OPTIONAL, dimension(:,:,:) :: mc0
  real, intent(out), OPTIONAL, dimension(:,:,:) :: dl0, di0, da0
  real,  intent(out), dimension(:,:,:,:), optional :: qtrras
!---------------------------------------------------------------------
integer :: klon, imax, jmax, kmax, klev, j,k

real, dimension(SIZE(temp0,1),SIZE(temp0,3)) :: plitot, zdgeoh, zdph,&
                                                ptenta,ptenqa,ptenrhoq,&
                                                ptenrhol,ptenrhoi,&
                                                ptenrhor,ptenrhos,&
real, dimension(SIZE(temp0,1),SIZE(temp0,2),SIZE(temp0,3)) :: rho_full

  ! --- Set dimensions
  imax  = size( temp0, 1 )
  jmax  = size( temp0, 2 )
  kmax  = size( temp0, 3 )
  klon = jmax
  klev = kmax
  ! Calculate air density
  rho_full(:,:,:) = pres0(:,:,:) / (rdgas * temp0(:,:,:))
  DO j=1,jmax
    ! Set these output variables to zero just in case
    ptenrhol = 0.0
    ptenrhoi = 0.0
    ptenrhor = 0.0
    ptenrhos = 0.0
    ! If the dynamic correction of cape closure tuning factor 
    ! is positive these will be used.
    ! If llo1 is ever true ptenqa will be used.
    ! These are supposed to be the dynamics moisture tendency and
    ! temperature tendency respectively but these look like a pain to
    ! get in Isca so just setting them to zero.
    ptenqa = 0.0
    ptenta = 0.0


    IF (PRESENT(ql0) .and. PRESENT(qi0)) THEN
      plitot = ql0(:,j,:) + qi0(:,j,:)
    ELSE
      plitot = 0.0
    ENDIF
    DO k=1,kmax
      zdgeoh(:,k) = ABS(zhalf0(:,k)-zhalf0(:,k+1))
      zdph(:,k) = ABS(pres0_int(:,k)-pres0_int(:,k+1))
    END DO
    CALL cumastrn &
  & (  kidia,    kfdia,    klon,   ktdia,   klev, &
  & ldland, ldlake, dtime, phy_params, k950,      &
  & trop_mask, mtnmask,  paer_ss,                 &
  & temp0(:,j,:),qvap0(:,j,:),                    &
  & uwnd0(:,j,:), vwnd0(:,j,:), plitot,           &
  & omega0(:,j,:),                                &
  & ql0(:,j,:), qi0(:,j,:), flux_t0(:,j),         &
  & flux_q0(:,j), flux_q0(:,j),    flux_t0(:,j),  &
  & pres0(:,j,:), pres0_int(:,j,:), zfull0(:,j,:),&
  & zhalf0(:,j,:),                                &
  & zdph,               zdgeoh,                   &
  & dtent0(:,j,:),duwnd0(:,j,:),dvwnd0(:,j,:),    &
  & ptenta, ptenqa,                               &
  & dqvap0(:,j,:), ptenrhoq, ptenrhol, ptenrhoi,  &
  & ptenrhor, ptenrhos,                           &
  & ldcum,      ktype , kcbot,    kctop,          &
  & LDSHCV,   fac_entrorg, fac_rmfdeps,           &
  & pmfu,     pmfd,                               &
  & pmfude_rate,        pmfdde_rate,              &
  & ptu,      pqu,      plu,  pcore,              &
  & pmflxr,   pmflxs,   prain, pdtke_con,         &
  & pcape,    pvddraf,                            &
  & pcen, ptenrhoc,                               &
  & l_lpi, l_lfd, lpi, mlpi, koi, lfd, peis,      &
  & pertb,                                        &
  & lspinup, k650,k700, temp_s,                   &
  & cell_area,iseed,                              &
  & mf_bulk,mf_perturb,mf_num,p_cloud_ensemble,   &
  & pclnum_a, pclmf_a, pclnum_p, pclmf_p,         &
  & pclnum_d, pclmf_d, lacc                       )


! Code Description:
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KSTART*       FIRST STEP OF MODEL

!    *k650*         LEVEL INDEX AT 650hPa
!    *iseed*        SEED FOR RANDOM NUMBER GENERATOR

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDLAKE*       LAKE MASK (.TRUE. FOR LAKE)
!    *lspinup*      SPINUP CONVECTIVE CLOUD ENSEMBLE
!    *L_LPI*        COMPUTE LPI, MLPI, KOI
!    *L_LFD*        COMPUTE LFD

!     INPUT PARAMETERS (REAL)

!    paer_ss    monthly aerosol climatology sea salt (optical thickness)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *SHFL_S*       SENSIBLE HEAT FLUX (never halo-averaged)      W/M2
!    *QHFL_S*       MOISTURE FLUX (never halo-averaged)           KG/(SM2)
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *zdgeoh*       geopot thickness on full levels               M2/S2
!    *zdph*         pressure thickness on full levels               PA
!
!    *PSSTRU*       SURFACE MOMENTUM FLUX U                - not used presently
!    *PSSTRV*       SURFACE MOMENTUM FLUX V                - not used presently
!
!    *PTENTA*       TEMPERATURE TENDENCY DYNAMICS=TOT ADVECTION    K/S
!    *PTENQA*       MOISTURE    TENDENCY DYNAMICS=TOT ADVECTION    1/S

!    *temp_s*       TEMPERATURE IN LOWEST MODEL LEVEL                K
!    *cell_area*    GRID CELL AREA                                  M2?
!!!  FOR SPP
!    *pertb*        STOCHASTIC PATTERN FOR PERTURBATION
!!!  ALLOCATED ONLY IF lstoch_sde=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_a*     ACTIVE CLOUD NUMBER (T)               M-2
!    *pclmf_a*      ACTIVE MASS FLUX (T)                KG/(M2*S)
!    *pclnum_p*     PASSIVE CLOUD NUMBER (T)              M-2
!    *pclmf_p*      PASSIVE  MASS FLUX (T)              KG/(M2*S)

!   !!!  ALLOCATED ONLY IF lstoch_deep=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_d* PROGNOSTIC DEEP CLOUD NUMBER (T)               M-2
!    *pclmf_d*  PROGNOSTIC DEEP MASS FLUX (T)                KG/(M2*S)

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!    *PTENRHOC*     TENDENCY OF CHEMICAL TRACERS                  KG/(M3*S)

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS
!!!  ALLOCATED ONLY IF lstoch_sde=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_a*     PROGNOSTIC ACTIVE CLOUD NUMBER (T+1)               M-2
!    *pclmf_a*      PROGNOSTIC ACTIVE MASS FLUX (T+1)                KG/(M2 S)
!    *pclnum_p*     PROGNOSTIC PASSIVE CLOUD NUMBER (T+1)              M-2
!    *pclmf_p*      PROGNOSTIC PASSIVE  MASS FLUX (T+1)              KG/(M2 S)
!!!  ALLOCATED ONLY IF lstoch_deep=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_d*     PROGNOSTIC DEEP CLOUD NUMBER (T+1)                 M-2
!    *pclmf_d*      PROGNOSTIC DEEP MASS FLUX (T+1)                  KG/(M2 S)

!    OUTPUT PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PCORE*        UPDRAFT CORE FRACTION                         0-1
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M2*S)
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PDTKE_CON     CONV. BUOYANT TKE-PRODUCTION AT HALF LEVELS   M2/S**3)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S)
!    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
!    *PWMEAN*       VERTICALLY AVERAGED UPDRAUGHT VELOCITY         M/S
!    *pvddraf*      convective gust at surface                     M/S
!    *PTENRHOQ*     MOISTURE MASS DENSITY TENDENCY                KG/(M3*S)
!    *PTENRHOL*     LIQUID WATER MASS DENSITY TENDENCY            KG/(M3*S)
!    *PTENRHOI*     ICE CONDENSATE MASS DENSITY TENDENCY          KG/(M3*S)
!    *PTENRHOR*     DETRAINED RAIN MASS DENSITY TENDENCY          KG/(M3*S)
!    *PTENRHOS*     DETRAINED SNOW MASS DENSITY TENDENCY          KG/(M3*S)
!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG
!    *MLPI*         MODIFIED LPI USING KOI
!    *KOI*          KONVEKTIONS INDEX                             K
!    *LFD*          LIGHTNING FLASH DENSITY AS IN LOPEZ(2016)     1/(KM2*DAY)
!    *mf_bulk*      CLOUD BASE MASS FLUX FROM T-B SCHEME          KG/(M2*S)
!    *mf_perturb*   CLOUD BASE MASS FLUX FROM STOCHASTIC SCHEME   KG/(M2*S)
!    *mf_num*       NUMBER OF SHALLOW CLOUDS STOCHASTIC SCHEME    1


!     EXTERNALS.
!     ----------

!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V

!     SWITCHES.
!     --------

!          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFIT=.TRUE.    UPDRAUGHT ITERATION
!          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
!          LMFTRAC=.false. TRACER TRANSPORT

!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     RMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     RMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     RMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN

!     REFERENCE.
!     ----------

!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)
!          Bechtold et al. (2008 QJRMS 134,1337-1351), Rooy et al. (2012 QJRMS)
!          Bechtold et al. (2013 JAS)

!     AUTHOR.
!     -------
!      M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
END DO
  end subroutine bechtold_alg
end module bechtold_alg_mod

