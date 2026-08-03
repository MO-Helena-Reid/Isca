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
                  klzbs,  klcls, rad_lat,                         &
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
!     rad_lat   - latitudes in radians
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
  real,            intent(in), dimension(:,:)   :: rad_lat
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
integer :: klon, imax, jmax, kmax, klev, i,j,k
logical :: LDSHCV, lacc
logical :: l_lpi = .FALSE.
real, dimension(SIZE(temp0,1),SIZE(temp0,3)) :: plitot, zdgeoh, zdph,&
                                                ptenta,ptenqa,ptenrhoq,&
                                                ptenrhol,ptenrhoi,&
                                                ptenrhor,ptenrhos,&
                                                pmfu, pmfd,&
                                                pmfude_rate, pmfdde_rate,&
                                                ptu,pqu,plu,pcore
integer, dimension(SIZE(temp0,1),SIZE(temp0,2)) :: k650, k700, k950
real, dimension(SIZE(temp0,1),SIZE(temp0,2)) :: tropics_mask
real, dimension(SIZE(temp0,1),SIZE(temp0,2),SIZE(temp0,3)) :: rho_full
logical, dimension(SIZE(temp0,1)) :: ldcum
integer, dimension(SIZE(temp0,1)) :: ktype, kcbot, kctop, iseed
real,    dimension(SIZE(temp0,1)) :: fac_entrorg, fac_rmfdeps,&
                                     pcape, pvddraf, peis, cell_area,&
                                     mtnmask
real,    dimension(SIZE(zhalf0,1),SIZE(zhalf0,3)) :: pmflxr,pmflxs,pdtke_con
REAL, POINTER, DIMENSION(:) :: pertb => NULL ()
REAL, DIMENSION(2) :: mf_bulk,mf_perturb,mf_num ! setting dimension to silly size because won't use these
real :: abs_lat_degrees
real, parameter :: rad2deg = 57.2957795131
  ! --- Set dimensions
  imax  = size( temp0, 1 )
  jmax  = size( temp0, 2 )
  kmax  = size( temp0, 3 )
  klon = jmax
  klev = kmax
  ! Set this hardware switch to false
  lacc = .FALSE.
  ! Set these output variables
  ldcum = .FALSE.
  ktype = 0
  kctop = 0
  kcbot = 0
  peis  = 0.0
  ! Set mountain mask to 0 due to not knowing the required conversion.
  ! This means some conditional-on-orography behaviour will be switched off.
  mtnmask = 0.0
  ! Initialise shallow convection switch
  LDSHCV = .TRUE.
  ! Initialise ensemble perturbation config options
  fac_entrorg = 1.0
  fac_rmfdeps = 1.0
  ! Calculate air density
  rho_full(:,:,:) = pres0(:,:,:) / (rdgas * temp0(:,:,:))

  ! Find 650hPa and 700hPa levels
  DO i=1,imax
    DO j=1,jmax
      DO k=1,kmax
        IF (pres0(i,j,k) < 65000.0) THEN
          k650(i,j) = k
        ENDIF
        IF (pres0(i,j,k) < 70000.0) THEN
          k700(i,j) = k
        ENDIF
        IF (pres0(i,j,k) < 95000.0) THEN
          k950(i,j) = k
        ENDIF
      END DO
      k650(i,j) = MIN(k650(i,j),kmax-1)
      k700(i,j) = MIN(k700(i,j),kmax-1)
      k950(i,j) = MIN(k950(i,j),kmax-1)
    END DO
  END DO

  ! TODO: move this to an initialisation routine
  ! Set tropical mask
  DO i=1,imax
    DO j=1,jmax
      abs_lat_degrees = ABS(rad_lat(i,j)*rad2deg)
      IF (abs_lat_degrees < 25.0) THEN
        tropics_mask(i,j) = 1.0
      ELSE IF (abs_lat_degrees > 30.0) THEN
        tropics_mask(i,j) = 0.0
      ELSE
        tropics_mask(i,j) = (30.0 - abs_lat_degrees)/5.0
      ENDIF
    END DO
  END DO

  ! Begin main loop
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
    ! set cell area to negative so the model hopefully crashes if this is used
    cell_area = -1.0
    ! set iseed to 0 everywhere. This should not be used either.
    iseed = 0
    
    ! set total liquid and ice if present.
    IF (PRESENT(ql0) .and. PRESENT(qi0)) THEN
      plitot = ql0(:,j,:) + qi0(:,j,:)
    ELSE
      plitot = 0.0
    ENDIF

    ! set thickness of geopotential and pressure levels
    DO k=1,kmax
      zdgeoh(:,k) = ABS(zhalf0(:,k)-zhalf0(:,k+1))
      zdph(:,k) = ABS(pres0_int(:,k)-pres0_int(:,k+1))
    END DO

    ! Call the Bechtold scheme
    CALL cumastrn &
  & (  kidia,
  & kfdia, 
  & klon,   
  & ktdia,                              
  & klev,                                         &
  & ldland,                                       &
  & ldlake,                                       & 
  & ptsphy = dtime,                               & ! in: timestep (s)
  & phy_params = ,                                & ! in: TODO: add this input
  & k950 = k950,                                  & ! in: level at which 950hPa
  & trop_mask = tropics_mask,                     & ! in: 1 if lat less than 25, 0 if more than 30, and linear between.
  & mtnmask = mtnmask,                            & ! in: mountain mask
  & pten = temp0(:,j,:),                          & ! in: initial temperature
  & pqen = qvap0(:,j,:),                          & ! in: initial specific humidity
  & puen = uwnd0(:,j,:),                          & ! in: initial u wind
  & pven = vwnd0(:,j,:),                          & ! in: initial v wind
  & plitot = plitot,                              & ! in: initial total liquid and ice
  & pvervel = omega0(:,j,:),                      & ! in: initial vertical velocity
  & plen = ql0(:,j,:),                            & ! in: initial liquid
  & pien = qi0(:,j,:),                            & ! in: initial ice
  & shfl_s = flux_t0(:,j),                        & ! in: sensible heat flux
  & qhfl_s = flux_q0(:,j),                        & ! in: surface moisture flux
  & pqhfl = flux_q0(:,j),                         & ! in: surface moisture flux again (don't ask)
  & pahfs = flux_t0(:,j),                         & ! in: sensible heat flux again (don't ask)
  & pap = pres0(:,j,:),                           & ! in: pressure on full levels
  & paph = pres0_int(:,j,:),                      & ! in: pressure on half levels
  & pgeo = zfull0(:,j,:),                         & ! in: geopotential on full levels
  & pgeoh = zhalf0(:,j,:),                        & ! in: geopotential on half levels
  & zdph = zdph,                                  & ! in: pressure thickness on full levels
  & zdgeoh = zdgeoh,                              & ! in: geopotential thickness on full levels
  & ptent = dtent0(:,j,:),                        & ! out: temperature tendency
  & ptenu = duwnd0(:,j,:),                        & ! out: u wind tendency
  & ptenv = dvwnd0(:,j,:),                        & ! out: v wind tendency
  & ptenta = ptenta,                              & ! in: temperature tendency due to advection
  & ptenqa = ptenqa,                              & ! in: moisture tendency due to advection
  & ptenq = dqvap0(:,j,:),                        & ! out: moisture tendency kg/(kg*s)
  & ptenrhoq = ptenrhoq,                          & ! out: MOISTURE MASS DENSITY TENDENCY                KG/(M3*S)
  & ptenrhol = ptenrhol,                          & ! out: LIQUID WATER MASS DENSITY TENDENCY            KG/(M3*S)
  & ptenrhoi = ptenrhoi,                          & ! out: ICE CONDENSATE MASS DENSITY TENDENCY          KG/(M3*S)
  & ptenrhor = ptenrhor,                          & ! out: DETRAINED RAIN MASS DENSITY TENDENCY          KG/(M3*S)
  & ptenrhos = ptenrhos,                          & ! out: DETRAINED SNOW MASS DENSITY TENDENCY          KG/(M3*S)
  & ldcum = ldcum,                                & ! out: .TRUE. for convective points
  & ktype = ktype,                                & ! out: convection type (1: penetrative, 2: shallow, 3: midlevel)
  & kcbot = kcbot,                                & ! out: level of cloud bottom
  & kctop = kctop,                                & ! out: level of cloud top
  & LDSHCV=LDSHCV,                                & ! in: shallow convection indicator, not clear what this does
  & fac_entrorg = fac_entrorg,                    & ! in: tuning factor for ensemble perturbations for entrainment parameter
  & fac_rmfdeps = fac_rmfdeps,                    & ! in: tuning factor for ensemble perturbations for downdraft mass flux
  & pmfu = pmfu,                                  & ! out: mass flux updrafts
  & pmfd = pmfd,                                  & ! out: mass flux downdrafts
  & pmfude_rate,                                  & ! out: mass flux updraft detrainment rate
  & pmfdde_rate,                                  & ! out: mass flux downdraft detrainment rate
  & ptu = ptu,                                    & ! out: temperature in updrafts,
  & pqu = pqu,                                    & ! out: humidity in updrafts,
  & plu = plu,                                    & ! out: liquid in updrafts,
  & pcore = pcore,                                & ! out: updraft core fraction
  & pmflxr = pmflxr,                              & ! out: conv. rain flux
  & pmflxs = pmflxs,                              & ! out: conv. snow flux
  & prain = prain,                                & ! out: tot. prec. in updrafts without evap in downdrafts
  & pdtke_con = pdtke_con,                        & ! out: buoyant tke production at half levels
  & pcape = pcape,                                & ! out: cape
  & pvddraf = pvddraf,                            & ! out: convective gust at surface
  & peis = peis,                                  & ! out: estimated inversion strength
  & k650 = k650,                                  & ! in: level at which 650hPa
  & k700 = k700,                                  & ! in: level at which 700hPa
  & lacc = lacc                                   ) ! in: a hardware option we'll always have at .FALSE.

END DO
  end subroutine bechtold_alg
end module bechtold_alg_mod

