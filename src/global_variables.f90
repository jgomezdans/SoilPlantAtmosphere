! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module gv_clim

  implicit none

  ! met-drivers used/updated each timestep..
  real :: atmos_press, & ! Surface atmospheric pressure (Pa)
          coa,         & ! Ambient CO2 concentration (ppm)
          par_top,     & ! PAR at top canopy layer
          ppt,         & ! Precipitation (mm)
          sw_rad,      & ! Incident short wave radiation (Wm-2)
          temp_bot,    & ! Surface temperature at bottom and..
          temp_top,    & !  ..at top of canopy
          vpd_bot,     & ! Vapour Pressure Deficit at bottom..
          vpd_top,     & !  ..and at top of canopy.
          wind_spd       ! Wind speed (m/s)

  ! other meteorology-related variables.. 
  real    :: avtemp,     & ! average daily temperature (Celcius); used in deciduous
             daypar,     & ! accumulation of PAR over a day (umol.m-2.s-1)
             dayppt,     & ! accumulation of precipitation over a day (mm.d-1)
             daytempsum, & ! daily-sum of temperatures (Celcius)
             gdd,        & ! growing degree days for decidous forests (Celcius) 
             max_fol,    & ! deciduous switch to indicate that maxfoliage has been reached
             mint,       & ! minimum daily temperature
             rnet,       & ! ?
             wdbot,      & ! wind speed at bottom of canopy (m.s-1)
             wdtop         ! wind speed at top of canopy (m.s-1)

  real,parameter :: par_ratio = 4.6  ! PAR energy convertion ratio with SW radiation (umol.J-1)

  real,dimension(:),allocatable :: &
                    gbw, & ! boundary layer conductance for each canopy layer (m.s-1)
                  wetev    ! Wet evaporation from canopy (mm.t-1) per time step

end module gv_clim
!
!------------------------------------------------------------------------
!
module gv_hourscale

  implicit none

  integer :: hour            ! ?
  real    :: canopy_store = 0, & ! ?
                  dayevap, & ! soil moisture flux '=-Qe'
                discharge, & ! ?
               evap_store, & ! ?
                      gaw, & ! boundary layer conductance (m.s-1), calculated each time step 
                      gws, & ! soil conductance to water vapour 
                  hourppt, & ! Met -> precip
                hourpress, & ! Met -> surface pressure
                  hourrad, & ! Soil net radiation balance including both long and shortwave (W.m-2),
                             !  used in long wave determination.
                 hourrnet, & ! Met -> sw rad
                 hourtemp, & ! Met -> temperature (oC) converted to (K)
                 hourtime, & ! ?
                   hourts, & ! ?
                  hourvpd, & ! Met -> VPD
                 hourwind, & ! Met -> wind speed (m.s-1)
                 overflow, & ! over flow from water input into the soil layer. Set as a fixed 
                             !  proportion of surface_watermm (mm)
                       Qc, & ! ?
                       Qe, & ! Latent energy flux (soil) within the model structure (W.m-2)
                       Qh, & ! ?
                       Qn, & ! Net LW emissions (W.m-2) - based upon the emissivity and the total radation penetration
                   runoff, & ! ?
          surface_watermm, & ! timestep calculated surface water. i.e. the water that cannot be
                             !  infiltrated within the given timestep or surface layer saturation (mm) 
                    totet, & ! ?
                underflow, & ! ?
            unintercepted    ! ?

  real,parameter :: freeze = 273.15 ! freezing point of water in Kelvin

end module gv_hourscale
!
!------------------------------------------------------------------------
!
module gv_Hydrol

  implicit none

  real,dimension(:),allocatable :: &
                soil_frac_clay, & ! Percentage of soil that is clay.
                soil_frac_sand    ! Percentage of soil that is sand.

end module gv_Hydrol
!
!------------------------------------------------------------------------
!
module gv_Irradiance_Sunshade

  implicit none

  real :: check,  & !
          soilnet   ! absorbed radiation by soil (W.m-2); no long emissions, isothermal or otherwise

end module gv_Irradiance_Sunshade
!
!------------------------------------------------------------------------
!
module gv_metab

  implicit none

  real :: an, & ! GPP for given canopy layer in the shade or light leaf area loops (umolC.m-2.timestep-1)
          ci, & ! intra leaf CO2 concentration (ppm)
          ht, & ! Height of specific canopy layer, used in canopy loop (m)
 layer_capac, & ! Canopy layer specific canopy capacitance, based on canopy area
          rn, & ! Respiration constant for leaf area under photosynthetic analysis (umol CO2.gN.m-2 at 10oC)
      rplant, & ! Plant hydrolic resistance for each canopy layer (MPa.s-1.mmol-1)
       rsoil, & ! Root and soil resistence for a given canopoy layer (MPa.s-1.mmol-1)
         vcm, & ! Maximum rate of carboxylation (umol CO2.gN.m-2.s-1); based on kappaV coefficient
         vjm    ! Maximum rate of electrob transport (umol.m-2.s-1)

  real,parameter :: metabolic_opt_temp = 30.0 ! metabolic temperature optimum

end module gv_metab
!
!------------------------------------------------------------------------
!
module gv_meteo

  implicit none

  real :: gbb, & ! boundary layer conductance of specific canopy layer in spa_canopy.F loop (m3.s-1)
           la, & ! Leaf area in current photosynthetic pass, i.e. sun or shade leaf areas (m.canopy_layer-1))
          nit, & ! Canopy layer specific nitrogen content (gN.canopy_layer-1)
          par, & ! PAR penetrating given canopy layer in spa_canopy.F loop (umol.m-2.s-1)
         psil, & ! canopy layer specific leaf water potential in spa_canopy.F loop (MPa)
         psis, & ! Weighted soil water potential, by total evaporation estimate (MPa)
          rad, & ! net radiation penetration to canopy layer (kW.m-2)
         temp, & ! temperature of a given canopy layer (oC); see spa_canopy.F
         wdef    ! wind speed at given canopy layer within spa_canopy loop (m.s-1)

  real,parameter :: gi   = 1.0,      & ! mesophyll conductance (m s-1) set high, and thus ignored, /0.0025/ old value
                                       ! v large mesophyll conductance - ACi curves have not been Epron adjusted
                    head = 0.009807, & ! head of pressure  (MPa/m)
                    Rcon = 8.3144      ! gas constant

end module gv_meteo
!
!------------------------------------------------------------------------
!
module gv_scale_declarations

  !! These declarations define the size of the soil  !!
  !!  and canopy profiles, & the time resolution.    !!
  !! They also declare constructs that are needed at !!
  !!  all levels, such as the user's config options, !!
  !!  the time holder, and the met-drivers.          !!
  !! Physical parameters, such as pi, are also here. !!

  implicit none

  ! Internal model lengths (fixed)..
  integer,parameter :: fname_length = 100, & ! length of filename variables
                       nmax         = 2      ! number of iterations in math-loops


  ! Model grid sizing info
  type grid_holder
    integer :: canopy    = 10,    & ! number of canopy layers
               core      = 21,    & ! number of soil layers + 1
               soil      = 20,    & ! number of soil layers
               wetting   = 10       ! number of layers to use for wetting calcs
    real    :: latitude  = 50.00, & ! +ve == nth, -ve == sth
               longitude = 00.00    ! 0--360. -ve not recognised.
  end type
  type(grid_holder),save :: grid


  ! All the time information in one place..
  type time_holder
    integer :: nos_of_years  = 1,   & ! number of years to simulate
               days_in_year  = 365, & ! number of (whole!) days per year
               steps_per_day = 24,  & ! number of timesteps per day
               year = 0, & ! Current year
               day  = 0, & ! Current day
               step = 0, & ! Current step
               steps_count = 0             ! count of number of steps completed so far.
    real    :: seconds_per_step = 3600., & ! timesteps (3600=>fixed at 24 steps per day)
               daytime = 0.                ! date+time in days, e.g. noon on day 30 of yr 2
                                           ! (if yr=365days) == 365 + 30 + 0.5 = 395.5
  end type time_holder
  type(time_holder),save :: time


  ! meteorological inputs/drivers
  type met_drivers
    real,allocatable,dimension(:) :: ambient_co2  ! (ppm) Ambient Carbon Dioxide concentration
    real,allocatable,dimension(:) :: par          ! (umol.m-2.s-1) Photosynthetically active radiation
    real,allocatable,dimension(:) :: precip       ! (mm.t-1) precipitation
    real,allocatable,dimension(:) :: sfc_pressure ! (pa)  Atmospheric surface pressure
    real,allocatable,dimension(:) :: sw_rad       ! (W.m-2) surface downward short-wave radiation
    real,allocatable,dimension(:) :: temp         ! (C) temperature
    real,allocatable,dimension(:) :: vpd          ! (kPa) vapour pressure deficit
    real,allocatable,dimension(:) :: wind_spd     ! (ms-1) wind strength
    ! If user does not provide variables, but still wants to alter co2/pressure from defaults
    ! then they can supply these in the 'met parameters' section of the config file..
    real :: const_ambient_co2  = 360    ! (ppm) Ambient Carbon Dioxide concentration
    real :: const_sfc_pressure = 100000 ! (pa)  Atmospheric surface pressure
  end type met_drivers
  type(met_drivers),save :: met


  ! All of the configuration info for SPA...
  type user_config_holder
     ! Input files --------------------------------------------------------
     character(fname_length) :: met_filename     = ''  ! Met file name/path
     logical                 :: met_file_is_nc = .false. ! Is met file in netcdf format?
     character(fname_length) :: soils_filename   = ''  ! Soils file name/path
     character(fname_length) :: veg_filename     = ''  ! Veg file name/path
     character(fname_length) :: restart_filename = ''  ! file to use for restart
     logical                 :: load_restart = .false. ! Should SPA start from a restart dump?
     ! Options ----------------------------------------------------------
     logical :: veg_is_deciduous  = .False.
     logical :: loop_met_driver   = .False. ! repeat the met-driver as needed?
     logical :: use_co2_from_met_file  = .False. ! Whether to use the CO2 values in the input met file
     logical :: use_par_from_met_file  = .False. ! likewise for PAR.
     logical :: met_file_has_sfc_press = .False. ! likewise for surface atmospheric pressure
     logical :: precip_is_rate         = .False. ! Is precip in rate (mm/sec) or volume-per-timestep (mm)?
     ! Output ----------------------------------------------------------
     integer :: print_msg_at_severity = 5 ! print-to-screen messages of this severity or higher
     character(fname_length) :: output_directory =  ''    ! directory to put output data in
     logical :: std_csv_output          = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: make_restart_file       = .false. ! Should SPA make a restart dump?
     integer :: restart_write_frequency = 1000    ! frequency (in years) with which to produce restart dumps
  end type user_config_holder
  type(user_config_holder),save :: user_opts


  ! General physical constants..
  real,parameter :: boltz = 5.670400e-8, & ! (W m-2 K-4)
                    pi    = 3.14159265     ! (-)


  ! Internal store of what version of the model this is,
  ! where (23,0) equates to "2.3.0"
  integer,parameter :: model_sci_version = 24, & ! 20=>2.0, 23=> 2.3, etc.
                       model_bug_version = 0     ! 0=>0, 1=>1, etc.  


end module gv_scale_declarations
!
!------------------------------------------------------------------------
!
module gv_soil_structure

  implicit none

  real,parameter :: abovebelow = 1. ! saxton water retention equation are off by default

  integer :: rooted_layers     ! number of root layers penetrated into the soil layers

  real ::   drythick = 0.1,  & ! Thickness of dry soil layer above water table (m); minimum level is set to 0.001 (m)
                 max_depth,  & ! PARAMETER maximum rooting depth (m)
               max_storage,  & ! PARAMETER maximum canopy water storage (mm)
      resp_rate_temp_coeff,  & ! response of autotrophs to change in temperature
              root_biomass,  & ! root biomass (g Biomass.m-2) in total is determined as twice the root C (gC.m-2)
               root_radius,  & ! PARAMETER root radius for specific land cover type
                root_reach,  & ! maximum depth of root system based on the available biomass. i.e. does not have to be the
                               !   max_rooting depth parameter
                    root_k,  & ! PARAMETER mass of roots for reaching 50% maximum depth (g.m-2!)
                      snow,  & ! Used in calculation of change in soil profile wetting space
              surf_biomass,  & ! calculation of surface biomass (g), assumes that 50 % of biomass is in top 25 % of the soil profile
              through_fall,  & ! PARAMETER fraction of precip as throughfall 
              weighted_SWP     ! Soil water potential for whole soil profile, weighted by the total evaporation. This value is 
                               !   determined each time step for use in spa_canopy

  real :: rootresist = 400., & ! number of m of roots in a layer to give resistance of 1 MPa m2 s mmol-1
          thermal    = 1.58    ! conversion rate of thermal conductivity from (W m-1 K-1) to (J m-1 K-1 h-1)

  ! Key components of the Saxton soil water retention equations, these are re calculated on each timestep
  real,dimension(:),allocatable :: cond1, cond2, cond3, potA, potB

  real,dimension(:),allocatable :: &
                          conduc,  & ! Soil layer conductivity (m.s-1)
                 field_capacity,   & ! Field capacity of moisture for each layer (mm?), when soil water content at SWP = -10kPa
                fraction_uptake,   & ! fraction of evapotranspiration (i.e. root draw) from each layer. Reset at each timestep
                        iceprop,   & ! Soil ice proportion 
                    layer_depth,   & ! PARAMETER Soil layer depth (m)
                    mineralfrac,   & ! PARAMETER soil proportional inorganic content
                    organicfrac,   & ! PARAMETER soil proportional organic content
                       porosity,   & ! soil layer porosity
                        pptgain,   & ! soil water potential (MPa); used in determining vapour pressure (kPa) of soil air space
                                     !   for latent heat fluxes and water movement between between layers.
                    root_length,   & ! (m) based on root biomass and an exponential decline of mass with soil depth.
                                     !   Assumes root subroutine called each timestep.
                      root_mass,   & ! (g biomass) per soil layer are determined every timestep in root call.
                      soil_temp,   & !
               soil_temp_nplus1,   & ! temporary soil temperature (K) vector used for determining the next time step temperature
                                     !   soil profile beginning with the new surface temperature ure and the constant soil core value.
                          soilR,   & ! soil root hydraulics resistance (MPa.s-1.m-2.mmol-1); calculated from soilR1 and SoilR2 
                         soilR1,   & ! soil root hydraulics component; represents some function of the root mass and volume
                         soilR2,   & ! soil root hydraulics component; represents some function of root mass and the root resistance
                                     !   and ratio of mass above and below ground
                            SWP,   & ! Soil water potential (MPa); used in determining vapour pressure (kPa) of soil air space for
                                     !   latent heat fluxes and water movement between between layers.
                      thickness,   & ! PARAMETER Soil layer thickness (m)
                      waterfrac,   & ! Extracts soil moisture by volumetric mixing ratio (m3.m-3) from its array  
                      watergain,   & ! Water gained (mm) from a given soil layer; reset with each timestep
                      watericemm,  & ! Total soil moisture content per layer (mm)
                      waterloss,   & ! Water lost (mm) from a given soil layer; reset with each timestep
                      wettingbot,  & ! Depth to bottom of wet soil layers (m)
                      wettingtop     ! Depth to top of wet soil layers (m)

end module gv_soil_structure
!
!------------------------------------------------------------------------
!
module gv_Snow_Info

  implicit none

  real :: snow_watermm,  & ! 
          snowheight,    & ! 
          snowweight       !

end module gv_Snow_Info
!
!------------------------------------------------------------------------
!
module gv_veg

  implicit none

  integer :: nlink        = 1, &  !
             conductivity = 1     ! 1=>conductivity set, 0=>conductance set

  real    :: altitude, & !
                       avN,  & !
             canopy_height,  & ! canopy height (m) is land cover type specific, generated from top canopy layer
                     capac,  & !
                    co2amb,  & !
               dimen = 0.08, & ! HF leaf dimension
                    gplant,  & ! depends on conductivity switch. if 0 => plant hydraulic conductance  (mmol m-2 s-1 MPa-1),
					           !                                 if 1 => plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      iota,  & !
                    kappac,  & !
                    kappaj,  & !
                       lat,  & ! degrees input converted to radians in io
                       LMA,  & ! leaf carbon per unit area (g.m-2)
                    minlwp,  & !
                   modRnet,  & !
                     prevC,  & !
                     sensh,  & !
                    totass,  & !
                   totevap,  & !
                     totla,  & !
                      totn,  & !
                    totres,  & !
              tower_height     ! height of tower at which observation are made (m)

  ! Carbon pools (gC m-2)...
  real :: stock_foliage,       & ! Carbon stocks of..foliage pool
          stock_labile,        & ! ..labile pool
          stock_litter,        & ! ..litter pool
          stock_resp_auto,     & ! ..autotrophic respiration pool
          stock_roots,         & ! ..roots pool
          stock_soilOrgMatter, & ! ..soil organic matter pool
          stock_stem             ! ..stem pool (ie the trunk/wood/stem on a plant)

  ! Fluxes to/from carbon pools..
  real :: alloc_to_foliage,    & ! Carbon allocation..to foliage pool
          alloc_to_labile,     & ! ..to labile pool
          alloc_from_labile,   & ! ..from labile pool
          alloc_to_roots,      & ! ..to roots pool
          alloc_to_stem,       & ! ..to structure pool
          decomposition,       & ! Decomposition
          GPP,                 & ! Gross primary product
          litterfall_foliage,  & ! Loss of stock from...foliage
          litterfall_stem,     & ! ..stem
          litterfall_roots,    & ! ..roots
          resp_auto,           & ! Respiration - autotrophic
          resp_h_litter,       & ! Respiration - heterotrophic, in litter pool
          resp_h_soilOrgMatter   ! Respiration - heterotrophic, in soil organic matter pool

  ! Rates/fractions/thresholds that affect carbon pool calculations..
  real :: decomposition_rate,                & ! decomposition rate
          frac_alloc_foliage,                & ! fraction of npp allocated..to foliage
          frac_alloc_roots,                  & ! ..to fine roots
          frac_GPP_resp_auto,                & ! fraction of gpp respired
          frac_leafLoss_litter,              & ! fraction of leaf loss to litter
          GDD_thresh,                        & ! GDD threshhold
          max_stock_foliage,                 & ! maximum foliar carbon stock
          min_Temp_thresh,                   & ! minimum temperature threshold
          mineralisation_rate_litter,        & ! Rate of mineralisation of..litter pool
          mineralisation_rate_soilOrgMatter, & ! ..soil organic matter pool
          resp_cost_labile_trans,            & ! respiratory cost of labile transfers
          turnover_rate_foliage,             & ! turnover rate of..foliage pool
          turnover_rate_labile,              & ! ..labile pool
          turnover_rate_resp_auto,           & ! ..autotrophic respiration pool
          turnover_rate_roots,               & ! ..roots pool
          turnover_rate_stem                   ! ..stem pool

  ! Arrays..
  real,dimension(:),allocatable ::  &
                canopy_soil_resistance, & !
                                   ess, & !
                                  gppt, & !
                                lafrac, & !
                                   lai, & !
                          layer_height, & !
                              LWPstore, & ! initial LWP
                                 nfrac, & !
                                   Nla, & !
                                 respt, & !
                              soiletmm, & !
                                transt    !

  real,allocatable :: flux(:,:) !

end module gv_veg
!
!------------------------------------------------------------------------
!
