/**
 * @file glacier_physics.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief Calculate Glacier Energy Budget Items
 * @version 0.1
 * @date 2022-01-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <vic_run.h>

/**
 * @brief Calculate Sensible Heat Exchange Over Glacier Surface
 * Hock R, Noetzli C. Areal melt and discharge modelling of 
 * Storglaciären, Sweden[J]. Annals of glaciology, 1997, 24: 211-216.
 * 
 * @param air_den Air Density                   kg/m^3
 * @param air_tmp Air Temperature               C
 * @param srf_tmp Surface Temperature           C
 * @param Ra      Aerodynamics Resistance       No Scale
 * @return double SensibleHeat                  W/m^2
 */
double calc_glacier_sensible_heat(double air_den, double air_tmp, double srf_tmp, double Ra) {
    double sensible;
    sensible = CONST_CPMAIR * air_den * (air_tmp - srf_tmp) / Ra;
    return sensible;
}

/**
 * @brief Calculate Latent Heat Exchange Over Glacier Surface
 * 
 * @param air_den Air Density                   kg/m^3
 * @param air_pre Air Pressure                  Pa
 * @param srf_tmp Glacier Surface Temperature   C
 * @param vap_pre Vapor Pressure                Pa
 * @param Ra      Aerodynamics Resistance       No Scale
 * @return double LatentHeat                    W/m^2
 */
double calc_glacier_latent_heat(double air_den, double air_pre, double srf_tmp, double vap_pre, double Ra) {
    double latent = 0;
    latent = 0.623 * CONST_LATSUB * air_den * (vap_pre - svp(srf_tmp)) / air_pre / Ra;
    return latent;
}

/**
 * @brief Calculate outgoing Longwave Heat Flux
 * By Stefan-Boltzmann Law
 * @param srf_tmp Glacier Surface Temperature   C
 * @param emmis                                 Assumed as 1
 * @return double Outgoing Longwave Radiation   W/m^2
 */
double calc_glacier_outgoing_longwave_radiation(double srf_tmp, double emmis) {
    double longwave_radiation;
    double srf_tmp_k = srf_tmp + CONST_TKFRZ;
    longwave_radiation = emmis * CONST_STEBOL * srf_tmp_k * srf_tmp_k * srf_tmp_k * srf_tmp_k;
    return(longwave_radiation);
}

/**
 * @brief Calculate Heat from Rain
 * 
 * @param srf_tmp Glacier Surface Temperature   C
 * @param air_tmp Air Temperature               C
 * @param pcp_rat Precipitation Rate            m
 * @param dt      Step Time Delta               Second
 * @return double 
 */
double calc_glacier_rain_heat(double srf_tmp, double air_tmp, double pcp_rat, double dt) {
    double Qr;
    // the unit of precipitation should be transferred to M
    Qr = CONST_RHOFW*CONST_CPFW*(pcp_rat/MM_PER_M)*(air_tmp - srf_tmp)/dt;
    return(Qr);
}

/**
 * @brief Calculate Glacier Albedo Considering Snow Cover
 * 
 * @param min_abd Minimum Bare Ice albedo       No Scale
 * @param snw_abd Snow Surafce Albedo           No Scale
 * @param snw_dph Snow Depth                    m
 * @param d       Reference Depth               mm
 * @return double Glacir Surafce Albedo         No Scale
 */
double calc_glacier_albedo(double min_abd, double snw_abd, double snw_dph, double d) {
    double glc_abd;
    glc_abd = snw_abd*(1-exp(-snw_dph*MM_PER_M/d))+min_abd*exp(-snw_dph*MM_PER_M/d);
    // fprintf(LOG_DEST, "min_abd = %f\n", min_abd);
    // fprintf(LOG_DEST, "glc_abd = %f\n", glc_abd);
    return glc_abd;
}

double calc_glacier_ground_heat(int day) {
    double gh = 0.0;
    return gh;
}

double calc_saturated_vp_over_ice(temp) {
    double svp;
    return svp;
}

/**
 * @brief 占位函数 计算冰川能量平衡
 * 
 * @param srf_tmp    Glacier Surface Temperature    C
 * @param ...        

 * @return double 
 */
double calc_glacier_energy_balance(double srf_tmp, ...) {
    va_list ap;

    double netQ;
    va_start(ap, srf_tmp);

    netQ = glacier_energy_balance(srf_tmp, ap);

    va_end(ap);

    return netQ;
}

/**
 * @brief 实际函数 计算冰川能量平衡
 * 
 * @param srf_tmp    Glacier Surface Temperature    C
 * @param ap 
 * @return double 
 */
double glacier_energy_balance(double TSurf, va_list ap) {

    extern option_struct     options;
    extern parameters_struct param;

    /* Define Variable Argument List */

    /* General Model Parameters */
    double  Dt;                     /* Model time step (sec) */
    double  Ra;                     /* Aerodynamic resistance (s/m) */
    double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction */

    /* Vegetation Parameters */
    double  Z;                      /* Reference height (m) */
    double *Z0;                     /* surface roughness height (m) */

    /* Atmospheric Forcing Variables */
    double  AirDens;                /* Density of air (kg/m3) */
    double  EactAir;                /* Actual vapor pressure of air (Pa) */
    double  LongSnowIn;             /* Incoming longwave radiation (W/m2) */
    double  Lv;                     /* Latent heat of vaporization (J/kg3) */
    double  Press;                  /* Air pressure (Pa) */
    double  Rain;                   /* Rain fall (m/timestep) */
    double  NetShortUnder;          /* Net incident shortwave radiation (W/m2) */
    double  Vpd;                    /* Vapor pressure deficit (Pa) */
    double  Wind;                   /* Wind speed (m/s) */

    /* Snowpack Variables */
    double  OldTSurf;               /* Surface temperature during previous time step */
    double  SnowCoverFract;         /* Fraction of area covered by snow */
    double  SnowDepth;              /* Depth of snowpack (m) */
    double  SnowDensity;            /* Density of snowpack (kg/m^3) */
    double  SurfaceLiquidWater;     /* Liquid water in the surface layer (m) */
    double  SweSurfaceLayer;        /* Snow water equivalent in surface layer (m) */

    /* Energy Balance Components */
    double  Tair;                   /* Canopy air / Air temperature (C) */
    double  TGrnd;                  /* Ground surface temperature (C) */

    double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
    double *AdvectedSensibleHeat;   /* Sensible heat advected from snow-free area into snow covered area (W/m^2) */
    double *DeltaColdContent;       /* Change in cold content of surface layer (W/m2) */
    double *GroundFlux;             /* Ground Heat Flux (W/m2) */
    double *LatentHeat;             /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;          /* Latent heat of sublimation exchange at surface (W/m2) */
    double *NetLongUnder;           /* Net longwave radiation at snowpack surface (W/m^2) */
    double *RefreezeEnergy;         /* Refreeze energy (W/m2) */
    double *SensibleHeat;           /* Sensible heat exchange at surface (W/m2) */
    double *vapor_flux;             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
    double *blowing_flux;           /* Mass flux of water vapor from blowing snow. (m/timestep) */
    double *surface_flux;           /* Mass flux of water vapor from pack snow. (m/timestep) */

    /* Internal Routine Variables */

    double Density;                 /* Density of water/ice at TMean (kg/m3) */
    double NetRad;                  /* Net radiation exchange at surface (W/m2) */
    double RestTerm;                /* Rest term in surface energy balance (W/m2) */
    double TMean;                   /* Average temperature for time step (C) */
    double Tmp;                     /* Average temperature for time step (K) */
    double VaporMassFlux;           /* Mass flux of water vapor to or from the intercepted snow (kg/m2s) */
    double BlowingMassFlux;         /* Mass flux of water vapor from blowing snow. (kg/m2s) */
    double SurfaceMassFlux;         /* Mass flux of water vapor from pack snow. (kg/m2s) */

    /* General Model Parameters */
    Dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    Ra_used = (double *) va_arg(ap, double *);

    /* Vegetation Parameters */
    Z = (double) va_arg(ap, double);
    Z0 = (double *) va_arg(ap, double *);

    /* Atmospheric Forcing Variables */
    AirDens = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    LongSnowIn = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    NetShortUnder = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);

    /* Snowpack Variables */
    OldTSurf = (double) va_arg(ap, double);
    SnowCoverFract = (double) va_arg(ap, double);
    SnowDepth = (double) va_arg(ap, double);
    SnowDensity = (double) va_arg(ap, double);
    SurfaceLiquidWater = (double) va_arg(ap, double);
    SweSurfaceLayer = (double) va_arg(ap, double);

    /* Energy Balance Components */
    Tair = (double) va_arg(ap, double);
    TGrnd = (double) va_arg(ap, double);

    AdvectedEnergy = (double *) va_arg(ap, double *);
    AdvectedSensibleHeat = (double *)va_arg(ap, double *);
    DeltaColdContent = (double *) va_arg(ap, double *);
    GroundFlux = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    NetLongUnder = (double *) va_arg(ap, double *);
    RefreezeEnergy = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);
    blowing_flux = (double *) va_arg(ap, double *);
    surface_flux = (double *) va_arg(ap, double *);

    TMean = TSurf;
    Density = CONST_RHOFW;

    if (Wind > 0.0) {
        Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0[2]);
    } else {
        Ra_used[0] = param.HUGE_RESIST;
    }

    /* Calculate longwave exchange and net radiation */

    Tmp = TMean + CONST_TKFRZ;
    (*NetLongUnder) = LongSnowIn - calc_outgoing_longwave(Tmp, param.EMISS_SNOW);
    NetRad = NetShortUnder + (*NetLongUnder);

    // 显热
    *SensibleHeat = calc_sensible_heat(AirDens, Tair, TMean, Ra_used[0]);

    // 无雪地区输入显热
    (*AdvectedSensibleHeat) = 0.;

    /* Convert sublimation terms from m/timestep to kg/m2s */
    VaporMassFlux = *vapor_flux * Density / Dt;
    BlowingMassFlux = *blowing_flux * Density / Dt;
    SurfaceMassFlux = *surface_flux * Density / Dt;

    // 潜热 可能需要一部分修改
    latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra_used[0], TMean, Vpd, LatentHeat, LatentHeatSub, 
                          &VaporMassFlux,
                          &BlowingMassFlux,
                          &SurfaceMassFlux);

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * Dt / Density;
    *blowing_flux = BlowingMassFlux * Dt / Density;
    *surface_flux = SurfaceMassFlux * Dt / Density;

    // 降水热
    // 降水的温度即为气温
    *AdvectedEnergy = (CONST_CPFW * CONST_RHOFW * (Tair) * Rain) / Dt;

    // ColdContent的变化
    // 保留格式 但是冰川不需要
    *DeltaColdContent = 0.;

    // 地热 假设不存在
    /**
     * @todo
     * 赵求东方案
     */
    *GroundFlux = 0;

    // 计算冰雪一体表面的能量余项
    // 能量余项 = 净辐射 + 显热 + 潜热 + 升华潜热 + 降水热 + 地热 - Content变化 + 无雪地区输入显热
    // 其中 地热 无雪地区输入显热 被标记为0. 仅保留形式防止程序崩溃
    RestTerm = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy + *GroundFlux - *DeltaColdContent + *AdvectedSensibleHeat;

    // 基于冰雪一体的表面存在的液态水
    // 但是冰川不考虑
    *RefreezeEnergy = 0.;

    /**
     * @attention
     * 直接返回余项
     * 余项为正值 说明冰川表面获得能量
     * 余项为负值 说明冰川表面失去能量
     */
    return RestTerm;
}