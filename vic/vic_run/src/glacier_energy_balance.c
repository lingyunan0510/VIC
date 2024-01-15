// /**
//  * @file glacier_energy_balance.c
//  * @author Yunan Ling (lingyunan@outlook.com)
//  * @brief 
//  * @version 0.1
//  * @date 2022-01-25
//  * 
//  * @copyright Copyright (c) 2022
//  * 
//  */
// #include <vic_run.h>

// /**
//  * @brief Calculate Energy Balance of Glacier Surface
//  * 
//  * @param TSurf 
//  * @param ap 
//  * @return double 
//  */
// double glacier_energy_balance(double TSurf, va_list ap) {
//     // option setting file
//     extern option_struct options;
//     // 
//     extern parameters_struct param;

//     /* Define Variable Argument List */

//     /* General Model Parameters */
//     double  Dt;                     /* Model time step (sec) */
//     double  Ra;                     /* Aerodynamic resistance (s/m) Used */
//     double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction Used*/

//     /* Vegetation Parameters */
//     double  Z;                      /* Reference height (m) used */
//     double *Z0;                     /* surface roughness height (m) used */

//     /* Atmospheric Forcing Variables */
//     double  AirDens;                /* Density of air (kg/m3) used */
//     double  EactAir;                /* Actual vapor pressure of air (Pa) used */
//     double  LongSnowIn;             /* Incoming longwave radiation (W/m2) used */
//     double  Lv;                     /* Latent heat of vaporization (J/kg3) */
//     double  Press;                  /* Air pressure (Pa) used */
//     double  Rain;                   /* Rain fall (m/timestep) used */
//     double  NetShortUnder;          /* Net incident shortwave radiation (W/m2) used */
//     double  Wind;                   /* Wind speed (m/s) used */

//     /* Snowpack Variables */
//     double  OldTSurf;               /* Surface temperature during previous time step */
//     double  SnowCoverFract;         /* Fraction of area covered by snow */
//     double  SnowDepth;              /* Depth of snowpack (m) */
//     double  SnowDensity;            /* Density of snowpack (kg/m^3) */
//     double  SurfaceLiquidWater;     /* Liquid water in the surface layer (m) */
//     double  SweSurfaceLayer;        /* Snow water equivalent in surface layer (m) */

//     /* Glacier Variables */

//     /* Energy Balance Components */
//     double  Tair;                   /* Canopy air / Air temperature (C) used */
//     double  TGrnd;                  /* Ground surface temperature (C) used */

//     double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
//     double *AdvectedSensibleHeat;   /* Sensible heat advected from snow-free area into snow covered area (W/m^2) */
//     double *DeltaColdContent;       /* Change in cold content of surface layer (W/m2) */
//     double *GroundFlux;             /* Ground Heat Flux (W/m2) */
//     double *LatentHeat;             /* Latent heat exchange at surface (W/m2) */
//     double *LatentHeatSub;          /* Latent heat of sublimation exchange at surface (W/m2) */
//     double *NetLongUnder;           /* Net longwave radiation at snowpack surface (W/m^2) */
//     double *RefreezeEnergy;         /* Refreeze energy (W/m2) */
//     double *SensibleHeat;           /* Sensible heat exchange at surface (W/m2) */
//     double *vapor_flux;             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
//     double *blowing_flux;           /* Mass flux of water vapor from blowing snow. (m/timestep) */
//     double *surface_flux;           /* Mass flux of water vapor from pack snow. (m/timestep) */

//     /* Internal Routine Variables */

//     double Density;                 /* Density of water/ice at TMean (kg/m3) */
//     double NetRad;                  /* Net radiation exchange at surface (W/m2) */
//     double RestTerm;                /* Rest term in surface energy balance (W/m2) */
//     double TMean;                   /* Average temperature for time step (C) */
//     double Tmp;                     /* Average temperature for time step (K) */
//     double VaporMassFlux;           /* Mass flux of water vapor to or from the intercepted snow (kg/m2s) */
//     double BlowingMassFlux;         /* Mass flux of water vapor from blowing snow. (kg/m2s) */
//     double SurfaceMassFlux;         /* Mass flux of water vapor from pack snow. (kg/m2s) */

//     double GlacierAlbedo;           /* Glacier Surafce Albedo (No Scale) */
//     double SnowAlbedo;              /* Snow Surafce Albedo (No Scale)*/

//     /* General Model Parameters */
//     Dt = (double) va_arg(ap, double);
//     Ra = (double) va_arg(ap, double);
//     Ra_used = (double *) va_arg(ap, double *);

//     /* Vegetation Parameters */
//     Z = (double) va_arg(ap, double);
//     Z0 = (double *) va_arg(ap, double *);

//     /* Atmospheric Forcing Variables */
//     AirDens = (double) va_arg(ap, double);
//     EactAir = (double) va_arg(ap, double);
//     LongSnowIn = (double) va_arg(ap, double);
//     Lv = (double) va_arg(ap, double);
//     Press = (double) va_arg(ap, double);
//     Rain = (double) va_arg(ap, double);
//     NetShortUnder = (double) va_arg(ap, double);
//     Wind = (double) va_arg(ap, double);

//     /* Snowpack Variables */
//     OldTSurf = (double) va_arg(ap, double);
//     SnowCoverFract = (double) va_arg(ap, double);
//     SnowDepth = (double) va_arg(ap, double);
//     SnowDensity = (double) va_arg(ap, double);
//     SurfaceLiquidWater = (double) va_arg(ap, double);
//     SweSurfaceLayer = (double) va_arg(ap, double);

//     /* Energy Balance Components */
//     Tair = (double) va_arg(ap, double);
//     TGrnd = (double) va_arg(ap, double);

//     AdvectedEnergy = (double *) va_arg(ap, double *);
//     AdvectedSensibleHeat = (double *)va_arg(ap, double *);
//     DeltaColdContent = (double *) va_arg(ap, double *);
//     GroundFlux = (double *) va_arg(ap, double *);
//     LatentHeat = (double *) va_arg(ap, double *);
//     LatentHeatSub = (double *) va_arg(ap, double *);
//     NetLongUnder = (double *) va_arg(ap, double *);
//     RefreezeEnergy = (double *) va_arg(ap, double *);
//     SensibleHeat = (double *) va_arg(ap, double *);

//     vapor_flux = (double *) va_arg(ap, double *);
//     blowing_flux = (double *) va_arg(ap, double *);
//     surface_flux = (double *) va_arg(ap, double *);

//     /* Calculate active temp for energy balance as average of old and new  */

//     TMean = TSurf;
//     Density = CONST_RHOFW;

//     /* 校正空气动力系数 */
//     if (Wind > 0.0) {
//         Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0[2]);
//     } else {
//         Ra_used[0] = param.HUGE_RESIST;
//     }

//     /* 计算当前温度 */
//     Tmp = TMean + CONST_TKFRZ;
//     /* 净长波辐射 = 入射长波辐射 -出射长波辐射 */
//     (*NetLongUnder) = LongSnowIn - calc_glacier_outgoing_longwave_radiation(TSurf, 1.0);
//     /* 净辐射 = 净入射短波辐射 + 净入射长波辐射 */
//     NetRad = NetShortUnder + (*NetLongUnder);

//     /* 显热 */
//     *SensibleHeat = calc_glacier_sensible_heat(AirDens, Tair, TMean, Ra_used[0]);

//     /* 潜热 */
//     *LatentHeat = calc_glacier_latent_heat(AirDens, Press, TMean, EactAir, Ra_used[0]);

//     /* 地面 */
//     *GroundFlux = 0.0;

//     /* 降水 */
//     if (TMean >= 0.) {
//         *AdvectedEnergy = (CONST_CPFW * CONST_RHOFW * (Tair) * Rain) / Dt;
//     } else {
//         *AdvectedEnergy = 0;
//     }

//     /* 辐射余项 = 净辐射 + 显热 + 潜热 + 地面热通量 + 降水热量*/
//     RestTerm = NetRad + *SensibleHeat + *LatentHeat + *GroundFlux + *AdvectedEnergy;

//     *RefreezeEnergy = 0.0;

//     /**
//      * @brief 
//      * 由于需要采用root brent法计算
//      * 本代码基础结构参照snow/ice_energy_balance
//      */
//     // 仅在表面温度为0 且 能量平衡余项为正时 发生融化
//     if (TSurf == 0.0 && RestTerm > (*RefreezeEnergy)) {
//         *RefreezeEnergy = -RestTerm;
//         RestTerm = 0.0;
//     } else {
//         // 有两种情况
//         // 表面温度低于0 即冰川未达到熔化温度 即使有能量输入仅使冰川升温
//         // 或能量平衡余项为负 即冰川变冷 冰川停止融化或降温
//         RestTerm += *RefreezeEnergy;
//     }

//     /**
//      * @brief 计算之后需要注意的点
//      * 能量余项RestTerm在输出后作为指示变量
//      * 取值为0时意为融化 但RefreezeEnergy才真正表示可用能量
//      * 
//      * 在非0值情况下 RestTerm可作为冰川升温/降温的可用能量
//      * 但是与融化过程无关s
//      */

//     return RestTerm;
// }