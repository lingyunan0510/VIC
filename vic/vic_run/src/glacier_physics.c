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
 * @return double Glacir Surafce Albedo         No Scale
 */
double calc_glacier_albedo(double min_abd, double snw_abd, double snw_dph) {
    double glc_abd;
    glc_abd = snw_abd*(1-exp(-snw_dph*MM_PER_M/24))+min_abd*exp(-snw_dph*MM_PER_M/24);
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
 * @brief 计算冰川能量平衡项
 * 
 * @param srf_tmp        Glacier Surface Temperature    C
 * @param air_tmp        Air Temperature                C
 * @param glc_abd        Glacier Minimum Albedo         No Scale
 * @param pcp_rat        Precipitation Rate             mm
 * @param shortwave_in   Income Shortwave               W/m^2
 * @param longwave_in    Income Longwave                W/m^2
 * @param air_den        Air Density                    kg/m^3
 * @param air_pre        Air Pressure                   Pa 
 * @param vap_pre        Vapor Pressure                 Pa 
 * @param dt 
 * @param Ra             Aerodynamics Resistance        No Scale
 * @return double 
 */
double calc_glacier_energy_balance( double srf_tmp, 
                                    double air_tmp, 
                                    double glc_abd,  
                                    double pcp_rat, 
                                    double shortwave_in, 
                                    double longwave_in, 
                                    double air_den, 
                                    double air_pre, 
                                    double vap_pre, 
                                    double dt,
                                    double Ra) {
    double netQ;

    double shortwave_net;
    double longwave_net;

    double sensible_heat;
    double latent_heat;
    double rain_heat;
    double ground_heat;
    
    // 净短波辐射
    shortwave_net = (1 - glc_abd) * shortwave_in;
    // 净长波辐射
    longwave_net = longwave_in - calc_glacier_outgoing_longwave_radiation(srf_tmp, 1);
    // 显热通量
    sensible_heat = calc_glacier_sensible_heat(air_den, air_tmp, srf_tmp, Ra);
    // 潜热通量
    latent_heat = calc_glacier_latent_heat(air_den, air_pre, srf_tmp, vap_pre, Ra);
    // 降水热通量
    rain_heat = calc_glacier_rain_heat(srf_tmp, air_tmp, pcp_rat, dt);
    // 地热通量
    ground_heat = 0.0;
    
    // 辐射余项为以上变量之和
    netQ = shortwave_net + longwave_net + sensible_heat + latent_heat + rain_heat + ground_heat;

    log_info("%f %f %f %f %f %f %f %f %f %f", 
    srf_tmp, 
    air_tmp, 
    glc_abd,  
    pcp_rat, 
    shortwave_in, 
    longwave_in,
    air_den, 
    air_pre, 
    vap_pre, 
    dt,
    Ra);


    return netQ;
}
