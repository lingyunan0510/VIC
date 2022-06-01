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
 * @param air_den Atmospheric Density
 * @param air_tmp Atmospheric Temperature
 * @param srf_tmp Surface Temperature
 * @param Ra 
 * @return double 
 */
double calc_glacier_sensible_heat(double air_den, double air_tmp, double srf_tmp, double Ra) {
    double sensible;
    sensible = CONST_CPMAIR * air_den * (air_tmp - srf_tmp) / Ra;
    return sensible;
}

/**
 * @brief Calculate Latent Heat Exchange Over Glacier Surface
 * 
 * @param air_den 
 * @param air_pre 
 * @param srf_tmp 
 * @param vap_pre 
 * @param Ra 
 * @return double 
 */
double calc_glacier_latent_heat(double air_den, double air_pre, double srf_tmp, double vap_pre, double Ra) {
    double latent = 0;
    latent = 0.623 * CONST_LATSUB * air_den * (vap_pre - svp(srf_tmp)) / air_pre / Ra;
    return latent;
}

/**
 * @brief Calculate outgoing Longwave Heat Flux
 * By Stefan-Boltzmann Law
 * @param srf_tmp Glacier Surface Temperature C
 * @param emmis   Assumed as 1
 * @return double 
 */
double calc_glacier_outgoing_longwave_radiation(double srf_tmp, double emmis) {
    double longwave_radiation;
    double srf_tmp_k = srf_tmp + CONST_TKFRZ;
    longwave_radiation = emmis * CONST_STEBOL * srf_tmp_k * srf_tmp_k * srf_tmp_k * srf_tmp_k;
    return(longwave_radiation);
}

/**
 * @brief hao
 * 
 * @param srf_tmp 
 * @param air_tmp 
 * @param pcp_rat 
 * @param dt 
 * @return double 
 */
double calc_glacier_rain_heat(double srf_tmp, double air_tmp, double pcp_rat, double dt) {
    double Qr;
    Qr = CONST_RHOFW*CONST_CPFW*(pcp_rat/MM_PER_M)*(air_tmp - srf_tmp)/dt;
    return(Qr);
}

/**
 * @brief Calculate Glacier Albedo Considering Snow Cover
 * 
 * @param air_tmp 2-Meter Air Temperature C 
 * @param snw_abd Snow Surafce Albedo 
 * @param snw_dph Snow Depth mm 
 * @return double Glacir Surafce Albedo
 */
double calc_glacier_albedo(double air_tmp, double snw_abd, double snw_dph) {
    double ai;
    double gi;
    if (air_tmp <= -15.0) {
        air_tmp = -15.0;
    }
    if (air_tmp > 10.0) {
        air_tmp = 10.0;
    }
    ai = 0.324-0.018*air_tmp;
    // double tdew;
    // tdew = (116.91+237.3*log(air_tmp*(1/PA_PER_KPA)))/(16.78-log(air_tmp*(1/PA_PER_KPA)));
    // ai = 0.2577+(-0.031)*tdew;
    gi = snw_abd*(1-exp(-snw_dph*MM_PER_M/24))+ai*exp(-snw_dph*MM_PER_M/24);
    return gi;
}

double calc_glacier_ground_heat(int day) {
    double gh = 0.0;
    return gh;
}

double calc_saturated_vp_over_ice(temp) {
    double svp;
    return svp;
}
