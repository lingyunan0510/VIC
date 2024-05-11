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
 * @param ...        double air_tmp, double glc_abd, double pcp_rat, double shortwave_in, double longwave_in, double air_den, double air_pre, double vap_pre, double dt, double Ra,

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

// /**
//  * @brief 计算冰川能量平衡项
//  * 
//  * @param srf_tmp        Glacier Surface Temperature    C
//  * @param air_tmp        Air Temperature                C
//  * @param glc_abd        Glacier Minimum Albedo         No Scale
//  * @param pcp_rat        Precipitation Rate             mm
//  * @param shortwave_in   Income Shortwave               W/m^2
//  * @param longwave_in    Income Longwave                W/m^2
//  * @param air_den        Air Density                    kg/m^3
//  * @param air_pre        Air Pressure                   Pa 
//  * @param vap_pre        Vapor Pressure                 Pa 
//  * @param dt 
//  * @param Ra             Aerodynamics Resistance        No Scale
//  * @return double 
//  */

/**
 * @brief 实际函数 计算冰川能量平衡
 * 
 * @param srf_tmp    Glacier Surface Temperature    C
 * @param ap 
 * @return double 
 */
double glacier_energy_balance(double srf_tmp, va_list ap) {
    
    double air_tmp; 
    double glc_abd; 
    double pcp_rat; 
    double shortwave_in; 
    double longwave_in; 
    double air_den; 
    double air_pre; 
    double vap_pre; 
    double dt; 
    double Ra; 
    // bool is_snow;

    air_tmp = (double) va_arg(ap, double);
    glc_abd = (double) va_arg(ap, double);
    pcp_rat = (double) va_arg(ap, double);
    shortwave_in = (double) va_arg(ap, double);
    longwave_in = (double) va_arg(ap, double);
    air_den = (double) va_arg(ap, double);
    air_pre = (double) va_arg(ap, double);
    vap_pre = (double) va_arg(ap, double);
    dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    // is_snow = (bool) va_arg(ap, bool);

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
    rain_heat = 0.0;
    // if (is_snow) {
    //     rain_heat = 0.0;
    // } else {
    //     rain_heat = calc_glacier_rain_heat(srf_tmp, air_tmp, pcp_rat, dt);
    // }
    // 地热通量
    ground_heat = 0.0;
    
    // 辐射余项为以上变量之和
    netQ = shortwave_net + longwave_net + sensible_heat + latent_heat + rain_heat + ground_heat;

    // log_info("Srf_Tmp=%f Air_Tmp=%f abd=%f pcp=%f SR=%f LR=%f Air_Den=%f Air_Pre=%f Vap_Pre=%f dt=%f Ra=%f", 
    // srf_tmp, 
    // air_tmp, 
    // glc_abd,  
    // pcp_rat, 
    // shortwave_in, 
    // longwave_in,
    // air_den, 
    // air_pre, 
    // vap_pre, 
    // dt,
    // Ra);

    // log_info("netQ=%f SWnet=%f LWnet=%f SH=%f LH=%f LR=%f", netQ, shortwave_net, longwave_net, sensible_heat, latent_heat, rain_heat);

    return netQ;
}

void update_acc_glacier_melt(all_vars_struct *all_vars, dmy_struct *dmy, size_t Nbands) {

    // 计数变量
    unsigned short band;

    // 在每年的10月1日清空冰川累积融水
    if ((dmy->month==10)&&(dmy->day==1)&&(dmy->dayseconds==6*3600)) {

        // 冰川数据
        glacier_data_struct *glacier;
        // 积雪数据
        snow_data_struct *snow;

        for (band = 0; band < Nbands; band++) {
            glacier = &(all_vars->glacier[band]);
            // 清空冰川累积融水变量
            glacier->acc_melt = 0.0;
        }
    }
}


// /**
//  * @brief 
//  * 
//  * @param TSurf             冰川表面温度 CC
//  * @param NetRad            净辐射通量 为净长波与净短波之和 W/m^2
//  * @param SensibleHeat      显热通量 W/m^2
//  * @param LatentHeat        潜热通量 W/m^2
//  * @param GroundFlux        地热通量 W/m^2
//  * @param AdvectedEnergy    降水热通量 W/m^2
//  * @param RefreezeEnergy    在融化情景中 为融化所耗能量 W/m^2
//  * @return double RestTerm  能量余项 仅在非融化情景中 作为温变所耗能量 W/m^2
//  */
// double calc_glacier_energy_balance(double TSurf, 
//                                     double NetRad, 
//                                     double SensibleHeat, 
//                                     double LatentHeat, 
//                                     double GroundFlux, 
//                                     double AdvectedEnergy, 
//                                     double *RefreezeEnergy) {
//     double RestTerm = 0.0;
    
//     *RefreezeEnergy = 0.0;
//     RestTerm = NetRad + *SensibleHeat + *LatentHeat + *GroundFlux + *AdvectedEnergy;
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
//      * 但是与融化过程无关
//      */

//     return RestTerm;
// }

// void update_annual_glacier(all_vars_struct *all_vars, soil_con_struct soil_con) {

//     extern option_struct options;

//     // 计数变量
//     int b;
//     // 分带数
//     size_t Nbands;
//     // 格网面积
//     double cell_area;
//     // 冰川LUCC面积占比
//     double cv;
//     // 各分带面积
//     double *band_area;
//     // 各分带冰川面积
//     double glacier_area;
//     double *band_glacier_area;

//     // 初始化变量
//     Nbands = options.SNOW_BAND;
//     band_area = calloc(Nbands, sizeof(*(band_area)));
//     check_alloc_status(band_area, "Memory allocation error.");
//     band_glacier_area = calloc(Nbands, sizeof(*(band_glacier_area)));
//     check_alloc_status(band_glacier_area, "Memory allocation error.");

//     // 分带面积
//     for (b = 0; b < options.SNOW_BAND; b++) {
//         band_area[b] = cell_area * soil_con.AreaFract[b];
//     }

//     // 冰川面积
//     glacier_area = cv * cell_area;

//     // 冰川分带面积
//     for (b = 0; b < options.SNOW_BAND; b++) {
//         glacier_data_struct glacier = all_vars->glacier[b];
//         band_glacier_area[b] = glacier.coverage * cv * cell_area;
//     }

//     // 取得累计年融化垂直水通量 取得冰川表面积雪垂直水通量

//     // 计算各分带累计融化体积

//     // 计算分带各冰川新面积

//     // 计算各分带冰川新面积占总冰川面积的百分比 更新glacier

//     // 计算新的

//     // 
// }