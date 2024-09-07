/**
 * @file solve_glacier.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief Core Code of Glacier Energy Balance
 * @version 0.1
 * @date 2022-02-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <vic_run.h>

double solve_glacier(char               overstory,
                     double             BareAlbedo,
                     double             LongUnderOut,          // LW from understory
                     double             MIN_RAIN_TEMP,
                     double             MAX_SNOW_TEMP,
                     double             new_snow_albedo,
                     double             Tcanopy,                // canopy air temperature
                     double             Tgrnd,                  // soil surface temperature
                     double             air_temp,               // air temperature
                     double             prec,
                     double             snow_grnd_flux,
                     double            *AlbedoUnder,
                     double            *Le,
                     double            *LongUnderIn,            // surface incomgin LW
                     double            *NetLongSnow,            // net LW at snow surface
                     double            *NetShortGrnd,           // net SW reaching ground
                     double            *NetShortSnow,           // net SW at snow surface
                     double            *ShortUnderIn,           // surfave incoming SW
                     double            *Torg_snow,
                     double            *aero_resist,
                     double            *aero_resist_used,
                     double            *coverage,               // best guess snow coverage
                     double            *delta_coverage,         // cover fract change
                     double            *delta_snow_heat,        // change in pack heat
                     double            *displacement,
                     double            *gauge_correction,
                     double            *melt_energy,
                     double            *out_prec,
                     double            *out_rain,
                     double            *out_snow,
                     double            *ppt,
                     double            *rainfall,
                     double            *ref_height,
                     double            *roughness,
                     double            *snow_inflow,
                     double            *snowfall,
                     double            *surf_atten,
                     double            *wind,
                     double            *root,
                     int                INCLUDE_SNOW,
                     size_t             Nveg,
                     unsigned short     iveg,
                     unsigned short     band,
                     double             dt,
                     size_t             hidx,
                     int                veg_class,
                     int               *UnderStory,
                     double            *CanopLayerBnd,
                     double            *dryFrac,
                     dmy_struct        *dmy,
                     force_data_struct *force,
                     energy_bal_struct *energy,
                     layer_data_struct *layer,
                     snow_data_struct  *snow,
                     soil_con_struct   *soil_con,
                     veg_var_struct    *veg_var,
                     glacier_data_struct *glacier) {

    // 
    extern option_struct options;
    extern parameters_struct param;

    // 从Solve Snow中借用的代码块 定义中间变量
    int                         ErrorFlag;
    double                      ShortOverIn;
    double                      melt;               // 输出的冰川径流 
    double                      old_depth;
    double                      old_swq;
    double                      rainonly;
    double                      tmp_grnd_flux;
    double                      store_snowfall;
    int                         month;
    int                         day_in_year;
    double                      density;
    double                      longwave;
    double                      pressure;
    double                      shortwave;
    double                      vp;
    double                      vpd;

    // 能量平衡本地变量 用于取值
    double                      glacier_albedo;     // 冰川表面反照率
    double                      glacier_snow_albedo;// 冰川表面的积雪反照率

    double                      tsurf;              // 冰雪表面温度
    double                      old_tsurf;          // 同上 
    double                      tbrent;             // 同上
    double                      ra;                 // 表面粗糙度

    double                      ShortwaveIn;        // 入射短波
    double                      ShortwaveNet;       // 净短波
    double                      LongwaveIn;         // 入射长波
    double                      LongwaveOut;        // 出射长波
    double                      LongwaveNet;        // 净长波
    double                      SensibleHeat;       // 显热
    double                      LatentHeat;         // 潜热
    double                      PcpHeat;            // 降水热
    double                      GroundHeat;         // 地热

    // 冰雪一体演进冰川表面物质-能量平衡
    // 在冰川LUCC中 solve_glacier取代了solve_snow
    // 需要按照solve_snow的样式定义物质变量和能量变量

    // 时间
    month = dmy->month;
    day_in_year = dmy->day_in_year;
    // 驱动数据
    density = force->density[hidx];
    longwave = force->longwave[hidx];
    pressure = force->pressure[hidx];
    shortwave = force->shortwave[hidx];
    vp = force->vp[hidx];
    vpd = force->vpd[hidx];

    melt = 0.; // 融化/出流量
    *ppt = 0.; // 下渗量 恒为0.

    // 融化能量
    (*melt_energy) = 0.;

    // 积雪热能delta
    (*delta_snow_heat) = 0.;

    // 雨雪分离
    rainonly = calc_rainonly(air_temp, prec, MAX_SNOW_TEMP, MIN_RAIN_TEMP);
    *snowfall = gauge_correction[SNOW] * (prec - rainonly);
    *rainfall = gauge_correction[RAIN] * rainonly;
    if (*snowfall < 1e-5) {
        *snowfall = 0.;
    }
    (*out_prec) = *snowfall + *rainfall;
    (*out_rain) = *rainfall;
    (*out_snow) = *snowfall;
    store_snowfall = *snowfall;

    // 潜热
    (*Le) = calc_latent_heat_of_vaporization(air_temp);

    // 垂直层标号
    *UnderStory = 2;
    // 表面覆盖率
    (*surf_atten) = 1.;
    // 输入短波 指针赋值
    (*ShortUnderIn) = shortwave;
    // 输入长波 指针赋值
    (*LongUnderIn) = longwave;
    // 积雪过程标识 
    glacier->snow = true; // 恒为正 表明有雪/降雪 即涉及积雪过程
    // 冰雪面积 恒为1.
    *coverage = glacier->coverage;
    // 冠层净长波
    energy->NetLongOver = 0;
    // 冠层入射长波
    energy->LongOverIn = 0;
    // 地表净短波 
    (*NetShortGrnd) = 0.;
    // 输入水通量 
    (*snow_inflow) += *rainfall + *snowfall;
    // 旧雪水当量
    old_swq = glacier->swq;

    if (glacier->swq > 0 && store_snowfall == 0) {
        // 无降雪 但是有积雪
        // VIC积雪反照率演算
        glacier->last_snow++;
        glacier_snow_albedo = snow_albedo(*snowfall, new_snow_albedo, glacier->swq, glacier->albedo_snow, glacier->coldcontent, dt, glacier->last_snow, glacier->MELTING);
    } else {
        // 新雪
        glacier->last_snow = 0;
        glacier_snow_albedo = new_snow_albedo;
    }
    /***
     * @attention 冰川表面反照率演算
     */
    glacier_albedo = calc_glacier_albedo(glacier->albedo_min, glacier_snow_albedo, glacier->swq, options.d_star);
    (*AlbedoUnder) = glacier_albedo;
    // 表面净短波
    (*NetShortSnow) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

    if ((dmy->year==2003)&&(dmy->month==7)&&(dmy->day==15)) {
        log_info("\n-------Before\nSWE %f\nSurfW %f\nPackW %f", glacier->swq, glacier->surf_temp, glacier->pack_water);
        // log_info("\n----Band %d\nTime %d\nSnw_Mlt1 %f\nSnw_Mlt2 %f", band, dmy->dayseconds, glacier->snow_melt, snow->melt);
        // log_info("\n----Band %d\nTime %d\nNetSW %f\nInSW %f\nAlbedo %f", band, dmy->dayseconds, (*NetShortSnow), (*ShortUnderIn), (*AlbedoUnder));
        // log_info("\n----Before\nBand %d\nTime %d\nSWE1 %f\nSWE2 %f", band, dmy->dayseconds, glacier->swq, snow->swq);
        // log_info("\nSnowFall %f\nRainFall %f\nTemp %f", *snowfall, *rainfall, air_temp);
        // log_info("\nNetSW %f\nNetLW %f", *NetShortSnow, *NetLongSnow);
        // log_info("\nMelt %f\nSnowMelt %f\nGlacierMelt %f", melt, glacier->snow_melt, glacier->glacier_melt);
    }

    // 计算
    ErrorFlag = glacier_melt((*Le), (*NetShortSnow), Tcanopy, Tgrnd, roughness, aero_resist[*UnderStory], 
                            aero_resist_used, air_temp, *coverage, dt, density, snow_grnd_flux, *LongUnderIn, 
                            pressure, *rainfall, *snowfall, vp, vpd, wind[*UnderStory], ref_height[*UnderStory], 
                            NetLongSnow, Torg_snow, &melt, &energy->error, &energy->advected_sensible, 
                            &energy->advection, &energy->deltaCC, &tmp_grnd_flux, &energy->latent, &energy->latent_sub, 
                            &energy->refreeze_energy, &energy->sensible, INCLUDE_SNOW, iveg, band, snow, glacier);
    // 演算失败 仅出现在积雪过程中 表明薄雪无法被单独作为一层进行能量平衡演算
    // 冰川过程中大概不会出现 大概不会 嗯
    if (ErrorFlag == ERROR) {
        return (ERROR);
    }

    if ((dmy->year==2003)&&(dmy->month==7)&&(dmy->day==15)) {
        // log_info("\n----Band %d\nTime %d\nSnw_Mlt1 %f\nSnw_Mlt2 %f", band, dmy->dayseconds, glacier->snow_melt, snow->melt);
        // log_info("\n----Band %d\nTime %d\nNetSW %f\nInSW %f\nAlbedo %f", band, dmy->dayseconds, (*NetShortSnow), (*ShortUnderIn), (*AlbedoUnder));
        // log_info("\n----After\nBand %d\nTime %d\nSWE1 %f\nSWE2 %f", band, dmy->dayseconds, glacier->swq, snow->swq);
        log_info("\nBand %d\nTime %d\nSensiH %f\nLatenH %f\nLatenHS %f\nNetSw %f\nNetLw %f\nNetQ %f\nTsurf %f\nTair %f\nMlt %f\nSnw_Mlt %f\nGlc_Mlt %f\nSWE %f\nSurfW %f\nPackW %f", 
                band, dmy->dayseconds, energy->sensible, energy->latent, energy->latent_sub, (*NetShortSnow), (*NetLongSnow), energy->refreeze_energy, 
                glacier->surf_temp, air_temp, melt, glacier->snow_melt, glacier->glacier_melt, glacier->swq*1000, glacier->surf_water, glacier->pack_water);
        // log_info("\nSnowFall %f\nRainFall %f\nTemp %f", *snowfall, *rainfall, air_temp);
        // log_info("\nNetSW %f\nNetLW %f", *NetShortSnow, *NetLongSnow);
        // log_info("\nMelt %f\nSnowMelt %f\nGlacierMelt %f", melt, glacier->snow_melt, glacier->glacier_melt);
    }

    // 表面反照率
    energy->AlbedoUnder = *AlbedoUnder;

    /***
     * @attention 
     * 在有雪的情况下 更新积雪属性
     * 在无雪的情况下 积雪属性初始化
     * MELTING仅标识积雪是否融化 当不存在积雪的情况下 此标识失效
     * coverage恒为1.0 以屏蔽地表能量平衡
     */
    if (glacier->swq > 0.) {
        // 雪密度
        if (glacier->surf_temp <= 0) {// 演进积雪密度
            glacier->density = snow_density(snow, *snowfall, old_swq, air_temp, dt);
        } else if (glacier->last_snow == 0) { // 新雪密度
            glacier->density = new_snow_density(air_temp);
        }
        // 更新雪深
        old_depth = glacier->depth;
        glacier->depth = CONST_RHOFW * glacier->swq / glacier->density;
        // 融化标识
        if (glacier->coldcontent >= 0 && ((soil_con->lat >= 0 && (day_in_year > 60 && day_in_year < 273)) || (soil_con->lat < 0 && (day_in_year < 60 || day_in_year > 273)))) {
            glacier->MELTING = true;
        } else if (glacier->MELTING && *snowfall > param.SNOW_TRACESNOW) {
            glacier->MELTING = false;
        }
    } else {
        // 重置反照率 雪深 雪密度
        glacier->albedo_snow = new_snow_albedo;
        glacier->density = 0.;
        glacier->depth = 0.;
        glacier->swq = 0.;
        // 清空表层/底层状态
        glacier->surf_water = 0;
        glacier->pack_water = 0;
        glacier->surf_temp = 0;
        glacier->pack_temp = 0;
        // 
        glacier->store_coverage = 1.0;
        glacier->store_swq = 0.0;
        glacier->store_snow = true;
        glacier->MELTING = false;
    }

    /***
     * @attention 无论积雪是否存在
     * 冰川都会阻隔向下的能量余项
     * 能量总是用于融雪/融冰 不考虑能量的向下/向外传递
     * 不用进行delta_coverage能量平衡校正
     */
    (*delta_coverage) = 0.;

    /***
     * @brief 当无雪的时候 重置积雪参数
     */
    if (glacier->swq == 0) {


        glacier->density = 0.;
        glacier->depth = 0.;
        // 
        glacier->surf_water = 0;
        glacier->pack_water = 0;
        glacier->surf_temp = 0;
        glacier->pack_temp = 0;
        // 
        snow->snow_distrib_slope = 0;
        snow->store_snow = true;
        snow->MELTING = false;
    }

    // 为了防止程序不正常输出 将glacier对象的属性拷贝至snow对象中
    // State
    snow->albedo = glacier->albedo;
    snow->coldcontent = glacier->coldcontent;
    snow->coverage = glacier->coverage;
    snow->density = glacier->density;
    snow->last_snow = glacier->last_snow;
    snow->MELTING = glacier->MELTING;
    snow->snow = glacier->snow;
    snow->depth = glacier->depth;
    // 表层
    snow->surf_temp = glacier->surf_temp;
    snow->surf_water = glacier->surf_water;
    snow->surf_temp_fbcount = glacier->surf_temp_fbcount;
    snow->surf_temp_fbflag = glacier->surf_temp_fbflag;
    // 深层
    snow->pack_temp = glacier->pack_temp;
    snow->pack_water = glacier->pack_water;
    // 输出 将被同步到snow中
    snow->store_coverage = glacier->store_coverage;
    snow->store_swq = glacier->store_swq;
    snow->store_snow = glacier->store_snow;
    // 积雪物理特性
    snow->swq = glacier->swq;
    snow->density = glacier->density;
    snow->depth = glacier->depth;
    snow->max_snow_depth = glacier->max_snow_depth;
    // 通量
    snow->blowing_flux = glacier->blowing_flux;
    snow->mass_error = glacier->mass_error;
    snow->Qnet = glacier->Qnet;
    snow->surface_flux = glacier->surface_flux;
    snow->transport = glacier->transport;
    snow->vapor_flux = glacier->vapor_flux;
    /***
     * @attention 融雪和融冰分开
     */
    snow->melt = glacier->snow_melt;

    // 降雨和降雪都被冰川/积雪消耗了
    // 通量不再向下传递
    (*rainfall) = 0.;
    (*snowfall) = 0.;

    energy->melt_energy *= -1.;

    return (melt);
}