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
                     snow_data_struct  *snow,
                     glacier_data_struct *glacier) {

    // 
    extern option_struct options;
    extern parameters_struct param;

    // Error Flag
    int                      ErrorFlag;
    // Day Of The Year
    int                      day_in_year;
    // Melt Water
    double                   melt;

    /**
     * @brief Metrological Forcing
     */
    double                   tair;
    double                   density;
    double                   longwavein;
    double                   pressure;
    double                   shortwavein;
    double                   vp;
    double                   rain;

    double                   cos_theta;

    /**
     * @brief Snow Related Parameters
     */
    double                   snow_albedo;
    double                   snow_depth;
    double                   snow_tsurf;

    /**
     * @brief 
     * 
     */
    double                   ra;

    /**
     * @brief Glacier Related Paramters
     */
    double                   tsurf;
    double                   old_tsurf;
    double                   tbrent;
    double                   glacier_melt;
    double                   glacier_albedo;

    /**
     * @brief Energy Balance Items
     */
    double                   longwaverout;
    double                   NetRadiation;
    double                   SensibleHeat;
    double                   LatentHeat;
    double                   GroundHeat;
    double                   PcpHeat;
    double                   Qm;
    double                   new_Qm;

    // log_info("%d", veg_class);
    // fprintf(LOG_DEST, "veg_class = %d\n", veg_class);

    if ((veg_class != 17) || (glacier->coverage <= 0.00)) {
        // No Glacier LUCC
        // No Glacier Melt
        glacier_melt = 0.0;
    } else {
        // Glacier LUCC

        // Date Struct
        day_in_year = dmy->day_in_year;
        // fprintf(LOG_DEST, "DOY = %d\n", day_in_year);
        // fprintf(LOG_DEST, "Band = %d\n", band);

        // Snow 
        snow_albedo = snow->albedo;
        snow_depth = snow->depth;
        snow_tsurf = snow->surf_temp;

        // Tair
        tair = air_temp + glacier->adjust_tmp; // 校正冰川表面温度
        // fprintf(LOG_DEST, "air_temp = %f\n", tair);
        // fprintf(LOG_DEST, "bnd_temp = %f\n", air_temp);
        // TSurf
        tsurf = glacier->surf_tmp;
        old_tsurf = glacier->surf_tmp;

        // // 短波辐射校正因子赋值
        cos_theta = force->cos_theta[band];

        // fprintf(LOG_DEST, "glacier_surf_temp = %f\n", tsurf);

        // Air Density (kg/m^3)
        density = force->density[hidx];
        // Incoming Longwave Radiation (W/m^2)
        longwavein = force->longwave[hidx];
        // Air Pressure (Pa)
        pressure = force->pressure[hidx];
        // Incoming Shortwave Radiation (W/m^2)
        /**
         * @brief 
         * 
         */
        shortwavein = force->shortwave[hidx];
        // shortwavein *= cos_theta;
        // if (shortwavein <= 0.0) {
        //     shortwavein = 0.0;
        // }
        // Vapor Pressure (Pa)
        vp = force->vp[hidx];
        // 
        rain = *rainfall;

        if (wind[*UnderStory] > 0.0) {
            ra = aero_resist[*UnderStory] / StabilityCorrection(ref_height[*UnderStory], 0.f, tsurf, Tcanopy, wind[*UnderStory], roughness[2]);
        } else {
            ra = param.HUGE_RESIST;
        }

        glacier_albedo = calc_glacier_albedo(glacier->albedo, snow_albedo, snow_depth, options.d_star);

        Qm = calc_glacier_energy_balance(tsurf, tair, glacier_albedo, rain, 
                                         shortwavein, longwavein, density, 
                                         pressure, vp, dt, ra);

        if ((tsurf >= 0.0) && (Qm >= 0.0)) {
            /**
             * 当表面温度为0 且 能量平衡余项为正时
             * 发生融化
             */
            glacier->METTING = true;
            glacier->surf_tmp = 0.0;
            glacier_melt = Qm / (CONST_LATICE * CONST_RHOFW) * dt;
            // fprintf(LOG_DEST, "glc_mlt = %f\n", glacier_melt);
        } else {
            /**
             * 否则 仅涉及温度变化
             * 采用Brent Root方法计算温度变化
             */
            glacier->METTING = false;
            tbrent = root_brent((double) (tsurf - param.SNOW_DT), 
                                (double) (tsurf + param.SNOW_DT), 
                                glacier_energy_balance, 
                                tair, glacier_albedo, rain, shortwavein, 
                                longwavein, density, pressure, vp, dt, ra);
            if (tbrent <= -998) { // 计算错误
                glacier->surf_tmp = old_tsurf;
                glacier->METTING = false;
                glacier_melt = 0.0;
            } else if (tbrent > 0.0) { // 冰川表面温度非正值
                glacier->surf_tmp = 0.0;
                glacier->METTING = false;
                glacier_melt = 0.0;
            } else {
                glacier->surf_tmp = tbrent;
                glacier->METTING = false;
                glacier_melt = 0.0;
            }
        }
        // if (band == 0) {
        //     fprintf(LOG_DEST, "after_gsf = %f\n", glacier->surf_tmp);
        // }

        // // /**
        // //  * @brief May Be Modified Later
        // //  * Marked By Yunan Ling in 2022-03-01
        // //  */
        // // if (snow->depth > 0.0) {
        // //     tsurf = snow->surf_temp;
        // // } else {
        // //     if (Tgrnd <= 0.0) {
        // //         tsurf = Tgrnd;
        // //     } else {
        // //         tsurf = 0.0;
        // //     }
        // // }
        // // // tsurf = 0.0;



        // // Calculate Glacier Albedo
        // glacier_albedo = calc_glacier_albedo(tair, snow_albedo, snow_depth);
        // // Calculate Glacier Outgoing Longwave
        // longwaverout = calc_glacier_outgoing_longwave_radiation(tsurf, 1.0);
        // // Calculate Glacier Net Shortwave
        // NetRadiation = shortwavein * (1 - glacier_albedo) + longwavein - longwaverout;
        // // Calculate Glacier Sensible Heat
        // SensibleHeat = calc_glacier_sensible_heat(density, tair, tsurf, ra);
        // // Calculate Glacier Latent Heat
        // LatentHeat = calc_glacier_latent_heat(density, pressure, tsurf, vp, ra);
        // // Calculate Glacier Ground Heat
        // GroundHeat = 0.0;
        // // Calculate Glacier Precipitation Heat 
        // PcpHeat = calc_glacier_rain_heat(tsurf, tair, (*snowfall+*rainfall), dt);

        // // Q net
        // Qm = NetRadiation - SensibleHeat - LatentHeat + GroundHeat + PcpHeat;

        // if ((Qm>0)&(tsurf==0.0)) {
            
        // } else {
        //     // 
        //     glacier_melt = 0.0;
        // }

        // log_info("%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f", 
        //     band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, dt,
        //     snow_depth, snow_albedo, glacier_albedo, ra, tsurf, 
        //     NetRadiation, SensibleHeat, LatentHeat, Qm, glacier_melt*MM_PER_M);
        
        // log_info("%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
        //     band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, 
        //     tair, density, pressure, vp, shortwavein, longwavein, longwaverout, 
        //     snow_depth, snow_albedo, glacier_albedo, ra, tsurf, 
        //     NetRadiation, SensibleHeat, LatentHeat, Qm, glacier_melt*MM_PER_M);

        // log_info("%d %d %d %d %d %f", 
        // band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, glacier_melt);

        /**
         * @brief 如果融水太少 就相当于没有融水
         * Marked By Yunan Ling In 2022-03-05
         */
        if (glacier_melt < 1e-5) {
            glacier_melt = 0.00;
        }
        // glacier_melt = 0.00;
    }
    /**
     * @brief Modify Before Submit It
     * Marked By Yunan Ling In 2022-03-05
     */
    ErrorFlag = 0;
    /**
     * @brief Transfer Unit
     */
    return glacier_melt;
}