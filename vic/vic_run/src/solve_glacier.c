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
    // Month
    double                   month;
    // 
    double                   d_g;

    /**
     * @brief Metrological Forcing
     */
    double                   tair;

    /**
     * @brief Glacier Related Paramters
     */
    double                   glacier_melt;

    if (veg_class != 17) {
        /**
         * @brief If Not In Glacier Area, No Melt Occur
         */
        glacier_melt = 0.0;
    } else {
        /**
         * @brief If In Glacier Area, Do The Math
         */
        if (dmy->month >= 8) {
            month = (double) dmy->month-7;
        } else {
            month = (double) dmy->month+5;
        }
        // log_info("The Current Month is %f", month);
        d_g = options.b_g*((options.a*month*month*month*month)+(options.b*month*month*month)+(options.c*month*month)+(options.d*month)+(options.e));

        // Tair
        tair = air_temp;

        if (tair >= 0.0) {
            glacier_melt = (1-options.f_r)*d_g*(tair-0.0);
        } else {
            glacier_melt = 0.0;
        }
        
    }
    ErrorFlag = 0; 
    /**
     * @brief Transfer Unit
     */
    return glacier_melt;
}