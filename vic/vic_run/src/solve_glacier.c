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
    // Melt Water
    double                   glacier_melt;

    double bg;
    double xa =  1.66707295e-03;
    double xb = -4.35292948e-02;
    double xc =  4.07314960e-01;
    double xd = -1.59978653e+00;
    double xe =  2.19952775e+00;

    if ((veg_class != 17) || (glacier->coverage <= 0.00)) {
        // No Glacier LUCC
        // No Glacier Melt
        glacier_melt = 0.0;
    } else {
        // Tair
        double tair;
        tair = air_temp;
        double x;
        if (tair > 0) {
            if (dmy->month >= 8) {
                x = dmy->month - 7;
            } else {
                x = dmy->month + 5;
            }
            bg = options.DD*((xa*x*x*x*x)+(xb*x*x*x)+(xc*x*x)+(xd*x)+xe);
            glacier_melt = tair * bg * ( 1.0 - (snow->coverage));
        } else {
            glacier_melt = 0.0;
        }
        if (glacier_melt < 1e-3) {
            glacier_melt = 0.00;
        }
    }
    /**
     * @brief Success Log
     */
    ErrorFlag = 0;
    return glacier_melt;
}