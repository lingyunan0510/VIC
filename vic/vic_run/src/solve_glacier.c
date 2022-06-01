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

    /**
     * @brief Snow Related Parameters
     */
    double                   snow_albedo;
    double                   snow_depth;

    /**
     * @brief 
     * 
     */
    double                   ra;

    /**
     * @brief Glacier Related Paramters
     */
    double                   tsurf;
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

    // log_info("%d", veg_class);

    if (veg_class != 17) {
        /**
         * @brief If Not In Glacier Area, No Melt Occur
         */
        glacier_melt = 0.0;
    } else {
        /**
         * @brief If In Glacier Area, Do The Math
         */
        day_in_year = dmy->day_in_year;

        //
        snow_albedo = snow->albedo;
        //
        snow_depth = snow->depth;

        // Tair ()
        tair = air_temp;

        // Air Density (kg/m^3)
        density = force->density[hidx];
        // Incoming Longwave Radiation (W/m^2)
        longwavein = force->longwave[hidx];
        // Air Pressure (kPa)
        pressure = force->pressure[hidx];
        // Incoming Shortwave Radiation (W/m^2)
        shortwavein = force->shortwave[hidx];
        // Vapor Pressure (kPa)
        vp = force->vp[hidx];

        /**
         * @brief May Be Modified Later
         * Marked By Yunan Ling in 2022-03-01
         */
        if (snow->depth > 0.0) {
            tsurf = snow->surf_temp;
        } else {
            if (Tgrnd <= 0.0) {
                tsurf = Tgrnd;
            } else {
                tsurf = 0.0;
            }
        }
        // tsurf = 0.0;

        if (wind[*UnderStory] > 0.0) {
            ra = aero_resist[*UnderStory] / StabilityCorrection(ref_height[*UnderStory], 0.f, tsurf, Tcanopy, wind[*UnderStory], roughness[2]);
        } else {
            ra = param.HUGE_RESIST;
        }

        // Calculate Glacier Albedo
        glacier_albedo = calc_glacier_albedo(tair, snow_albedo, snow_depth);
        // log_info("air temperature is %f", tair);
        // log_info("snow depth is %f", snow_depth);
        // log_info("glacier albedo is %f", glacier_albedo);
        // Calculate Glacier Outgoing Longwave
        longwaverout = calc_glacier_outgoing_longwave_radiation(tsurf, 1.0);
        // log_info("longwaverout is %f", longwaverout);
        // Calculate Glacier Net Shortwave
        NetRadiation = shortwavein * (1 - glacier_albedo) + longwavein - longwaverout;
        // log_info("netradiation is %f", NetRadiation);
        // Calculate Glacier Sensible Heat
        SensibleHeat = calc_glacier_sensible_heat(density, tair, tsurf, ra);
        // log_info("sensibleheat is %f", SensibleHeat);
        // Calculate Glacier Latent Heat
        // log_info("pressure is %f", pressure);
        // log_info("vapor pressure is %f", vp);
        LatentHeat = calc_glacier_latent_heat(density, pressure, tsurf, vp, ra);
        // log_info("latentheat is %f", LatentHeat);
        // Calculate Glacier Ground Heat
        GroundHeat = 0.0;
        // Calculate Glacier Heat 
        if (snow->depth > 0.0) {
            PcpHeat = 0.0;
        } else {
            PcpHeat = calc_glacier_rain_heat(tsurf, tair, (*snowfall+*rainfall), dt);
        }
        // Q net
        Qm = NetRadiation - SensibleHeat - LatentHeat + GroundHeat + PcpHeat;

        if ((Qm>0)&(tsurf==0.0)) {
            glacier_melt = Qm / (CONST_LATICE * CONST_RHOFW) * dt; 
        } else {
            // 
            glacier_melt = 0.0;
        }

        // log_info("%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f", 
        //     band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, dt,
        //     snow_depth, snow_albedo, glacier_albedo, ra, tsurf, 
        //     NetRadiation, SensibleHeat, LatentHeat, Qm, glacier_melt*MM_PER_M);
        
        // log_info("%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
        //     band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, 
        //     tair, density, pressure, vp, shortwavein, longwavein, longwaverout, 
        //     snow_depth, snow_albedo, glacier_albedo, ra, tsurf, 
        //     NetRadiation, SensibleHeat, LatentHeat, Qm, glacier_melt*MM_PER_M);

        // log_info("%d %d %d %d %d %10f %10f %10f", 
        // band, dmy->year, dmy->month, dmy->day, dmy->dayseconds, *out_prec, *out_rain, *out_snow);

        /**
         * @brief Modify Before Submit It 
         * Marked By Yunan Ling In 2022-03-05
         */
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