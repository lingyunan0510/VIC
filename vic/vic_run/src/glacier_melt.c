/**
 * @file glacier_melt.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2022-01-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <vic_run.h>

int glacier_melt(double Le,
                    double            NetShortSnow,     // net SW at absorbed by snow
                    double            Tcanopy,
                    double            Tgrnd,
                    double           *Z0,               // roughness
                    double            aero_resist,      // aerodynamic resistance
                    double           *aero_resist_used, // stability-corrected aerodynamic resistance
                    double            air_temp,         // air temperature
                    double            coverage,         // snowpack cover fraction
                    double            delta_t,          // time step in secs
                    double            density,          // atmospheric density
                    double            grnd_flux,        // ground heat flux
                    double            LongSnowIn,       // incoming longwave radiation
                    double            pressure,
                    double            rainfall,
                    double            snowfall,
                    double            vp,
                    double            vpd,
                    double            wind,
                    double            z2,
                    double           *NetLongSnow,
                    double           *OldTSurf,
                    double           *melt,
                    double           *save_Qnet,
                    double           *save_advected_sensible,
                    double           *save_advection,
                    double           *save_deltaCC,
                    double           *save_grnd_flux,
                    double           *save_latent,
                    double           *save_latent_sub,
                    double           *save_refreeze_energy,
                    double           *save_sensible,
                    int               UNSTABLE_SNOW,
                    int               iveg,
                    int               band,
                    glacier_data_struct *glacier) {

    double Qnet = 0;
    double Tsurf;
    double *glacier_melt;

    double SnowFall = snowfall / MM_PER_M; /* convet to m */
    double RainFall = rainfall / MM_PER_M; /* convet to m */

    // To Be Filled
    // Qnet = glacier_energy_balance();

    if (Qnet > 0) {
        // melt
        glacier->surf_tmp = 0.0;
        *glacier_melt = Qnet / (CONST_LATICE * CONST_RHOFW) * delta_t; 
    } else {
        // don't melt
        glacier->surf_tmp = 0.0;
        *glacier_melt = 0.0;
    }

    return(0);
}