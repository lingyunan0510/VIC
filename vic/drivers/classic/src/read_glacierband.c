/**
 * @file read_glacierband.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2022-01-30
 * Modified in 2023-10-12
 * Modified in 2022-02-10
 * Checked in 2022-02-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <vic_driver_classic.h>

void read_glacierband(FILE *glacierband, soil_con_struct *soil_con, glacier_con_struct *glacier_con) {
    
    extern option_struct     options;
    
    char ErrStr[MAXSTRING];
    size_t band;
    size_t Nbands;
    unsigned int cell;
    double total;
    double area_fract;
    double min_albedo;
    // 
    // double cell_slope;
    // double cell_aspect;
    double band_slope;
    double band_aspect;
    double band_elev;

    // double area;

    Nbands = options.SNOW_BAND;

    glacier_con->AreaFract = calloc(Nbands, sizeof(*(glacier_con->AreaFract)));
    check_alloc_status(glacier_con->AreaFract, "Memory allocation error.");

    glacier_con->BandElev = calloc(Nbands, sizeof(*(glacier_con->BandElev)));
    check_alloc_status(glacier_con->BandElev, "Memory allocation error.");

    glacier_con->MinAlbedo = calloc(Nbands, sizeof(*(glacier_con->MinAlbedo)));
    check_alloc_status(glacier_con->MinAlbedo, "Memory allocation error.");

    // glacier_con->Area = calloc(Nbands, sizeof(*(glacier_con->Area)));
    // check_alloc_status(glacier_con->Area, "Memory allocation error.");

    if (Nbands > 1) {
        // log_info("More than 1 Glacier Band");
        fscanf(glacierband, "%d", &cell);
        while (cell != soil_con->gridcel && !feof(glacierband)) {
            fgets(ErrStr, MAXSTRING, glacierband);
            fscanf(glacierband, "%d", &cell);
        }

        /* When arrives at the end of the file without anyfound */
        if (feof(glacierband)) {
            // /** 1 band is the default; no action necessary **/
            log_warn("Cannot find current gridcell (%i) in glacier band file; setting cell to have one elevation band.", soil_con->gridcel);
            return;
        }

        total = 0.;
        for (band = 0; band < Nbands; band++) {
            fscanf(glacierband, "%lf", &area_fract);
            if (area_fract < 0) {
                log_err("Negative glacier band area fraction (%f) read from file", area_fract);
            }
            // glacier_con->BandElev[band] = soil_con->BandElev[band];
            glacier_con->AreaFract[band] = area_fract;
            // log_info("%f", glacier_con->AreaFract[band]);
            total += area_fract;
            log_info("GBA Read From File %f", area_fract);
        }

        // log_info("Total Area of Glacier is %f", total);

        // if (total != 1.) {
        //     log_warn("Sum of the glacier band area fractions does not equal 1 (%f), dividing each fraction by the sum", total);
        //     for (band = 0; band < options.SNOW_BAND; band++) {
        //         glacier_con->AreaFract[band] /= total;
        //     }
        // } else {
        //     log_info("Glacier Band Area Correct");
        // }
        
        /**
         * @brief Read Minimum Albedo
         * 
         */
        for (band = 0; band < Nbands; band++) {
            fscanf(glacierband, "%lf", &min_albedo);
            glacier_con->MinAlbedo[band] = min_albedo;
        }
        // GSave Grid Scale Slope and Aspect
        // fscanf(glacierband, "%lf", &cell_slope);
        // soil_con->slope = cell_slope;
        // fscanf(glacierband, "%lf", &cell_aspect);
        // soil_con->aspect = cell_aspect;
        for (band = 0; band < Nbands; band++) {
            fscanf(glacierband, "%lf", &band_slope);
            soil_con->BandSlope[band] = band_slope;
        }
        for (band = 0; band < Nbands; band++) {
            fscanf(glacierband, "%lf", &band_aspect);
            soil_con->BandAspect[band] = band_aspect;
        }
        for (band = 0; band < Nbands; band++) {
            fscanf(glacierband, "%lf", &band_elev);
            glacier_con->BandElev[band] = (band_elev - soil_con->BandElev[band])*(-0.00498);
            if (band_elev == 0.) {
                glacier_con->BandElev[band] = 0.;
            }
        }
        // for (band = 0; band < Nbands; band++) {
        //     fscanf(glacierband, "%lf", &area);
        //     glacier_con->Area[band] = area;
        // }
    } else {
        log_info("Only 1 Glacier Band");
    }
}
