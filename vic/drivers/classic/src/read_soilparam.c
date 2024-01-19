/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads soil parameters for each grid cell.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine reads soil parameters for each grid cell.
 *****************************************************************************/
void
read_soilparam(FILE            *soilparam,
               soil_con_struct *temp,
               bool            *RUN_MODEL,
               bool            *MODEL_DONE)
{
    void ttrim(char *string);
    extern option_struct     options;
    extern veg_lib_struct   *veg_lib;
    extern parameters_struct param;

    char                     line[MAXSTRING];
    char                     tmpline[MAXSTRING];
    const char               delimiters[] = " \t";
    char                    *token;
    size_t                   layer;
    int                      i, tempint, j;
    double                   off_gmt;
    double                   tempdbl;
    size_t                   length;
    int                      Nbands, band;
    int                      flag;
    double                   tmp_depth;
    double                   tmp_depth2, tmp_depth2_save;
    double                   b, b_save;
    double                   bubble, bub_save;
    double                   tmp_max_moist;
    double                   tmp_resid_moist;
    double                   zwt_prime, zwt_prime_eff;
    double                   tmp_moist;
    double                   w_avg;
    double                   Zsum, dp;
    double                   tmpdp, tmpadj, Bexp;
    size_t                   k;
    size_t                   Nnodes;

    /** Read plain ASCII soil parameter file **/
    if ((fscanf(soilparam, "%d", &flag)) != EOF) {
        if (flag) {
            *RUN_MODEL = true;
        }
        else {
            *RUN_MODEL = false;
        }

        if (fgets(line, MAXSTRING, soilparam) == NULL) {
            log_err("Unexpected EOF while reading soil file");
        }
    }
    else {
        *MODEL_DONE = true;
        *RUN_MODEL = false;
    }

    if (!(*MODEL_DONE) && (*RUN_MODEL)) {
        strcpy(tmpline, line);
        ttrim(tmpline);
        if ((token = strtok(tmpline, delimiters)) == NULL) {
            log_err("Can't find values for CELL NUMBER in soil file");
        }
        sscanf(token, "%d", &(temp->gridcel));
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL LATITUDE in soil file");
        }
        sscanf(token, "%lf", &(temp->lat));
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL LONGITUDE in soil file");
        }
        sscanf(token, "%lf", &(temp->lng));

        /* read infiltration parameter */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for INFILTRATION in soil file");
        }
        sscanf(token, "%lf", &(temp->b_infilt));
        if (temp->b_infilt <= 0) {
            log_err("b_infilt (%f) in soil file is <= 0; b_infilt must "
                    "be positive", temp->b_infilt);
        }

        /* read fraction of baseflow rate */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for FRACTION OF BASEFLOW RATE "
                    "in soil file");
        }
        sscanf(token, "%lf", &(temp->Ds));

        /* read maximum baseflow rate */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for MAXIMUM BASEFLOW RATE in "
                    "soil file");
        }
        sscanf(token, "%lf", &(temp->Dsmax));

        /* read fraction of bottom soil layer moisture */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for FRACTION OF BOTTOM SOIL LAYER "
                    "MOISTURE in soil file");
        }
        sscanf(token, "%lf", &(temp->Ws));

        /* read exponential */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for EXPONENTIAL in soil file");
        }
        sscanf(token, "%lf", &(temp->c));

        /* read expt for each layer */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for EXPT for layer %zu in "
                        "soil file", layer);
            }
            sscanf(token, "%lf", &(temp->expt)[layer]);
            if (temp->expt[layer] < 3.0) {
                log_err("Exponent in layer %zu is %f < 3.0; This must be "
                        "> 3.0", layer, temp->expt[layer]);
            }
        }

        /* read layer saturated hydraulic conductivity */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SATURATED HYDRAULIC "
                        "CONDUCTIVITY for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Ksat)[layer]);
        }

        /* read layer phi_s */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for PHI_S for layer %zu in "
                        "soil file", layer);
            }
            sscanf(token, "%lf", &(temp->phi_s)[layer]);
        }

        /* read layer initial moisture */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for INITIAL MOISTURE for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->init_moist)[layer]);
            if (temp->init_moist[layer] < 0.) {
                log_err("Initial moisture for layer %zu cannot be "
                        "negative (%f)", layer, temp->init_moist[layer]);
            }
        }

        /* read cell mean elevation */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL MEAN ELEVATION in soil "
                    "file");
        }
        sscanf(token, "%lf", &(temp->elevation));

        /* soil layer thicknesses */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for LAYER THICKNESS for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->depth)[layer]);
        }
        /* round soil layer thicknesses to nearest mm */
        for (layer = 0; layer < options.Nlayer; layer++) {
            temp->depth[layer] =
                round(temp->depth[layer] * MM_PER_M) / MM_PER_M;
        }

        /* read average soil temperature */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for AVERAGE SOIL TEMPERATURE in "
                    "soil file");
        }
        sscanf(token, "%lf", &(temp->avg_temp));
        if ((options.FULL_ENERGY || options.LAKES) &&
            (temp->avg_temp > 100. || temp->avg_temp < -50)) {
            log_err("Need valid average soil temperature in degrees C to "
                    "run.  Full Energy model, %f is not acceptable.",
                    temp->avg_temp);
        }

        /* read soil damping depth */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for SOIL DAMPING DEPTH in soil "
                    "file");
        }
        sscanf(token, "%lf", &(temp->dp));

        /* read layer bubbling pressure */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for BUBBLING PRESSURE for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->bubble)[layer]);
            if ((options.FULL_ENERGY ||
                 options.FROZEN_SOIL) && temp->bubble[layer] < 0) {
                log_err("Bubbling pressure in layer %zu is %f < 0; "
                        "This must be positive for FULL_ENERGY = true or "
                        "FROZEN_SOIL = true", layer, temp->bubble[layer]);
            }
        }

        /* read layer quartz content */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for QUARTZ CONTENT for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->quartz)[layer]);
            if (options.FULL_ENERGY &&
                (temp->quartz[layer] > 1. || temp->quartz[layer] < 0)) {
                log_err("Need valid quartz content as a fraction to run "
                        "Full Energy model, %f is not acceptable.",
                        temp->quartz[layer]);
            }
        }

        /* read layer bulk density */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for mineral BULK DENSITY "
                        "for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->bulk_dens_min)[layer]);
            if (temp->bulk_dens_min[layer] <= 0) {
                log_err("layer %zu mineral bulk density (%f) must "
                        "be > 0", layer, temp->bulk_dens_min[layer]);
            }
        }

        /* read layer soil density */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for mineral SOIL DENSITY "
                        "for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->soil_dens_min)[layer]);
            if (temp->soil_dens_min[layer] <= 0) {
                log_err("layer %zu mineral soil density (%f) must "
                        "be > 0", layer, temp->soil_dens_min[layer]);
            }
            if (temp->bulk_dens_min[layer] >= temp->soil_dens_min[layer]) {
                log_err("layer %zu mineral bulk density (%f) must "
                        "be less than mineral soil density (%f)", layer,
                        temp->bulk_dens_min[layer],
                        temp->soil_dens_min[layer]);
            }
        }

        if (options.ORGANIC_FRACT) {
            /* read layer organic content */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for ORGANIC CONTENT for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->organic)[layer]);
                if (temp->organic[layer] > 1. || temp->organic[layer] < 0) {
                    log_err("Need valid volumetric organic soil "
                            "fraction when options.ORGANIC_FRACT is set "
                            "to true.  %f is not acceptable.",
                            temp->organic[layer]);
                }
            }

            /* read layer bulk density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for organic BULK "
                            "DENSITY for layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->bulk_dens_org)[layer]);
                if (temp->bulk_dens_org[layer] <= 0 && temp->organic[layer] >
                    0) {
                    log_warn("layer %zu organic bulk density (%f) must "
                             "be > 0; setting to mineral bulk density "
                             "(%f)", layer, temp->bulk_dens_org[layer],
                             temp->bulk_dens_min[layer]);
                    temp->bulk_dens_org[layer] = temp->bulk_dens_min[layer];
                }
            }

            /* read layer soil density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for organic SOIL DENSITY for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->soil_dens_org)[layer]);
                if (temp->soil_dens_org[layer] <= 0 && temp->organic[layer] >
                    0) {
                    log_warn("layer %zu organic soil density (%f) must be "
                             "> 0; setting to mineral soil density (%f)",
                             layer, temp->soil_dens_org[layer],
                             temp->soil_dens_min[layer]);
                    temp->soil_dens_org[layer] = temp->soil_dens_min[layer];
                }
                if (temp->organic[layer] > 0 && temp->bulk_dens_org[layer] >=
                    temp->soil_dens_org[layer]) {
                    log_err("layer %zu organic bulk density (%f) "
                            "must be less than organic soil density (%f)",
                            layer, temp->bulk_dens_org[layer],
                            temp->soil_dens_org[layer]);
                }
            }
        }
        else {
            for (layer = 0; layer < options.Nlayer; layer++) {
                temp->organic[layer] = 0.0;
                temp->bulk_dens_org[layer] = MISSING;
                temp->soil_dens_org[layer] = MISSING;
            }
        }

        if (options.BULK_DENSITY_COMB) {
            /* read soil bulk density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for SOIL BULK DENSITY for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->bulk_density)[layer]);
            }
        }
        else {
            for (layer = 0; layer < options.Nlayer; layer++) {
                temp->bulk_density[layer] =
                    (1 -
                     temp->organic[layer]) * temp->bulk_dens_min[layer] +
                    temp->organic[layer] * temp->bulk_dens_org[layer];
            }
        }

        /* read cell gmt offset */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for GMT OFFSET in soil file");
        }
        sscanf(token, "%lf", &off_gmt);

        /* read layer critical point */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for CRITICAL POINT for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Wcr[layer]));
        }

        /* read layer wilting point */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for WILTING POINT for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Wpwp[layer]));
        }

        /* read soil roughness */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for SOIL ROUGHNESS in soil file");
        }
        sscanf(token, "%lf", &(temp->rough));

        /* Overwrite default bare soil aerodynamic resistance parameters
           with the values taken from the soil parameter file */
        for (j = 0; j < MONTHS_PER_YEAR; j++) {
            veg_lib[veg_lib[0].NVegLibTypes].roughness[j] = temp->rough;
            veg_lib[veg_lib[0].NVegLibTypes].displacement[j] = temp->rough *
                                                               0.667 /
                                                               0.123;
        }

        /* read snow roughness */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for SNOW ROUGHNESS in soil file");
        }
        sscanf(token, "%lf", &(temp->snow_rough));

        /* read cell annual precipitation */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for ANNUAL PRECIPITATION in soil file");
        }
        sscanf(token, "%lf", &(temp->annual_prec));

        /* read layer residual moisture content */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for RESIDUAL MOISTURE CONTENT for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->resid_moist[layer]));
        }

        /* read frozen soil active flag */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for FROZEN SOIL ACTIVE FLAG in "
                    "soil file");
        }
        sscanf(token, "%d", &tempint);
        temp->FS_ACTIVE = (char)tempint;

        /* read minimum snow depth for full coverage */
        if (options.SPATIAL_SNOW) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SPATIAL SNOW in soil file");
            }
            sscanf(token, "%lf", &tempdbl);
            temp->max_snow_distrib_slope = tempdbl;
        }
        else {
            temp->max_snow_distrib_slope = 0;
        }

        /* read slope of frozen soil distribution */
        if (options.SPATIAL_FROST) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SPATIAL FROST in soil file");
            }
            sscanf(token, "%lf", &tempdbl);
            temp->frost_slope = tempdbl;
        }
        else {
            temp->frost_slope = 0;
        }

        /* If specified, read cell average July air temperature in the final
           column of the soil parameter file */
        if (options.JULY_TAVG_SUPPLIED) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for average July Tair in "
                        "soil file");
            }
            sscanf(token, "%lf", &tempdbl);
            temp->avgJulyAirTemp = tempdbl;
        }

        /*******************************************
           End of soil parameters for this grid cell
        *******************************************/

        /*******************************************
           Compute Soil Layer Properties
        *******************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            temp->soil_density[layer] =
                (1 -
                 temp->organic[layer]) * temp->soil_dens_min[layer] +
                temp->organic[layer] * temp->soil_dens_org[layer];
            if (temp->resid_moist[layer] == MISSING) {
                temp->resid_moist[layer] = param.SOIL_RESID_MOIST;
            }
            temp->porosity[layer] = 1.0 - temp->bulk_density[layer] /
                                    temp->soil_density[layer];
            temp->max_moist[layer] = temp->depth[layer] *
                                     temp->porosity[layer] * MM_PER_M;
        }

        /**********************************************
           Validate Soil Layer Thicknesses
        **********************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            if (temp->depth[layer] < MINSOILDEPTH) {
                log_err("Model will not function with layer %zu "
                        "depth %f < %f m.", layer, temp->depth[layer],
                        MINSOILDEPTH);
            }
        }
        if (temp->depth[0] > temp->depth[1]) {
            log_err("Model will not function with layer %d depth"
                    "(%f m) > layer %d depth (%f m).", 0, temp->depth[0],
                    1, temp->depth[1]);
        }

        /**********************************************
           Compute Maximum Infiltration for Upper Layers
        **********************************************/
        if (options.Nlayer == 2) {
            temp->max_infil = (1.0 + temp->b_infilt) * temp->max_moist[0];
        }
        else {
            temp->max_infil =
                (1.0 +
                 temp->b_infilt) * (temp->max_moist[0] + temp->max_moist[1]);
        }

        /****************************************************************
           Compute Soil Layer Critical and Wilting Point Moisture Contents
        ****************************************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            temp->Wcr[layer] *= temp->max_moist[layer];
            temp->Wpwp[layer] *= temp->max_moist[layer];
            temp->resid_moist[layer] *= temp->depth[layer] * MM_PER_M;
            if (temp->Wpwp[layer] > temp->Wcr[layer]) {
                log_err("Calculated wilting point moisture (%f mm) is "
                        "greater than calculated critical point moisture "
                        "(%f mm) for layer %zu.\n\tIn the soil parameter "
                        "file, Wpwp_FRACT MUST be <= Wcr_FRACT.",
                        temp->Wpwp[layer], temp->Wcr[layer], layer);
            }
            if (temp->Wpwp[layer] < temp->resid_moist[layer]) {
                log_err("Calculated wilting point moisture (%f mm) is "
                        "less than calculated residual moisture (%f mm) "
                        "for layer %zu.\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be >= resid_moist / (1.0 - "
                        "bulk_density/soil_density).", temp->Wpwp[layer],
                        temp->resid_moist[layer], layer);
            }
        }

        /**********************************************
           Validate Spatial Snow/Frost Params
        **********************************************/
        if (options.SPATIAL_SNOW) {
            if (temp->max_snow_distrib_slope < 0.0) {
                log_err("max_snow_distrib_slope (%f) must be positive.",
                        temp->max_snow_distrib_slope);
            }
        }

        if (options.SPATIAL_FROST) {
            if (temp->frost_slope < 0.0) {
                log_err("frost_slope (%f) must be positive.",
                        temp->frost_slope);
            }
        }
        for (k = 0; k < options.Nfrost; k++) {
            if (options.Nfrost == 1) {
                temp->frost_fract[k] = 1.;
            }
            else if (options.Nfrost == 2) {
                temp->frost_fract[k] = 0.5;
            }
            else {
                temp->frost_fract[k] = 1. / (options.Nfrost - 1);
                if (k == 0 || k == options.Nfrost - 1) {
                    temp->frost_fract[k] /= 2.;
                }
            }
        }

        /*************************************************
           If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
           parameters d1, d2, d3, and d4 to ARNO baseflow
           parameters Ds, Dsmax, Ws, and c
        *************************************************/
        if (options.BASEFLOW == NIJSSEN2001) {
            layer = options.Nlayer - 1;
            temp->Dsmax = temp->Dsmax *
                          pow(
                (double) (1. / (temp->max_moist[layer] - temp->Ws)),
                -temp->c) +
                          temp->Ds * temp->max_moist[layer];
            temp->Ds = temp->Ds * temp->Ws / temp->Dsmax;
            temp->Ws = temp->Ws / temp->max_moist[layer];
        }

        // Soil thermal node thicknesses and positions
        Nnodes = options.Nnode;
        dp = temp->dp;
        if (options.QUICK_FLUX) {
            /* node thicknesses */
            temp->dz_node[0] = temp->depth[0];
            temp->dz_node[1] = temp->depth[0];
            temp->dz_node[2] = 2. * (dp - 1.5 * temp->depth[0]);

            /* node depths (positions) */
            temp->Zsum_node[0] = 0;
            temp->Zsum_node[1] = temp->depth[0];
            temp->Zsum_node[2] = dp;
        }
        else {
            if (!options.EXP_TRANS) { 
                /* Compute soil node thicknesses
                   Nodes set at surface, the depth of the first layer,
                   twice the depth of the first layer, and at the
                   damping depth.  Extra nodes are placed equal distance
                   between the damping depth and twice the depth of the
                   first layer. */

                temp->dz_node[0] = temp->depth[0];
                temp->dz_node[1] = temp->depth[0];
                temp->dz_node[2] = temp->depth[0];
                temp->Zsum_node[0] = 0;
                temp->Zsum_node[1] = temp->depth[0];
                Zsum = 2. * temp->depth[0];
                temp->Zsum_node[2] = Zsum;
                tmpdp = dp - temp->depth[0] * 2.5;
                tmpadj = 3.5;
                for (k = 3; k < Nnodes - 1; k++) {
                    temp->dz_node[k] = tmpdp / (((double) Nnodes - tmpadj));
                    Zsum += (temp->dz_node[k] + temp->dz_node[k - 1]) / 2.;
                    temp->Zsum_node[k] = Zsum;
                }
                temp->dz_node[Nnodes -
                              1] =
                    (dp - Zsum - temp->dz_node[Nnodes - 2] / 2.) *
                    2.;
                Zsum +=
                    (temp->dz_node[Nnodes - 2] + temp->dz_node[Nnodes - 1]) /
                    2.;
                temp->Zsum_node[Nnodes - 1] = Zsum;
                if ((int) (Zsum * MM_PER_M + 0.5) !=
                    (int) (dp * MM_PER_M + 0.5)) {
                    log_err("Sum of thermal node thicknesses (%f) "
                            "in initialize_model_state do not "
                            "equal dp (%f), check initialization "
                            "procedure", Zsum, dp);
                }
            }
            else {
                // exponential grid transformation, EXP_TRANS = TRUE
                // calculate exponential function parameter
                // to force Zsum=dp at bottom node
                Bexp = logf(dp + 1.) / (double) (Nnodes - 1);
                // validate Nnodes by requiring that there be at
                // least 3 nodes in the top 50cm
                if (Nnodes < 5 * logf(dp + 1.) + 1) {
                    log_err("The number of soil thermal nodes (%zu) "
                            "is too small for the supplied damping "
                            "depth (%f) with EXP_TRANS set to "
                            "TRUE, leading to fewer than 3 nodes "
                            "in the top 50 cm of the soil column.  "
                            "For EXP_TRANS=TRUE, Nnodes and dp "
                            "must follow the relationship:\n"
                            "5*ln(dp+1)<Nnodes-1\n"
                            "Either set Nnodes to at least %d in "
                            "the global param file or reduce "
                            "damping depth to %f in the soil "
                            "parameter file.  Or set EXP_TRANS to "
                            "FALSE in the global parameter file.",
                            Nnodes, dp, (int) (5 * logf(dp + 1.)) + 2,
                            exp(0.2 * (Nnodes - 1)) + 1);
                }
                for (k = 0; k <= Nnodes - 1; k++) {
                    temp->Zsum_node[k] = expf(Bexp * k) - 1.;
                }
                if (temp->Zsum_node[0] > temp->depth[0]) {
                    log_err("Depth of first thermal node (%f) in "
                            "initialize_model_state is greater "
                            "than depth of first soil layer (%f); "
                            "increase the number of nodes or "
                            "decrease the thermal damping depth "
                            "dp (%f)", temp->Zsum_node[0], temp->depth[0], dp);
                }

                // top node
                k = 0;
                temp->dz_node[k] = temp->Zsum_node[k + 1] - temp->Zsum_node[k];
                // middle nodes
                for (k = 1; k < Nnodes - 1; k++) {
                    temp->dz_node[k] =
                        (temp->Zsum_node[k + 1] - temp->Zsum_node[k]) / 2. +
                        (temp->Zsum_node[k] - temp->Zsum_node[k - 1]) / 2.;
                }
                // bottom node
                k = Nnodes - 1;
                temp->dz_node[k] = temp->Zsum_node[k] - temp->Zsum_node[k - 1];
            } // end if !EXP_TRANS
        }

        /*******************************************************************
           Calculate grid cell area.
         ******************************************************************/
        compute_cell_area(temp);

        /*************************************************
           Allocate and Initialize Snow Band Parameters
        *************************************************/
        Nbands = options.SNOW_BAND;
        temp->AreaFract = calloc(Nbands, sizeof(*(temp->AreaFract)));
        check_alloc_status(temp->AreaFract, "Memory allocation error.");
        temp->BandElev = calloc(Nbands, sizeof(*(temp->BandElev)));
        check_alloc_status(temp->BandElev, "Memory allocation error.");
        temp->Tfactor = calloc(Nbands, sizeof(*(temp->Tfactor)));
        check_alloc_status(temp->Tfactor, "Memory allocation error.");
        temp->Pfactor = calloc(Nbands, sizeof(*(temp->Pfactor)));
        check_alloc_status(temp->Pfactor, "Memory allocation error.");
        temp->AboveTreeLine = calloc(Nbands, sizeof(*(temp->AboveTreeLine)));
        check_alloc_status(temp->AboveTreeLine, "Memory allocation error.");

        /** Set default values for factors to use unmodified forcing data **/
        for (band = 0; band < Nbands; band++) {
            temp->AreaFract[band] = 0.;
            temp->BandElev[band] = temp->elevation;
            temp->Tfactor[band] = 0.;
            temp->Pfactor[band] = 1.;
        }
        temp->AreaFract[0] = 1.;

        /*************************************************
           Compute soil moistures for various values of water table depth
        *************************************************/
        soil_moisture_from_water_table(temp, options.Nlayer);

        /*************************************************
           Compute soil albedo in PAR range (400-700nm)
           following eqn 122 in Knorr 1997
        *************************************************/
        if (options.CARBON) {
            temp->AlbedoPar = 0.92 * param.ALBEDO_BARE_SOIL - 0.015;
            if (temp->AlbedoPar < param.PHOTO_ALBSOIPARMIN) {
                temp->AlbedoPar = param.PHOTO_ALBSOIPARMIN;
            }
        }

        /*************************************************
           Miscellaneous terms for MTCLIM disaggregation
        *************************************************/
        /* Central Longitude of Current Time Zone */
        temp->time_zone_lng = off_gmt * 360. / HOURS_PER_DAY;
        /* Assume flat grid cell for radiation calculations */
        temp->slope = 0;
        temp->aspect = 0;
        temp->whoriz = 0;
        temp->ehoriz = 0;
    } // end if(!(*MODEL_DONE) && (*RUN_MODEL))
}
