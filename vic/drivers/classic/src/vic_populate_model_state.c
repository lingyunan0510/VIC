/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components).
 *
 * If a state file is provided to the model then its
 * contents are checked to see if it agrees with the current simulation set-up,
 * if so it is used to initialize the model state.  If no state file is
 * provided the model initializes all variables with defaults and the user
 * should expect to throw out the beginning of the simulation period as model
 * start-up.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize the model state (energy balance, water balance, and
 *           snow components).
 *****************************************************************************/
void
vic_populate_model_state(all_vars_struct *all_vars,
                         filep_struct     filep,
                         size_t           cellnum,
                         soil_con_struct *soil_con,
                         veg_con_struct  *veg_con,
                         lake_con_struct  lake_con,
                         dmy_struct      *dmy_current)
{
    extern option_struct options;

    size_t               Nveg;
    int                  tmp_lake_idx;

    cell_data_struct   **cell;
    energy_bal_struct  **energy;
    lake_var_struct     *lake;
    snow_data_struct   **snow;
    /**
     * @brief Assign Glacier Data Struct
     * Added in 2022-01-28
     * Checked in 2022-02-11
     */
    glacier_data_struct *glacier;
    veg_var_struct     **veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    lake = &all_vars->lake_var;
    snow = all_vars->snow;
    /**
     * @brief Get Glacier Varibales in All Variables
     * Modified in 2022-02-10
     * Checked in 2022-02-11
     */
    glacier = all_vars->glacier;
    veg_var = all_vars->veg_var;

    Nveg = veg_con[0].vegetat_type_num;

    // Initialize all data structures to 0
    initialize_soil(cell, Nveg);
    initialize_snow(snow, Nveg);
    /**
     * @brief Initialize Glacier Variables
     * Modified in 2022-02-10
     * Checked in 2022-02-10
     */
    initialize_glacier(glacier);
    initialize_veg(veg_var, Nveg);
    if (options.LAKES) {
        tmp_lake_idx = lake_con.lake_idx;
        if (tmp_lake_idx < 0) {
            tmp_lake_idx = 0;
        }
        initialize_lake(lake, lake_con, soil_con, &(cell[tmp_lake_idx][0]),
                        false);
    }
    initialize_energy(energy, Nveg);

    // Read initial state from a file if provided
    if (options.INIT_STATE) {
        read_initial_model_state(filep.init_state, all_vars, Nveg,
                                 options.SNOW_BAND, cellnum, soil_con,
                                 lake_con);
    }
    else {
        // else generate a default state
        generate_default_state(all_vars, soil_con, veg_con, dmy_current);
        if (options.LAKES) {
            generate_default_lake_state(lake, soil_con, lake_con);
        }
    }

    // compute those state variables that are derived from the others
    compute_derived_state_vars(all_vars, soil_con, veg_con);
    if (options.LAKES) {
        compute_derived_lake_dimensions(lake, lake_con);
    }
}
