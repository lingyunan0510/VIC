/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes all filefilenames before they are called by
 * the model.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize all filenames before they are called by the
 *           model.
 *****************************************************************************/
void
initialize_filenames()
{
    extern filenames_struct filenames;

    size_t                  i;

    strcpy(filenames.init_state, "MISSING");
    strcpy(filenames.statefile, "MISSING");
    strcpy(filenames.constants, "MISSING");
    strcpy(filenames.soil, "MISSING");
    strcpy(filenames.veg, "MISSING");
    strcpy(filenames.veglib, "MISSING");
    strcpy(filenames.snowband, "MISSING");
    /**
     * @brief Initialize Glacier File Before Being Called
     * Added in 2022-02-10
     * Checked in 2022-02-10
     */
    strcpy(filenames.glacierband, "MISSING");
    log_info("Function initialize_filenames %s", filenames.glacierband);
    strcpy(filenames.lakeparam, "MISSING");
    strcpy(filenames.result_dir, "MISSING");
    strcpy(filenames.log_path, "MISSING");
    for (i = 0; i < 2; i++) {
        strcpy(filenames.f_path_pfx[i], "MISSING");
    }
}

/******************************************************************************
 * @brief    Initialize all file pointers
 *****************************************************************************/
void
initialize_fileps()
{
    extern filep_struct filep;

    size_t              i;

    filep.globalparam = NULL;
    filep.constants = NULL;
    filep.init_state = NULL;
    filep.lakeparam = NULL;
    filep.snowband = NULL;
    /**
     * @brief Initialize Glacier Band File
     * Added in 2022-02-01
     * Cheched in 2022-02-10
     */
    filep.glacierband = NULL;
    printf("Function initialize_fileps %s", filep.glacierband);
    filep.soilparam = NULL;
    filep.statefile = NULL;
    filep.veglib = NULL;
    filep.vegparam = NULL;
    filep.logfile = NULL;
    for (i = 0; i < 2; i++) {
        filep.forcing[i] = NULL;
    }
}
