/**
 * @file initialize_glacier.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2022-01-28
 * Checked in 2022-02-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <vic_driver_shared_all.h>

/**
 * @brief 
 * 
 * @param glacier 
 */
void initialize_glacier(glacier_data_struct *glacier) {
    extern option_struct options;
    size_t i;
    for (i = 0; i < options.SNOW_BAND; i++) {
        glacier[i].coverage = 0.0;
        glacier[i].albedo = 0.0;
        glacier[i].albedo_min = 0.0;
        glacier[i].albedo_base = 0.0;
        glacier[i].surf_tmp = 0.0;
        glacier[i].adjust_tmp  = 0.0;
        glacier[i].METTING = false;
        glacier[i].qm = 0.0;
        glacier[i].melt = 0.0;
        glacier[i].acc_melt = 0.0;
    }
}