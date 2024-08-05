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
        // State
        glacier[i].albedo = 0.0;
        glacier[i].coldcontent = 0.0;
        glacier[i].swq = 0.;
        glacier[i].coverage = 1.0;
        glacier[i].band_coverage = 0.0;
        glacier[i].last_snow = 0;
        glacier[i].METING = false;
        glacier[i].snow = false;
        // 表层
        glacier[i].surf_temp = 0.;
        glacier[i].surf_water = 0.;
        glacier[i].surf_temp_fbcount = 0;
        glacier[i].surf_temp_fbflag = false;
        // 深层
        glacier[i].pack_temp = 0.;
        glacier[i].pack_water = 0.;
        // 冰川
        glacier[i].albedo_min = 0.;
        glacier[i].gwe = 0.;
        glacier[i].adjust_tmp = 0.;
        // 输出 将被同步到snow中
        glacier[i].store_coverage = 1.0;
        glacier[i].store_swq = 0.;
        glacier[i].stroe_snow = true;
        // 积雪物理特性
        glacier[i].albedo_snow = 0.85;
        glacier[i].density = 0.;
        glacier[i].depth = 0.;
        glacier[i].max_snow_depth = 0.;
        // Flux
        glacier[i].blowing_flux = 0.;
        glacier[i].mass_error = 0.;
        glacier[i].melt = 0.;
        glacier[i].snow_melt = 0.;
        glacier[i].glacier_melt = 0.;
        glacier[i].Qnet = 0.;
        glacier[i].surface_flux = 0.;
        glacier[i].transport = 0.;
        glacier[i].vapor_flux = 0.;
    }
}