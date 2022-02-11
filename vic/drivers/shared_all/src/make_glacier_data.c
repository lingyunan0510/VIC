/**
 * @file make_glacier_data.c
 * @author Yunan Ling (lingyunan@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2022-01-28
 * Modified in 2022-02-10
 * Checked in 2022-02-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <vic_driver_shared_all.h>

glacier_data_struct *
make_glacier_data() {
    extern option_struct options;

    glacier_data_struct *temp = NULL;

    temp = calloc(options.SNOW_BAND, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    return temp;
}