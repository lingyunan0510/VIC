/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the cosine of the solar zenith angle, given the
 * current location and date.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine computes the cosine of the solar zenith angle.
 *****************************************************************************/
double
compute_coszen(double         lat,
               double         lng,
               double         time_zone_lng,
               unsigned short day_in_year,
               unsigned       second)
{
    double coslat;
    double sinlat;
    double decl;
    double cosdecl;
    double sindecl;
    double cosegeom;
    double sinegeom;
    double coshss;
    double hour_offset;
    double cosh;
    double coszen;
    double hour;

    hour = second / SEC_PER_HOUR;

    /* calculate cos and sin of latitude */
    coslat = cos(lat * CONST_PI / 180);
    sinlat = sin(lat * CONST_PI / 180);

    /* calculate cos and sin of declination */
    decl = CONST_MINDECL * cos(((double) day_in_year + CONST_DAYSOFF) *
                               CONST_RADPERDAY);
    cosdecl = cos(decl);
    sindecl = sin(decl);

    /* calculate daylength as a function of lat and decl */
    cosegeom = coslat * cosdecl;
    sinegeom = sinlat * sindecl;
    coshss = -(sinegeom) / cosegeom;
    if (coshss < -1.0) {
        coshss = -1.0; /* 24-hr daylight */
    }
    if (coshss > 1.0) {
        coshss = 1.0; /* 0-hr daylight */
    }

    /* calculate cos of hour angle */
    hour_offset = (time_zone_lng - lng) * HOURS_PER_DAY / 360;
    cosh = cos((hour + hour_offset - 12) * CONST_PI / 12);

    /* calculate cosine of solar zenith angle */
    coszen = cosegeom * cosh + sinegeom;

    return coszen;
}

/**
 * @brief 计算太阳天顶角
 * 
 * @param lat           纬度 deg
 * @param lng           经度 deg        
 * @param time_zone_lng 时区中心经线 deg
 * @param day_in_year   儒略日 \
 * @param second        模型时刻秒 s
 * @return double       太阳天顶角 rad (0, pi)
 */
double compute_solar_zen(double         lat,
                         double         lng,
                         double         time_zone_lng,
                         unsigned short day_in_year,
                         unsigned       second) {
    // solar zenith (rad) ranges (0, pi)
    double solar_zen;
    // cos(solar_zen) ranges (-1, 1)
    double cos_solar_zen;

    cos_solar_zen = compute_coszen(lat, lng, time_zone_lng, day_in_year, second);
    solar_zen = acos(cos_solar_zen);
    return solar_zen;
}

/**
 * @brief 计算太阳方位角
 * 
 * @param lat           纬度 deg
 * @param lng           经度 deg
 * @param time_zone_lng 时区中心经线 deg
 * @param day_in_year   儒略日 \
 * @param second        模型时刻秒 s
 * @return double       太阳方位角 rad (0, 2pi)
 */
double compute_solar_azi(double         lat,
                         double         lng,
                         double         time_zone_lng,
                         unsigned short day_in_year,
                         unsigned       second) {
    
    // Solar Azimuth (rad) ranges (0, 2pi)
    double solar_azi;
    double cos_solar_azi;
    // Solar Zenith (rad) ranges (0, pi)
    double solar_zen;
    double cos_solar_zen;
    double sin_solar_zen;
    // Solar Declination (rad)
    double decl;
    double cos_decl;
    double sin_decl;
    // Local Time (H)
    double hour;
    double hour_offset;

    // Calculate Solar Zenith
    solar_zen = compute_solar_zen(lat, lng, time_zone_lng, day_in_year, second);
    cos_solar_zen = compute_coszen(lat, lng, time_zone_lng, day_in_year, second);
    sin_solar_zen = sin(solar_zen);
    
    // Calculate Solar Declination
    decl = CONST_MINDECL * cos(((double) day_in_year + CONST_DAYSOFF) * CONST_RADPERDAY);
    cos_decl = cos(decl);
    sin_decl = sin(decl);

    // Compute Loacl Time
    hour = second / SEC_PER_HOUR;
    hour_offset = (time_zone_lng - lng) * HOURS_PER_DAY / 360;
    hour = hour + hour_offset;

    // Calculate Cos Solar Azimuth
    cos_solar_azi = (sin_decl - (sin(lat*CONST_PI/180)*cos_solar_zen))/(cos(lat*CONST_PI/180)*sin_solar_zen);
    // log_warn("cos_sa %f", cos_solar_azi);
    // Get Solar Azimuth ranges (0, pi)
    solar_azi = acos(cos_solar_azi);
    // Adjust Solar Azimuth to (0, 2pi)
    if (hour > 12) {
        solar_azi = 2*CONST_PI - solar_azi;
    }

    return solar_azi;
}

/**
 * @brief 计算短波辐射地形-太阳校正因子
 * 
 * @param lat           纬度 deg
 * @param lng           经度 deg
 * @param time_zone_lng 时区中心经线 deg
 * @param day_in_year   儒略日 \
 * @param second        模型时刻秒 s 
 * @param slope         高程分带平均坡度 rad
 * @param aspect        高程分带平均坡向 rad
 * @return double       短波辐射地形-太阳校正因子 \
 */
double compute_cos_theta(double         lat,
                         double         lng,
                         double         time_zone_lng,
                         unsigned short day_in_year,
                         unsigned       second, 
                         double         slope, 
                         double         aspect) {
    // Output
    double cos_theta;
    // Mean Slope of Elevation Band (rad)
    double beta;
    double cos_beta;
    double sin_beta;
    // Mean Aspect of Elevation Band (rad)
    double phi;
    // Solar Zenith (rad)
    double solar_zenith;
    double cos_solar_zenith;
    double sin_solar_zenith;
    // Solar Azimuth (rad)
    double solar_azimuth;

    // Convert Slope to rad
    beta = slope;
    cos_beta = cos(beta);
    sin_beta = sin(beta);
    // Convert Aspect to rad
    phi = aspect;

    solar_zenith = compute_solar_zen(lat, lng, time_zone_lng, day_in_year, second);
    cos_solar_zenith = cos(solar_zenith);
    if (cos_solar_zenith <= 0) {
        return 0.0;
    }
    sin_solar_zenith = sin(solar_zenith);
    // log_warn("sz %f, cos_sz %f, sin_sz %f", solar_zenith, cos_solar_zenith, sin_solar_zenith);
    solar_azimuth = compute_solar_azi(lat, lng, time_zone_lng, day_in_year, second);
    // log_warn("sa %f", solar_azimuth);
    cos_theta = (cos_beta*cos_solar_zenith + sin_beta*sin_solar_zenith*cos(solar_azimuth-phi));

    if (cos_theta <= 0.0) {
        // 坡面法线角和太阳光线的夹角超过90°
        // 无直射 因子为0
        cos_theta = 0.0;
    }

    return cos_theta;
}

/**
 * @brief 计算时间段中间的短波辐射地形-太阳校正因子
 * 
 * @param lat           纬度 deg
 * @param lng           经度 deg
 * @param time_zone_lng 时区中心经线 deg
 * @param day_in_year   儒略日 \
 * @param second        模型时刻秒 s 
 * @param slope         高程分带平均坡度 rad
 * @param aspect        高程分带平均坡向 rad
 * @return double       短波辐射地形-太阳校正因子 \
 */
double compute_mean_cos_theta(double         lat,
                              double         lng,
                              double         time_zone_lng,
                              unsigned short day_in_year,
                              unsigned       second, 
                              double         slope, 
                              double         aspect) {
    unsigned int s0 = 0;
    unsigned int count = 0;
    double cz0 = 0;
    double sum_theta = 0.0;
    double mean_theta = 0.0;

    for (int i=0; i<12; i++) {
        s0 = second + (30*60)*i;
        cz0 = compute_coszen(lat, lng, time_zone_lng, day_in_year, s0);
        if (cz0 <= 0.0) {
            sum_theta += 0.0;
        } else {
            sum_theta += compute_cos_theta(lat, lng, time_zone_lng, day_in_year, s0, slope, aspect);
            count++;
        }
    }

    if (count==0) {
        mean_theta = 0.0;
    } else {
        mean_theta = (sum_theta/count);
    }

    return mean_theta;
}
