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
    sin_solar_zenith = sin(solar_zenith);

    solar_azimuth = compute_solar_azi(lat, lng, time_zone_lng, day_in_year, second);

    cos_theta = (cos_beta*cos_solar_zenith + sin_beta*sin_solar_zenith*cos(solar_azimuth-phi))/cos_solar_zenith;

    return cos_theta;
}

/**
 * @brief 计算逐Step短波辐射地形-太阳校正因子
 *        考虑到6h的Step步长可能会经历日出/日落
 *        选择模拟步长开始 中段 结束三个时刻中最高太阳高度计算短波辐射地形-太阳校正因子
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
double compute_max_cos_theta(double         lat,
                             double         lng,
                             double         time_zone_lng,
                             unsigned short day_in_year,
                             unsigned       second, 
                             double         slope, 
                             double         aspect) {

    // 在Step不同时间的时刻秒
    unsigned s1;
    unsigned s2;
    unsigned s3;
    // 太阳高度角最高时的时刻秒
    unsigned s0;
    // 在Step不同时间的太阳天顶角的cos
    double cz1;
    double cz2;
    double cz3;
    // 最高的太阳高度 此时太阳天顶角的cos也为最大
    double cz0;
    // 太阳高度角
    double solar_elevation;
    // 太阳方位角
    double solar_azimuth;

    // 短波辐射地形-太阳校正因子
    double cos_theta;

    s1 = second;
    s2 = second + SEC_PER_HOUR * 3;
    s3 = second + SEC_PER_HOUR * 6;

    cz1 = compute_coszen(lat, lng, time_zone_lng, day_in_year, s1);
    cz2 = compute_coszen(lat, lng, time_zone_lng, day_in_year, s2);
    cz3 = compute_coszen(lat, lng, time_zone_lng, day_in_year, s3);

    if ((cz1 <= 0.0) && (cz2 <= 0.0) && (cz3 <= 0.0)) {
        // 太阳高度始终低于地平线
        // 取时间步长中段
        s0 = s2;
    } else {
        // 太阳高度至少有一个高于地平线
        // 在时间步长的开始 中段 结尾之间选择太阳高度最高的时刻
        if (cz1 >= cz2) {
            cz0 = cz1;
            s0 = s1;
        } else {
            cz0 = cz2;
            s0 = s2;
        }
        if (cz3 >= cz0) {
            s0 = s3;
        }
    }

    cos_theta = compute_cos_theta(lat, lng, time_zone_lng, day_in_year, s0, slope, aspect);
    if (cos_theta <= 0.0) {
        solar_azimuth = (compute_solar_azi(lat, lng, time_zone_lng, day_in_year, s0))/CONST_PI*180;
        solar_elevation = 90.0-((compute_solar_zen(lat, lng, time_zone_lng, day_in_year, s0))/CONST_PI*180);
        log_warn("Invalid Cos_Theta %f in %d, %fE %fN, with Aspect: %f, Slope: %f, Solar Elevation: %f and: Solar Azimuth %f", 
        cos_theta, day_in_year, lng, lat, aspect, slope, solar_elevation, solar_azimuth);
        cos_theta = 0.0;
    }
    return cos_theta;
}