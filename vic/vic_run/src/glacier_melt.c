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

int glacier_melt(   double            Le,                       // 蒸发潜热 无用
                    double            NetShortSnow,             // 输入 净短波
                    double            Tcanopy,                  // 冠层温度 无用 冰雪面不存在冠层
                    double            Tgrnd,                    // 地表温度 无用 不计算地热通量
                    double           *Z0,                       // 地表粗糙度
                    double            aero_resist,              // 输入 空气动力学阻力 
                    double           *aero_resist_used,         // 输入 稳定性校正后的空气动力学阻力
                    double            air_temp,                 // 输入 气温
                    double            coverage,                 // 冰雪面面积占比 在冰川高程分带中恒为1
                    double            delta_t,                  // Ts对应的秒
                    double            density,                  // 输入 空气密度
                    double            grnd_flux,                // 输入 地热通量
                    double            LongSnowIn,               // 输入 入射长波
                    double            pressure,                 // 输入 气压
                    double            rainfall,                 // 输入 降雨
                    double            snowfall,                 // 输入 降雪
                    double            vp,                       // 输入 水汽压
                    double            vpd,                      // 输入 水汽压差
                    double            wind,                     // 输入 风速
                    double            z2,                       // 
                    double           *NetLongSnow,              // 输入 净长波
                    double           *OldTSurf,                 // 上个Ts的表面温度
                    double           *melt,                     // 输出 总融化量
                    double           *save_Qnet,                // 输出 能量余项
                    double           *save_advected_sensible,   // 输出 临近无雪地区显热
                    double           *save_advection,           // 输出 降雨热
                    double           *save_deltaCC,             // 输出 ColdContent变化
                    double           *save_grnd_flux,           // 输出 地热通量
                    double           *save_latent,              // 输出 蒸发潜热
                    double           *save_latent_sub,          // 输出 升华潜热
                    double           *save_refreeze_energy,     // 输出 重冻结
                    double           *save_sensible,            // 输出 显热
                    int               UNSTABLE_SNOW,            // 是否为稳定积雪
                    int               iveg,                     // LUCC类型编号
                    int               band,                     // 高程分带编号
                    snow_data_struct *snow,                     // 积雪 只用于存储与输出 保证程序不崩溃 不作实际演算
                    glacier_data_struct *glacier) {             // 冰川 冰雪一体实际变量存储位置
    
    extern option_struct     options;
    extern parameters_struct param;

    double                   error;                             // 
    double                   DeltaPackCC;                       /* Change in cold content of the pack                           */
    double                   DeltaPackSwq;                      /* Change in snow water equivalent of the pack (m)              */
    double                   Ice;                               /* Ice content of snow pack (m)                                 */
    double                   InitialSwq;                        /* Initial snow water equivalent (m)                            */
    double                   MassBalanceError;                  /* Mass balance error (m)                                       */
    double                   MaxLiquidWater;                    /* Maximum liquid water content of pack (m)                     */
    double                   PackCC;                            /* Cold content of snow pack (J)                                */
    double                   PackSwq;                           /* Snow pack snow water equivalent (m)                          */
    double                   Qnet;                              /* Net energy exchange at the surface (W/m2)                    */
    double                   RefreezeEnergy;                    /* refreeze/melt energy in surface layer (W/m2)                 */
    double                   PackRefreezeEnergy;                /* refreeze/melt energy in pack layer (W/m2)                    */
    double                   RefrozenWater;                     /* Amount of refrozen water (m)                                 */
    double                   SnowFallCC;                        /* Cold content of new snowfall (J)                             */
    double                   AllMelt;                           /* Amount of both snow melt & glacier melt (m water equivalent) */
    double                   SnowMelt;                          /* Amount of snow melt (m water equivalent)                     */
    double                   GlacierMelt;                       /* Amount of glacier melt (m water equivalent)                  */
    double                   SurfaceCC;                         /* Cold content of snow pack (J)                                */
    double                   SurfaceSwq;                        /* Surface layer snow water equivalent (m)                      */
    double                   SnowFall;
    double                   RainFall;
    double                   advection;
    double                   deltaCC;
    double                   latent_heat;
    double                   latent_heat_sub;
    double                   sensible_heat;
    double                   advected_sensible_heat;
    double                   melt_energy = 0.;

    SnowFall = snowfall / MM_PER_M;                             /* convet to m                                                  */
    RainFall = rainfall / MM_PER_M;                             /* convet to m                                                  */

    /***
     * @brief
     * 此处需要解决的问题
     * 1. 降雨和降雪进入冰川表面后各物质分量分配
     * 2. 表面ColdContent计算
     * 3. 能量平衡变量的赋值
     * 4. 废弃变量的赋值
     */

    // 分离积雪中的 固体 和 液体
    InitialSwq = glacier->swq;
    (*OldTSurf) = glacier->surf_temp;
    Ice = glacier->swq - glacier->pack_water - glacier->surf_water;
    // 确定积雪是否可以被分为两层
    if (Ice > param.SNOW_MAX_SURFACE_SWE) {
        SurfaceSwq = param.SNOW_MAX_SURFACE_SWE;
        PackSwq = Ice - SurfaceSwq;
    } else {
        SurfaceSwq = Ice;
        PackSwq = 0.;
    }

    // 分别计算表层和底层的CC CC均为非正
    SurfaceCC = CONST_VCPICE_WQ * SurfaceSwq * glacier->surf_temp;
    PackCC = CONST_VCPICE_WQ * PackSwq * glacier->pack_temp;
    // 计算降雪输入的CC
    if (air_temp > 0.0) {
        SnowFallCC = 0.0;
    } else {
        SnowFallCC = CONST_VCPICE_WQ * SnowFall * air_temp;
    }
    // 降雪进入积雪 在两层积雪层中分配CC
    if (SnowFall > (param.SNOW_MAX_SURFACE_SWE - SurfaceSwq) && (param.SNOW_MAX_SURFACE_SWE - SurfaceSwq) > DBL_EPSILON) {
        DeltaPackSwq = SurfaceSwq + SnowFall - param.SNOW_MAX_SURFACE_SWE;
        if (DeltaPackSwq > SurfaceSwq) {
            DeltaPackCC = SurfaceCC + (SnowFall - param.SNOW_MAX_SURFACE_SWE) / SnowFall * SnowFallCC;
        } else {
            DeltaPackCC = DeltaPackSwq / SurfaceSwq * SurfaceCC;
        }
        SurfaceSwq = param.SNOW_MAX_SURFACE_SWE;
        SurfaceCC += SnowFallCC - DeltaPackCC;
        PackSwq += DeltaPackSwq;
        PackCC += DeltaPackCC;
    } else {
        SurfaceSwq += SnowFall;
        SurfaceCC += SnowFallCC;
    }
    // 按照CC更新两层积雪的温度
    if (SurfaceSwq > 0.0) {
        glacier->surf_temp = SurfaceCC / (CONST_VCPICE_WQ * SurfaceSwq);
    } else {
        // 当无雪的时候
        // 无需更新表面温度 此时仅计算冰川
        // glacier->surf_temp = 0.0;
    }
    if (PackSwq > 0.0) {
        glacier->pack_temp = PackCC / (CONST_VCPICE_WQ * PackSwq);
    } else {
        glacier->pack_temp = 0.0;
    }
    // 分别·更新积雪的固体和液体部分
    Ice += SnowFall;
    glacier->surf_water += RainFall;

    // 当目标温度为0.0时 计算能量余项
    Qnet = calc_glacier_energy_balance((double) 0.0, delta_t, aero_resist, aero_resist_used, z2, Z0,
                                        density, vp, LongSnowIn, Le, pressure, RainFall, NetShortSnow, vpd,
                                        wind, (*OldTSurf), coverage, glacier->depth, glacier->density,
                                        glacier->surf_water, SurfaceSwq, Tcanopy, Tgrnd,
                                        &advection, &advected_sensible_heat, &deltaCC, &grnd_flux, &latent_heat,
                                        &latent_heat_sub, NetLongSnow, &RefreezeEnergy, &sensible_heat,
                                        &glacier->vapor_flux, &glacier->blowing_flux, &glacier->surface_flux);

    if (Qnet == 0.0) {
        // 此时表面温度恒为0. 但是需要确定是发生重冻结还是融化
        glacier->surf_temp = 0.0;
        if (RefreezeEnergy >= 0.0) { // RefreezeEnergy=(-RestTerm) 此时发生重冻结
            // 冻结水量
            RefrozenWater = RefreezeEnergy / (CONST_LATICE * CONST_RHOFW) * delta_t;
            // 冻结水量不能超过液态水量 如果重冻结放热仍然不能填补能量平衡缺口 则不会进入此分支
            if (RefrozenWater > glacier->surf_water) {
                RefrozenWater = glacier->surf_water;
                RefreezeEnergy = RefrozenWater * CONST_LATICE * CONST_RHOFW / (delta_t);
            }
            // 如果有积雪/双层积雪 那么重冻结到表层
            // 水当量从液体中转移到固体中
            SurfaceSwq += RefrozenWater;
            Ice += RefrozenWater;
            glacier->surf_water -= RefrozenWater;
            if (glacier->surf_water < 0.0) {
                glacier->surf_water = 0.0;
            }
            // 不发生任何融化
            AllMelt = 0.;
            SnowMelt = 0.;
            GlacierMelt = 0.;
        } else { // RefreezeEnergy=(-RestTerm) 此时发生融化
            AllMelt = fabs(RefreezeEnergy) / (CONST_LATICE * CONST_RHOFW) * delta_t;
            if (AllMelt < Ice) { // 只融雪
                SnowMelt = AllMelt;
                GlacierMelt = 0.;
                if (SnowMelt <= PackSwq) {
                    // 融雪较少 融雪量没有超过底层雪量
                    // 水当量从底层积雪的固体转移到液体中
                    glacier->surf_water += SnowMelt;
                    PackSwq -= SnowMelt;
                    Ice -= SnowMelt;
                } else {
                    // 融雪较多 融雪量超过了底层雪量
                    // 清空底层积雪 水当量全部转移到表层
                    glacier->surf_water += SnowMelt + glacier->pack_water;
                    glacier->pack_water = 0.0;
                    PackSwq = 0.0;
                    Ice -= SnowMelt;
                    SurfaceSwq = Ice;
                }
            } else {
                // 融化量大于积雪量 
                // 所有积雪都融化
                // 剩余的为冰川融化
                SnowMelt = Ice;
                GlacierMelt = AllMelt - SnowMelt;
                // 所有积雪固液水当量都转移到积雪表面 并清空存储
                glacier->surf_water += AllMelt + glacier->pack_water;
                SurfaceSwq = 0.0;
                PackSwq = 0.0;
                Ice = 0.0;
                glacier->surf_temp = 0.0;
                glacier->pack_temp = 0.0;
            }
        }
    } else {
        // 此时 不存在融化 仅涉及冰雪表面的温度变化
        AllMelt = 0.;
        SnowMelt = 0.;
        GlacierMelt = 0.;
        glacier->surf_temp = root_brent((double) (glacier->surf_temp - param.SNOW_DT),
                                        (double) (glacier->surf_temp + param.SNOW_DT),
                                        glacier_energy_balance,
                                        delta_t, aero_resist, aero_resist_used, z2, Z0,
                                        density, vp, LongSnowIn, Le, pressure,
                                        RainFall, NetShortSnow, vpd,
                                        wind, (*OldTSurf), coverage,
                                        glacier->depth, glacier->density,
                                        glacier->surf_water, SurfaceSwq,
                                        Tcanopy, Tgrnd,
                                       &advection, &advected_sensible_heat,
                                       &deltaCC, &grnd_flux, &latent_heat,
                                       &latent_heat_sub, NetLongSnow,
                                       &RefreezeEnergy, &sensible_heat,
                                       &glacier->vapor_flux, &glacier->blowing_flux,
                                       &glacier->surface_flux);
        if (glacier->surf_temp <= -998) {
            // 此时 Brent 法无解
            // 采取上一时段的表面温度 并进行标记
            glacier->surf_temp = *OldTSurf;
            glacier->surf_temp_fbflag = 1;
            glacier->surf_temp_fbcount++;
        } else if (glacier->surf_temp > 0.) {
            // 冰雪表面温度非负
            // 预期在升温过程中 才会进入此分支
            /***
             * @bug 暂时不考虑升温融化
             */
            glacier->surf_temp = 0.;
        }
        Qnet = calc_glacier_energy_balance(glacier->surf_temp, delta_t, aero_resist, aero_resist_used, z2, Z0,
                                           density, vp, LongSnowIn, Le, pressure, RainFall, NetShortSnow, vpd,
                                           wind, (*OldTSurf), coverage, glacier->depth, glacier->density,
                                           glacier->surf_water, SurfaceSwq, Tcanopy, Tgrnd,
                                          &advection, &advected_sensible_heat, &deltaCC, &grnd_flux, &latent_heat,
                                          &latent_heat_sub, NetLongSnow, &RefreezeEnergy, &sensible_heat,
                                          &glacier->vapor_flux, &glacier->blowing_flux, &glacier->surface_flux);
        // 冻结表层积雪的融水
        SurfaceSwq += glacier->surf_water;
        Ice += glacier->surf_water;
        glacier->surf_water = 0.0;
    }

    // 计算出流
    /***
     * @bug 更新表层积雪含水量
     */
    MaxLiquidWater = param.SNOW_LIQUID_WATER_CAPACITY * SurfaceSwq;
    if (glacier->surf_water > MaxLiquidWater) {
        // 此时所有的融水都在表层
        melt[0] = glacier->surf_water - MaxLiquidWater;
        // 如果液态水饱和 则出流
        glacier->surf_water = MaxLiquidWater;
    } else {
        // 如果液态水没有饱和 则没有出流
        melt[0] = 0.0;
    }
    /***
     * @bug 顶层出流 全部汇入底层
     */
    glacier->pack_water += melt[0]; 
    // 出流液态水进入 底层积雪所有液态水凝结 可以放出的最大热量
    PackRefreezeEnergy = glacier->pack_water * CONST_LATICE * CONST_RHOFW;

    if (PackCC < -PackRefreezeEnergy) {
        // 即使所有液态水冻结 底层积雪CC仍然为负 即低于0度
        // 所有液态水冻结
        PackSwq += glacier->pack_water;
        Ice += glacier->pack_water;
        glacier->pack_water = 0.0;
        // 更新底层积雪温度
        if (PackSwq > 0.0) {
            PackCC = PackSwq * CONST_VCPICE_WQ * glacier->pack_temp + PackRefreezeEnergy;
            glacier->pack_temp = PackCC / (CONST_VCPICE_WQ * PackSwq);
            if (glacier->pack_temp > 0.) {
                glacier->pack_temp = 0.;
            }
        } else {
            glacier->pack_temp = 0.0;
        }
    } else {
        // 没有全部冻结
        // 冰水混合物温度为 0.
        glacier->pack_temp = 0.0;
        DeltaPackSwq = -PackCC / (CONST_LATICE * CONST_RHOFW);
        glacier->pack_water -= DeltaPackSwq;
        PackSwq += DeltaPackSwq;
        Ice += DeltaPackSwq;
    }
    // 底层液态水出流
    MaxLiquidWater = param.SNOW_LIQUID_WATER_CAPACITY * PackSwq;
    if (glacier->pack_water > MaxLiquidWater) {
        melt[0] = glacier->pack_water - MaxLiquidWater;
        glacier->pack_water = MaxLiquidWater;
    } else {
        melt[0] = 0.0;
    }

    // 更新积雪物理特性
    Ice = PackSwq + SurfaceSwq;
    if (Ice > param.SNOW_MAX_SURFACE_SWE) {
        // 积雪较厚 可以分两层
        SurfaceCC = CONST_VCPICE_WQ * glacier->surf_temp * SurfaceSwq;
        PackCC = CONST_VCPICE_WQ * glacier->pack_temp * PackSwq;
        if (SurfaceSwq > param.SNOW_MAX_SURFACE_SWE) {
            PackCC += SurfaceCC * (SurfaceSwq - param.SNOW_MAX_SURFACE_SWE) / SurfaceSwq;
            SurfaceCC -= SurfaceCC * (SurfaceSwq - param.SNOW_MAX_SURFACE_SWE) / SurfaceSwq;
            PackSwq += SurfaceSwq - param.SNOW_MAX_SURFACE_SWE;
            SurfaceSwq -= SurfaceSwq - param.SNOW_MAX_SURFACE_SWE;
        } else if (SurfaceSwq < param.SNOW_MAX_SURFACE_SWE) {
            PackCC -= PackCC * (param.SNOW_MAX_SURFACE_SWE - SurfaceSwq) / PackSwq;
            SurfaceCC += PackCC * (param.SNOW_MAX_SURFACE_SWE - SurfaceSwq) / PackSwq;
            PackSwq -= param.SNOW_MAX_SURFACE_SWE - SurfaceSwq;
            SurfaceSwq += param.SNOW_MAX_SURFACE_SWE - SurfaceSwq;
        }
        // 仅当存在底层积雪时 更新底层温度
        glacier->pack_temp = PackCC / (CONST_VCPICE_WQ * PackSwq);
        // 仅当存在积雪时 更新表层温度
        glacier->surf_temp = SurfaceCC / (CONST_VCPICE_WQ * SurfaceSwq);
    } else {
        // 积雪较薄 只计算表层 底层置空
        PackSwq = 0.0;
        PackCC = 0.0;
        glacier->pack_temp = 0.0;
    }

    // 更新积雪固液总量
    glacier->swq = Ice + glacier->pack_water + glacier->surf_water;
    if (glacier->swq == 0.0) {
        /***
         * @attention
         * 无雪时
         * 表层温度连续演算
         * 底层温度置空
         */
        glacier->pack_temp = 0.0;
    }

    glacier->melt = melt[0]*MM_PER_M;
    glacier->snow_melt = SnowMelt*MM_PER_M;
    glacier->glacier_melt = GlacierMelt*MM_PER_M;

    // 物质平衡残差
    /***
     * @bug 没有冰川项
     */
    MassBalanceError = (InitialSwq - glacier->swq) + (RainFall + SnowFall) - melt[0] + glacier->vapor_flux;

    // 将流出量转化为mm
    melt[0] *= MM_PER_M;
    // 物质平衡残差
    glacier->mass_error = MassBalanceError;
    snow->mass_error = MassBalanceError;
    // CC
    glacier->coldcontent = SurfaceCC;
    snow->coldcontent = SurfaceCC;
    // 水汽通量 换向
    glacier->vapor_flux *= -1;
    snow->vapor_flux *= -1.;

    *save_advection = advection;                                // 降水热
    *save_deltaCC = deltaCC;                                    // 获得
    *save_grnd_flux = grnd_flux;                                // 地热
    *save_latent = latent_heat;                                 // 潜热
    *save_latent_sub = latent_heat_sub;                         // 升华潜热
    *save_sensible = sensible_heat;                             // 显热
    *save_advected_sensible = advected_sensible_heat;           // 无雪地区输入显热
    *save_refreeze_energy = RefreezeEnergy;                     // 重冻结
    *save_Qnet = Qnet;                                          // 能量平衡余项

    return (0);                                                 // 运算成功 输出成功标识符
}