//
// Created by admin on 2023/3/15.
// Time Synchronization for Base and Rover Station
#include "function.h"
/* 时间同步函数
 * b_range: 基站观测值
 * b_status: 基站是否具有观测值
 * r_range: 流动站观测值
 * r_status: 流动站是否具有观测值
 *
 * 返回: Sync::SYN-时间是同步的 Sync::ROV-流动站时间快于基站 Sync::BAS-基站时间快于流动站 Sync::UNK-无法知悉 */
Sync timeSync(const Range& b_range, bool b_status, const Range& r_range, bool r_status)
{
    if (!b_status && !r_status) return Sync::UNK;
    if (!b_status && r_status) return Sync::ROV;
    if (b_status && !r_status) return Sync::BAS;
    auto b_time = gpst2Sec(b_range.time);
    auto r_time = gpst2Sec(r_range.time);
    if (b_time - r_time > maxtimediff) {
        return Sync::BAS;
    } else if (r_time - b_time > maxtimediff) {
        return Sync::ROV;
    } else {
        return Sync::SYN;
    }
}
