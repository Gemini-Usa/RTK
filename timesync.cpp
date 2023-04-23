//
// Created by admin on 2023/3/15.
// Time Synchronization for Base and Rover Station
#include "function.h"
/* ʱ��ͬ������
 * b_range: ��վ�۲�ֵ
 * b_status: ��վ�Ƿ���й۲�ֵ
 * r_range: ����վ�۲�ֵ
 * r_status: ����վ�Ƿ���й۲�ֵ
 *
 * ����: Sync::SYN-ʱ����ͬ���� Sync::ROV-����վʱ����ڻ�վ Sync::BAS-��վʱ���������վ Sync::UNK-�޷�֪Ϥ */
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
