#include <set>
#include <algorithm>
#include <initializer_list>
#include "function.h"
//
// Created by admin on 2023/3/21.
// Double Differencing Positioning Model
template<class T>
std::set<T> Intersect(const std::initializer_list<std::set<T>>& sets)
{
    std::set<T> out_set = *sets.begin();
    if (sets.size() < 2) return out_set;
    std::set<T> temp_set;
    for (const auto &s: sets) {
        temp_set = out_set;
        out_set.clear();
        std::set_intersection(temp_set.begin(), temp_set.end(),
                              s.begin(), s.end(),
                              std::inserter(out_set, out_set.begin()));
    }
    return out_set;
}

bool
DetermineRefSat(const Range &b_range, const Range &r_range, SatRes b_sat_res[], SatRes r_sat_res[],
                const SDEpochObs &sd_obs, DDObs &dd_obs, const Config &config)
{
    // necessary(1): cycle slip valid = 1
    // necessary(2): ephemeris OK
    // necessary(3): parity flag = 1
    // sort by CN0
    std::set<unsigned short> n1[2];
    std::set<unsigned short> n2[2];
    std::set<unsigned short> n3[2];
    std::set<unsigned short> n4[2];
    std::set<unsigned short> n5[2];
    std::set<unsigned short> n[2];
    float sorter[2]{0.0, 0.0};
    for (int i = 0; i < sd_obs.sat_num; ++i) {
        auto& obs = sd_obs.sat_obs[i];
        if (obs.valid == 1) n1[obs.sys].insert(obs.prn);
    }
    for (int i = 0; i < b_range.sat_num; ++i) {
        auto& obs = b_range.obs[i];
        if ((obs.parity[0] & obs.parity[1]) == 1 && obs.SNR[0] > config.snr_mask_b) n2[obs.sys].insert(obs.prn);
    }
    for (int i = 0; i < r_range.sat_num; ++i) {
        auto& obs = r_range.obs[i];
        if ((obs.parity[0] & obs.parity[1]) == 1 && obs.SNR[1] > config.snr_mask_r) n3[obs.sys].insert(obs.prn);
    }
    for (int i = 0; i < MAXSATNUM; ++i) {
        auto& sat = b_sat_res[i];
        if (sat.Status) n4[sat.sys].insert(sat.prn);
    }
    for (int i = 0; i < MAXSATNUM; ++i) {
        auto& sat = r_sat_res[i];
        if (sat.Status) n5[sat.sys].insert(sat.prn);
    }
    for (int i = 0; i < 2; ++i) {
        n[i] = Intersect({ n1[i], n2[i], n3[i], n4[i], n5[i] });
        if (n[i].empty()) dd_obs.dd_sat_num[i] = 0;
        else dd_obs.dd_sat_num[i] = static_cast<int>(n[i].size()) - 1;
    }

    if (n[0].empty() && n[1].empty()) return false;
    for (int i = 0; i < r_range.sat_num; ++i) {
        auto& obs = r_range.obs[i];
        if (n[obs.sys].count(obs.prn)) {
            if (obs.SNR[0] > sorter[obs.sys]) {
                dd_obs.ref_prn[obs.sys] = obs.prn;
                sorter[obs.sys] = obs.SNR[0];
            }
        }
    }
    int rcv = 0;
    for (auto& sat_res : {b_sat_res, r_sat_res}) {
        for (int i = 0; i < MAXSATNUM; ++i) {
            auto& prn = sat_res[i].prn;
            auto& sys = sat_res[i].sys;
            if (sys == NavSys::UNK) continue;
            if (!n[sys].contains(prn)) sat_res[i].Status = false;
            if (prn == dd_obs.ref_prn[sys]) rcv == 0 ? dd_obs.b_ref_idx[sys] = i : dd_obs.r_ref_idx[sys] = i;
        }
        ++rcv;
    }
    return true;
}
