#include "function.h"
//
// Created by admin on 2023/3/16.
// Single Differencing Positioning Model
void FormSDEpochObs(const Range& b_obs, const Range& r_obs, SDEpochObs& sd_obs)
{
    short k = 0;
    SDSatObs obs;
    sd_obs.time = b_obs.time;
    for (short i = 0; i < b_obs.sat_num; ++i) {
        for (short j = 0; j < r_obs.sat_num; ++j) {
            if (!b_obs.obs[i].isDoubleFreq() || !r_obs.obs[j].isDoubleFreq()) continue;
            if (b_obs.obs[i].sys != r_obs.obs[j].sys) continue;
            if (b_obs.obs[i].prn != r_obs.obs[j].prn) continue;
            obs.sys = b_obs.obs[i].sys;
            obs.prn = b_obs.obs[i].prn;
            obs.dP[0] = r_obs.obs[j].P[0] - b_obs.obs[i].P[0];
            obs.dP[1] = r_obs.obs[j].P[1] - b_obs.obs[i].P[1];
            obs.dL[0] = r_obs.obs[j].L[0] - b_obs.obs[i].L[0];
            obs.dL[1] = r_obs.obs[j].L[1] - b_obs.obs[i].L[1];
            obs.valid = 1;
            obs.n_bas = i;
            obs.n_rov = j;
            sd_obs.com_obs[k].epoch = 1;
            sd_obs.com_obs[k].MW_smooth = sd_obs.com_obs[k].MW = obs.formMW();
            sd_obs.com_obs[k].P_GF = obs.formPGF();
            sd_obs.com_obs[k].L_GF = obs.formLGF();
            sd_obs.com_obs[k].P_IF = obs.formPIF();
            sd_obs.sat_obs[k++] = obs;
        }
    }
    sd_obs.sat_num = k;
}

void DetectCycleSlip(SDEpochObs& prev_obs, SDEpochObs& curr_obs)
{
    double d_time = gpst2Sec(curr_obs.time) - gpst2Sec(prev_obs.time);
    if (d_time > 1.1) return;
    double dMW, dGF;
    for (int i = 0; i < curr_obs.sat_num; ++i) {
        auto obs_prev = std::find_if(prev_obs.sat_obs, prev_obs.sat_obs + prev_obs.sat_num, [&curr_obs,i](const SDSatObs& o){
            return o.prn == curr_obs.sat_obs[i].prn;
        });
        if (!obs_prev) continue;
        auto j = std::distance(prev_obs.sat_obs, obs_prev);
        dMW = curr_obs.com_obs[i].MW - prev_obs.com_obs[j].MW_smooth;
        dGF = curr_obs.com_obs[i].L_GF - prev_obs.com_obs[j].L_GF;
        if (fabs(dMW) > 1.5 || fabs(dGF) > 0.5) {
            curr_obs.sat_obs[i].valid = -1;
            curr_obs.com_obs[i].epoch = 0;
            continue;
        }
        auto& epc = curr_obs.com_obs[i].epoch = prev_obs.com_obs[j].epoch + 1;
        curr_obs.com_obs[i].MW_smooth =  (1.0 / epc) * curr_obs.com_obs[i].MW + ((epc - 1.0) / epc) * prev_obs.com_obs[j].MW_smooth;
    }
}
