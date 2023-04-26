#include <map>
#include "function.h"
//
// Created by admin on 2023/3/27.
// Realtime Kinematic Float and Fixed Solution
static double lam[2][2]{{CLIGHT/FREQ_L1,CLIGHT/FREQ_L2},{CLIGHT/FREQ_B1,CLIGHT/FREQ_B3}};

template<class T>
concept SatType = requires { T().prn && T().sys; };

template<SatType T>
const T* find_by_sat(const T* sat_arr, unsigned int size, unsigned long prn, NavSys sys)
{
    return std::find_if(sat_arr, sat_arr + size, [=](const T& sat_obj){ return sat_obj.prn == prn && sat_obj.sys == sys; });
}

void LSE(const double B[], const double P[], const double W[], int obs, int par, double* x, double* Q)
{
    auto BT = std::make_unique<double[]>(par*obs);
    MatrixTranspose(B, BT.get(), obs, par);
    auto BTP = std::make_unique<double[]>(par*obs);
    MatrixMultiply(BT.get(), P, BTP.get(), par, obs, obs);
    auto N = std::make_unique<double[]>(par*par);
    MatrixMultiply(BTP.get(), B, N.get(), par, obs, par);
    auto L = std::make_unique<double[]>(par);
    MatrixMultiply(BTP.get(), W, L.get(), par, obs, 1);
    MatrixInverse(N.get(), Q, par);
    MatrixMultiply(Q, L.get(), x, par, par, 1);
}

void AccuracyAssess(const double W[], const double P[], int obs, int par, RcvRes& res)
{
    auto WT = std::make_unique<double[]>(obs);
    MatrixTranspose(W, WT.get(), obs, 1);
    auto WTP = std::make_unique<double[]>(obs);
    MatrixMultiply(WT.get(), P, WTP.get(), 1, obs, obs);
    auto D = std::make_unique<double[]>(1);
    MatrixMultiply(WTP.get(), W, D.get(), 1, obs, 1);
    int t = obs - par;
    res.sigma = sqrt(D[0] / t);
}

void CalculateDOP(const double Q[], int par, RcvRes& res)
{
    res.PDOP = res.HDOP = res.VDOP = 0;
    for (int i = 0; i < 3; ++i) {
        res.PDOP += Q[i * par + i];
        res.HDOP += i < 2 ? Q[i * par + i] : 0;
        res.VDOP += i == 2 ? Q[i * par + i] : 0;
    }
    res.PDOP = sqrt(res.PDOP);
    res.HDOP = sqrt(res.HDOP);
    res.VDOP = sqrt(res.VDOP);
}

bool RTKFloat(RcvRes &r_res, const SatRes b_sat_res[], const SatRes r_sat_res[], const SDEpochObs &sd_obs,
              DDObs &dd_obs, double *a, double *Q_a)
{
    XYZ r_b0;
    r_b0.dxyz[0] = -2267804.5263;
    r_b0.dxyz[1] = 5009342.3723;
    r_b0.dxyz[2] = 3220991.8632;
    XYZ r_r0{r_res.rcvpos};
    int sat_num = dd_obs.dd_sat_num[0] + dd_obs.dd_sat_num[1];
    double dd_rho, los_rs[3], los_r1[3], sd_los[3];
    double p[2];
    p[0] = dd_obs.dd_sat_num[0]/(dd_obs.dd_sat_num[0]+1.0);
    p[1] = dd_obs.dd_sat_num[1]/(dd_obs.dd_sat_num[1]+1.0);
    double p_major[2][4]{{p[0]/2.0,p[0]/2.0,p[0]/2.0/0.001,p[0]/2.0/0.001},
                         {p[1]/2.0,p[1]/2.0,p[1]/2.0/0.001,p[1]/2.0/0.001}};
    p[0] = -1.0/(dd_obs.dd_sat_num[0]+1.0);
    p[1] = -1.0/(dd_obs.dd_sat_num[1]+1.0);
    double p_cor[2][4]{{p[0]/2.0,p[0]/2.0,p[0]/2.0/0.001,p[0]/2.0/0.001},
                       {p[1]/2.0,p[1]/2.0,p[1]/2.0/0.001,p[1]/2.0/0.001}};
    auto W = std::make_unique<double[]>(4*sat_num);
    auto B = std::make_unique<double[]>(4*sat_num*(3+2*sat_num));
    auto P = std::make_unique<double[]>(4*sat_num*4*sat_num);
    auto x = std::make_unique<double[]>(3+2*sat_num);
    auto Q = std::make_unique<double[]>((3+2*sat_num)*(3+2*sat_num));
    auto N = std::make_unique<double[]>(2*sat_num);// Ambiguity(Cycle)
    SatRes b_sat_1[2]{b_sat_res[dd_obs.b_ref_idx[0]],b_sat_res[dd_obs.b_ref_idx[1]] };
    SatRes r_sat_1[2]{r_sat_res[dd_obs.r_ref_idx[0]],r_sat_res[dd_obs.r_ref_idx[1]] };
    SDSatObs sd_1[2]{*find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, dd_obs.ref_prn[0], NavSys::GPS),
                     *find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, dd_obs.ref_prn[1], NavSys::BDS)};
    std::vector<std::pair<unsigned long, NavSys>> dd_sat_view;
    for (int i = 0; i < MAXSATNUM; ++i) {
        auto& sat = b_sat_res[i];
        if (!sat.Status) continue;
        if (sat.prn == dd_obs.ref_prn[sat.sys]) continue;
        dd_sat_view.emplace_back(sat.prn, sat.sys);
    }
    for (int i = 0; i < sat_num; ++i) {
        for (int j = 0; j < sat_num; ++j) {
            if (i == j) {
                for (int k = 0; k < 4; ++k) P[(4*i+k)*4*sat_num+(4*i+k)] = p_major[dd_sat_view[i].second][k];
            } else {
                if (dd_sat_view[i].second != dd_sat_view[j].second) continue;
                for (int k = 0; k < 4; ++k) P[(4*i+k)*4*sat_num+(4*j+k)] = p_cor[dd_sat_view[i].second][k];
            }
        }
    }
    int row{0}, ite{0};
    x[0] = 1.0;
    while (Norm(x.get(), 3) > 1.0E-6 && ite < 10) {
        for (auto& [prn, sys] : dd_sat_view) {
            auto b_sat = find_by_sat(b_sat_res, MAXSATNUM, prn, sys);
            auto r_sat = find_by_sat(r_sat_res, MAXSATNUM, prn, sys);
            auto sd = find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, prn, sys);
            dd_rho = GetDist(r_r0, r_sat->satpos) - GetDist(r_b0, b_sat->satpos)
                     - GetDist(r_r0, r_sat_1[sys].satpos) + GetDist(r_b0, b_sat_1[sys].satpos);
            GetLOSVector(r_sat->satpos, r_r0, los_rs);
            GetLOSVector(r_sat_1[sys].satpos, r_r0, los_r1);
            VectorSub(los_rs, los_r1, sd_los, 3);
            if (ite == 0) {
                N[2*row] = (sd->dL[0] - sd_1[sys].dL[0]) - (sd->dP[0] - sd_1[sys].dP[0]) / lam[sys][0];
                N[2*row+1] = (sd->dL[1] - sd_1[sys].dL[1]) - (sd->dP[1] - sd_1[sys].dP[1]) / lam[sys][1];
            }
            W[4*row] = sd->dP[0] - sd_1[sys].dP[0] - dd_rho;
            W[4*row+1] = sd->dP[1] - sd_1[sys].dP[1] - dd_rho;
            W[4*row+2] = (sd->dL[0] - sd_1[sys].dL[0] - N[2*row]) * lam[sys][0] - dd_rho;
            W[4*row+3] = (sd->dL[1] - sd_1[sys].dL[1] - N[2*row+1]) * lam[sys][1] - dd_rho;
            for (int j = 0; j < 4; ++j) {
                B[(4*row+j)*(3+2*sat_num)] = -sd_los[0];
                B[(4*row+j)*(3+2*sat_num)+1] = -sd_los[1];
                B[(4*row+j)*(3+2*sat_num)+2] = -sd_los[2];
            }
            B[(4*row+2)*(3+2*sat_num)+3+2*row] = lam[sys][0];
            B[(4*row+3)*(3+2*sat_num)+4+2*row] = lam[sys][1];
            ++row;
        }
        LSE(B.get(), P.get(), W.get(), 4*sat_num, 3+2*sat_num, x.get(), Q.get());
        for (int i = 0; i < 3; ++i) r_r0.dxyz[i] += x[i];
        for (int i = 0; i < 2 * sat_num; ++i) N[i] += x[3+i];
        ++ite;
        row = 0;
    }
    if (ite == 10) return false;
    r_res.rcvpos = r_r0;
    r_res.rcvposblh = XYZ2BLH(r_r0, WGS84);
    for (int i = 0; i < 2*sat_num; ++i) {
        a[i] = N[i];
        for (int j = 0; j < 2 * sat_num; ++j)
            Q_a[i * 2 * sat_num + j] = Q[(i + 3) * (3 + 2 * sat_num) + j + 3];
    }
    AccuracyAssess(W.get(), P.get(), 4*sat_num, 3 + 2 * sat_num, r_res);
    CalculateDOP(Q.get(), 3 + 2 * sat_num, r_res);
    VectorSub(r_r0.dxyz, r_b0.dxyz, dd_obs.dpos, 3);
    return true;
}

bool RTKFixed(RcvRes& r_res, const SatRes b_sat_res[], const SatRes r_sat_res[], const SDEpochObs &sd_obs,
              DDObs &dd_obs)
{
    XYZ r_b0;
    r_b0.dxyz[0] = -2267804.5263;
    r_b0.dxyz[1] = 5009342.3723;
    r_b0.dxyz[2] = 3220991.8632;
    XYZ r_r0{r_res.rcvpos};
    int sat_num = dd_obs.dd_sat_num[0] + dd_obs.dd_sat_num[1];
    double dd_rho, los_rs[3], los_r1[3], sd_los[3];
    double p[2];
    p[0] = dd_obs.dd_sat_num[0]/(dd_obs.dd_sat_num[0]+1.0);
    p[1] = dd_obs.dd_sat_num[1]/(dd_obs.dd_sat_num[1]+1.0);
    double p_major[2][2]{{p[0]/2.0/0.001,p[0]/2.0/0.001},
                         {p[1]/2.0/0.001,p[1]/2.0/0.001}};
    p[0] = -1.0/(dd_obs.dd_sat_num[0]+1.0);
    p[1] = -1.0/(dd_obs.dd_sat_num[1]+1.0);
    double p_cor[2][2]{{p[0]/2.0/0.001,p[0]/2.0/0.001},
                       {p[1]/2.0/0.001,p[1]/2.0/0.001}};
    auto W = std::make_unique<double[]>(2*sat_num);
    auto B = std::make_unique<double[]>(2*sat_num*3);
    auto P = std::make_unique<double[]>(2*sat_num*2*sat_num);
    auto x = std::make_unique<double[]>(3);
    auto Q = std::make_unique<double[]>(9);
    SatRes b_sat_1[2]{b_sat_res[dd_obs.b_ref_idx[0]],b_sat_res[dd_obs.b_ref_idx[1]] };
    SatRes r_sat_1[2]{r_sat_res[dd_obs.r_ref_idx[0]],r_sat_res[dd_obs.r_ref_idx[1]] };
    SDSatObs sd_1[2]{*find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, dd_obs.ref_prn[0], NavSys::GPS),
                     *find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, dd_obs.ref_prn[1], NavSys::BDS)};
    std::vector<std::pair<unsigned long, NavSys>> dd_sat_view;
    for (int i = 0; i < MAXSATNUM; ++i) {
        auto& sat = b_sat_res[i];
        if (!sat.Status) continue;
        if (sat.prn == dd_obs.ref_prn[sat.sys]) continue;
        dd_sat_view.emplace_back(sat.prn, sat.sys);
    }
    for (int i = 0; i < sat_num; ++i) {
        for (int j = 0; j < sat_num; ++j) {
            if (i == j) {
                for (int k = 0; k < 2; ++k) P[(2*i+k)*2*sat_num+(2*i+k)] = p_major[dd_sat_view[i].second][k];
            } else {
                if (dd_sat_view[i].second != dd_sat_view[j].second) continue;
                for (int k = 0; k < 2; ++k) P[(2*i+k)*2*sat_num+(2*j+k)] = p_cor[dd_sat_view[i].second][k];
            }
        }
    }
    int row{0}, ite{0};
    x[0] = 1.0;
    while (Norm(x.get(), 3) > 1.0E-6 && ite < 10) {
        for (auto& [prn, sys] : dd_sat_view) {
            auto b_sat = find_by_sat(b_sat_res, MAXSATNUM, prn, sys);
            auto r_sat = find_by_sat(r_sat_res, MAXSATNUM, prn, sys);
            auto sd = find_by_sat(sd_obs.sat_obs, sd_obs.sat_num, prn, sys);
            dd_rho = GetDist(r_r0, r_sat->satpos) - GetDist(r_b0, b_sat->satpos)
                     - GetDist(r_r0, r_sat_1[sys].satpos) + GetDist(r_b0, b_sat_1[sys].satpos);
            GetLOSVector(r_sat->satpos, r_r0, los_rs);
            GetLOSVector(r_sat_1[sys].satpos, r_r0, los_r1);
            VectorSub(los_rs, los_r1, sd_los, 3);
            W[2*row] = (sd->dL[0] - sd_1[sys].dL[0] - dd_obs.fixed_amb[2*row]) * lam[sys][0] - dd_rho;
            W[2*row+1] = (sd->dL[1] - sd_1[sys].dL[1] - dd_obs.fixed_amb[2*row+1]) * lam[sys][1] - dd_rho;
            for (int j = 0; j < 2; ++j) {
                B[(2*row+j)*3] = -sd_los[0];
                B[(2*row+j)*3+1] = -sd_los[1];
                B[(2*row+j)*3+2] = -sd_los[2];
            }
            ++row;
        }
        LSE(B.get(), P.get(), W.get(), 2*sat_num, 3, x.get(), Q.get());
        for (int i = 0; i < 3; ++i) r_r0.dxyz[i] += x[i];
        ++ite;
        row = 0;
    }
    if (ite == 10) return false;
    r_res.rcvpos = r_r0;
    r_res.rcvposblh = XYZ2BLH(r_r0, WGS84);
    AccuracyAssess(W.get(), P.get(), 2*sat_num, 3, r_res);
    CalculateDOP(Q.get(), 3, r_res);
    dd_obs.fixed = true;
    VectorSub(r_r0.dxyz, r_b0.dxyz, dd_obs.dpos, 3);
    return true;
}
