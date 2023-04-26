//
// Created by admin on 2023/4/26.
//

#include <string>
#include <iostream>
#include <format>
#include "function.h"

using namespace std;

void OutRcvSol(const RcvRes& rcvres, const Range& range, std::ostream& os)
{
    std::string fmt;
    fmt += "{:10.3f} ";						// SOW
    fmt += "{:14.10f} {:14.10f} {:10.3f} ";	// posxyz
    fmt += "{:14.10f} {:14.10f} {:7.3f} ";	// posblh
    fmt += "{:7.3f} {:7.3f} ";				// clock bias
    fmt += "{:6.3f} {:6.3f} ";				// PDOP and sigma
    fmt += "{:7.3f} {:7.3f} {:7.3f} ";		// velocity
    fmt += "{:7.3f} {:3.3f} ";				// clock drift and velstd
    fmt += "{:3d} {:3d} {:3d}\n";			// satnum
    auto args = std::make_format_args(
            range.time.sow,
            rcvres.rcvpos.xyz.c1, rcvres.rcvpos.xyz.c2, rcvres.rcvpos.xyz.c3,
            rcvres.rcvposblh.blh.c1 * R2D, rcvres.rcvposblh.blh.c2 * R2D, rcvres.rcvposblh.blh.c3,
            rcvres.dtr[0], rcvres.dtr[1],
            rcvres.PDOP, rcvres.sigma,
            rcvres.vx, rcvres.vy, rcvres.vz,
            rcvres.vdtr, rcvres.velsigma,
            rcvres.satnum[0], rcvres.satnum[1], range.sat_num
    );
    if (rcvres.Status) os << std::vformat("success " + fmt, args);
    else os << std::format("failed in {:10.3f}\n", range.time.sow);
}

void OutSatSol(const RcvRes& rcvres, const SatRes& satres, const Range& range, std::ostream& os)
{
    std::string fmt;
    fmt += "{:10.3f} ";											// SOW
    fmt += satres.sys == GPS ? "G{:0>2d} " : "C{:0>2d} ";	// prn
    fmt += "{:14.3f} {:14.3f} {:14.3f} {:10.3f} ";				// posxyz and clock
    fmt += "{:10.3f} {:10.3f} {:10.3f} {:10.3f} ";				// velocity and drift
    fmt += "{:7.3f} {:10.3f} {:10.3f}\n";						// trop and azel
    auto args = std::make_format_args(
            range.time.sow,
            satres.prn,
            satres.satpos.xyz.c1, satres.satpos.xyz.c2, satres.satpos.xyz.c3, satres.dts,
            satres.vx, satres.vy, satres.vz, satres.dv,
            satres.trop, GetAzim(rcvres.rcvpos, satres.satpos, WGS84) * R2D, satres.elev * R2D
    );
    os << std::vformat(fmt, args);
}

void OutRTKSol(const GPSTime &b_time, const GPSTime &r_time, const RcvRes &rcv_res, const DDObs &dd_obs,
               const Config &config, ostream &ofs)
{
    std:string fmt;
    fmt += "{:10.3f},{:10.3f},";
    fmt += "{:14.3f},{:14.3f},{:14.3f},";
    fmt += "{:10.3f},{:10.3f},{:10.3f},";
    fmt += "{:6.3f},{:6.3f},{:6.3f},{:6.3f},";
    fmt += "{},{},{:7.3f}\n";
    double rover[3];
    if (config.sol_format == 0) {// llh
        for (int i = 0; i < 3; ++i) {
            rover[i] = rcv_res.rcvposblh.dblh[i];
        }
    } else if (config.sol_format == 1) {// xyz
        for (int i = 0; i < 3; ++i) {
            rover[i] = rcv_res.rcvpos.dxyz[i];
        }
    } else {// enu
        // get base xyz coord
        if (config.basepos_type == 0) {
            BLH blh;
            for (int i = 0; i < 3; ++i) {
                blh.dblh[i] = config.basepos[i];
            }
            for (int i = 0; i < 3; ++i) {
                rover[i] = XYZ2ENU(BLH2XYZ(blh, WGS84), rcv_res.rcvpos, WGS84).denu[i];
            }
        } else {
            XYZ xyz;
            for (int i = 0; i < 3; ++i) {
                xyz.dxyz[i] = config.basepos[i];
            }
            for (int i = 0; i < 3; ++i) {
                rover[i] = XYZ2ENU(xyz, rcv_res.rcvpos, WGS84).denu[i];
            }
        }
    }
    auto args = std::make_format_args(
            r_time.sow, b_time.sow,
            rover[0], rover[1], rover[2],
            dd_obs.dpos[0], dd_obs.dpos[1], dd_obs.dpos[2],
            rcv_res.sigma, rcv_res.PDOP, rcv_res.HDOP, rcv_res.VDOP,
            dd_obs.dd_sat_num[0], dd_obs.dd_sat_num[1], dd_obs.ratio
    );
    if (config.out_mode == 0) std::cout << std::vformat(fmt, args);
    else ofs << std::vformat(fmt, args);
}
