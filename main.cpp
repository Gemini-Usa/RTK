// main.cpp
// Single Point Positioning ()

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <format>
#include "constants.h"
#include "function.h"
#include "sockets.h"
#include "serial.h"
using namespace std;

void Test()
{
	int LenRead = 0;
	int LenRem = 0;
	bool Status = false;
	bool flag = false;
	unsigned char Buff[MAXNOVDLEN];
	Range range;
	Range temp_range;
	RcvRes res;
	RcvRes init;
	SatRes satres[MAXSATNUM];
	GPSEph geph[MAXGPSSAT];
	BDSEph beph[MAXBDSSAT];
	// BESTPOS={-2267804.5263,5009342.3723,3220991.8632}
	XYZ BESTPOS;
	BESTPOS.xyz.c1 = -2267804.5263;
	BESTPOS.xyz.c2 = 5009342.3723;
	BESTPOS.xyz.c3 = 3220991.8632;
	ENU rcvposenu;
	double ENUTransMat[9];
	GetENUTransMat(BESTPOS, WGS84, ENUTransMat);
	MatrixPrint(ENUTransMat, 3, 3);
	// 以二进制格式打开观测数据文件
	FILE* stream;
	FILE* ostream = fopen("400rcv.txt", "w");
	FILE* satstream = fopen("400sat.txt", "w");
	if (fopen_s(&stream, "oem719-202202080900-2.bin", "r+b") != 0)
	{
		printf("Problem opening the file\n");
		return;
	}
	fprintf(ostream, "FLAG\t SOW\t X\t Y\t Z\t B\t L\t H\t E\t N\t U\t GClk\t BClk\t PDOP\t SIGMA\t VX\t VY\t VZ\t VClk\t VSIGMA\t GNUM\t BNUM\t NUM\n");
	fprintf(satstream, "SOW\t prn\t X\t Y\t Z\t Clk\t VX\t VY\t VZ\t VClk\t TROP\t AZIM\t ELEV\n");
	

	// 循环读取数据进缓冲区
	LenRem = 0;
	while (feof(stream) == 0) 
	{
		LenRead = fread_s(Buff + LenRem, sizeof(unsigned char) * MAXNOVDLEN, sizeof(unsigned char), MAXNOVDLEN - LenRem, stream);
		if (LenRead < (MAXNOVDLEN - LenRem))
		{
			printf("finished!\n");
			fclose(stream);
			return;
		}
		// 调用解码模块，得到观测值和星历
		Status = DecodeNovOem7(Buff, LenRem, &range, geph, beph, 0);

		if (Status)
		{
			// 进行粗差探测
            DetectOutlier(&range, &temp_range);
			// SPP
			flag = SPP(&res, satres, &init, &range, geph, beph);
			if (flag)
			{
				if (res.Status == true)
				{
					SPV(&res, satres, &range);
					init = res;
					for (int i = 0; i < res.satnum[0]; i++)
					{
						fprintf(satstream, "%f G%02d %15.3f %15.3f %15.3f %10.3e %10.3f %10.3f %10.3f %10.3e %7.3f %10.3f %10.3f\n",
                                range.time.sow,
                                satres[i].prn, satres[i].satpos.xyz.c1, satres[i].satpos.xyz.c2, satres[i].satpos.xyz.c3, satres[i].dts,
                                satres[i].vx, satres[i].vy, satres[i].vz, satres[i].dv, satres[i].trop, GetAzim(res.rcvpos, satres[i].satpos, WGS84) * R2D, satres[i].elev * R2D);
					}
					for (int i = res.satnum[0]; i < res.satnum[0] + res.satnum[1]; i++)
					{
						fprintf(satstream, "%f C%02d %15.3f %15.3f %15.3f %10.3e %10.3f %10.3f %10.3f %10.3e %7.3f %10.3f %10.3f\n",
                                range.time.sow,
                                satres[i].prn, satres[i].satpos.xyz.c1, satres[i].satpos.xyz.c2, satres[i].satpos.xyz.c3, satres[i].dts,
                                satres[i].vx, satres[i].vy, satres[i].vz, satres[i].dv, satres[i].trop, GetAzim(res.rcvpos, satres[i].satpos, WGS84) * R2D, satres[i].elev * R2D);
					}
					//转成ENU坐标
					rcvposenu = XYZ2ENU(BESTPOS, res.rcvpos, WGS84);
					fprintf(ostream, "success %f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %10.3f %10.3f %3.3f %3.3f %7.3f %7.3f %7.3f %7.3f %3.3f %3d %3d %3d\n",
						range.time.sow,
						res.rcvpos.xyz.c1, res.rcvpos.xyz.c2, res.rcvpos.xyz.c3,
						res.rcvposblh.blh.c1, res.rcvposblh.blh.c2, res.rcvposblh.blh.c3,
						rcvposenu.denu[0], rcvposenu.denu[1], rcvposenu.denu[2],
						res.dtr[0], res.dtr[1], res.PDOP, res.sigma,
						res.vx, res.vy, res.vz, res.vdtr, res.velsigma,
						res.satnum[0], res.satnum[1], range.sat_num);
				}
				if (res.Status == false)
				{
					fprintf(ostream, "fail %f\n", range.time.sow);
				}
				res = RcvRes();
				for (int i = 0; i < MAXSATNUM; i++)
				{
					satres[i] = SatRes();
				}
			}
			else
			{
				fprintf(ostream, "%f SOW:SPP failed\n", range.time.sow);
			}
			temp_range = range;
			range = Range();
		}
	}
	fclose(stream);
	fclose(ostream);
	fclose(satstream);
}
void NetWork()
{
	SOCKET server;
	bool Status = false, flag = false;
	bool isOpen = OpenSocket(server, "47.114.134.129", 7190);
	int LenRead = 0;
	int LenRem = 0;
	unsigned char Buff[MAXNOVDLEN];
	Range range;
	Range temp_range;
	RcvRes res;
	RcvRes init;
	SatRes satres[MAXSATNUM];
	GPSEph geph[MAXGPSSAT];
	BDSEph beph[MAXBDSSAT];
	FILE* file = fopen("data.txt", "w");
	if (isOpen)
	{
		while (true)
		{
			LenRead = recv(server, (char*)Buff, MAXNOVDLEN, 0);
			if (LenRead > 0)
			{
				Status = DecodeNovOem7(Buff, LenRem, &range, geph, beph, 0);
				if (Status)
				{
					// 进行粗差探测
                    DetectOutlier(&range, &temp_range);
					// SPP
					flag = SPP(&res, satres, &init, &range, geph, beph);
					if (flag)
					{
						if (res.Status)
						{
							SPV(&res, satres, &range);
							init = res;
							printf("success %f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %10.3f %10.3f %3.3f %3.3f %7.3f %7.3f %7.3f %7.3f %3.3f %3d %3d %3d\n",
								range.time.sow,
								res.rcvpos.xyz.c1, res.rcvpos.xyz.c2, res.rcvpos.xyz.c3,
								res.rcvposblh.blh.c1 * R2D, res.rcvposblh.blh.c2 * R2D, res.rcvposblh.blh.c3,
								res.dtr[0], res.dtr[1], res.PDOP, res.sigma,
								res.vx, res.vy, res.vz, res.vdtr, res.velsigma,
								res.satnum[0], res.satnum[1], range.sat_num);
							fprintf(file, "success %f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %10.3f %10.3f %3.3f %3.3f %7.3f %7.3f %7.3f %7.3f %3.3f %3d %3d %3d\n",
								range.time.sow,
								res.rcvpos.xyz.c1, res.rcvpos.xyz.c2, res.rcvpos.xyz.c3,
								res.rcvposblh.blh.c1 * R2D, res.rcvposblh.blh.c2 * R2D, res.rcvposblh.blh.c3,
								res.dtr[0], res.dtr[1], res.PDOP, res.sigma,
								res.vx, res.vy, res.vz, res.vdtr, res.velsigma,
								res.satnum[0], res.satnum[1], range.sat_num);
						}
						if (!res.Status)
						{
							printf("%f SOW:SPP failed\n", range.time.sow);
							fprintf(file, "%f SOW:SPP failed\n", range.time.sow);
						}
						res = RcvRes();
						for (int i = 0; i < MAXSATNUM; i++) {
							satres[i] = SatRes();
						}
					}
					else
					{
						printf("%f SOW:SPP failed\n", range.time.sow);
						fprintf(file, "%f SOW:SPP failed\n", range.time.sow);
					}
					temp_range = range;
					range = Range();
				}
			}
		}
	}
	fclose(file);
}

void OutRcvSol(const RcvRes& rcvres, const Range& range, std::ostream& os = std::cout)
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
void OutSatSol(const RcvRes& rcvres, const SatRes& satres, const Range& range, std::ostream& os = std::cout)
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
void OutRTKSol(const GPSTime& b_time, const GPSTime& r_time, const RcvRes& rcv_res, const DDObs& dd_obs, std::ostream& os = std::cout)
{
    std:string fmt;
    fmt += "{:10.3f} {:10.3f} ";
    fmt += "{:14.3f} {:14.3f} {:14.3f} ";
    fmt += "{:10.3f} {:10.3f} {:10.3f} ";
    fmt += "{:6.3f} {:6.3f} ";
    fmt += "{:3d} {:3d}\n";
    auto args = std::make_format_args(
        r_time.sow, b_time.sow,
        rcv_res.rcvpos.dxyz[0], rcv_res.rcvpos.dxyz[1], rcv_res.rcvpos.dxyz[2],
        dd_obs.dpos[0], dd_obs.dpos[1], dd_obs.dpos[2],
        rcv_res.sigma, rcv_res.PDOP,
        dd_obs.dd_sat_num[0], dd_obs.dd_sat_num[1]
    );
    os << std::vformat(fmt, args);
}
void sppPos(int mode, std::ostream& rcvout = std::cout, std::ostream& satout = std::cout)
{
	int LenRead = 0, LenRem = 0;
	bool Status = false, flag = false;
	unsigned char Buff[MAXNOVDLEN];
	Range range, temp_range;
	RcvRes res, init;
	SatRes satres[MAXSATNUM];
	GPSEph geph[MAXGPSSAT];
	BDSEph beph[MAXBDSSAT];
	if (mode == 1) {// network
		SOCKET server;
		if (!OpenSocket(server, "47.114.134.129", 7190)) {
			std::cout << "server down" << std::endl;
			return;
		}
		while (true) {
			LenRead = recv(server, (char*)Buff, MAXNOVDLEN, 0);
			if (LenRead > 0) {
				Status = DecodeNovOem7(Buff, LenRem, &range, geph, beph, 0);// Decode
				if (Status) {
                    DetectOutlier(&range, &temp_range);// Outlier Detect
					flag = SPP(&res, satres, &init, &range, geph, beph);// SPP
					if (flag && res.Status) {
						SPV(&res, satres, &range);
						init = res;
						for (int i = 0; i < res.satnum[0] + res.satnum[1]; i++) {
							OutSatSol(res, satres[i], range, satout);
						}
					}
					OutRcvSol(res, range, rcvout);
					res = RcvRes();
					for (int i = 0; i < MAXSATNUM; i++) {
						satres[i] = SatRes();
					}
					temp_range = range;
					range = Range();
				}
			}
		}
	}
	else if (mode == 0) {// input file
		FILE* stream;
		if (fopen_s(&stream, "data\\oem719-202202080900-2.bin", "r+b")) {
			std::cout << "file open error" << std::endl;
			return;
		}
		while (feof(stream) == 0) {
			LenRead = fread_s(Buff + LenRem, sizeof(unsigned char) * MAXNOVDLEN, sizeof(unsigned char), MAXNOVDLEN - LenRem, stream);
			if (LenRead < (MAXNOVDLEN - LenRem)) {
				std::cout << "finished!" << std::endl;
				fclose(stream);
				return;
			}
			Status = DecodeNovOem7(Buff, LenRem, &range, geph, beph, 0);// Decode
			if (Status) {
                DetectOutlier(&range, &temp_range);// Outlier Detect
				flag = SPP(&res, satres, &init, &range, geph, beph);// SPP
				if (flag && res.Status) {
					SPV(&res, satres, &range);
					init = res;
					for (int i = 0; i < res.satnum[0] + res.satnum[1]; i++) {
						OutSatSol(res, satres[i], range, satout);
					}
				}
				OutRcvSol(res, range, rcvout);
				res = RcvRes();
				for (int i = 0; i < MAXSATNUM; i++)
				{
					satres[i] = SatRes();
				}
				temp_range = range;
				range = Range();
			}
		}
	}
	else {
		std::cout << "unknown mode" << std::endl;
		return;
	}
}
void rtkPos(int mode)
{
    int r_LenRead{0}, r_LenRem{0};
    int b_LenRead{0}, b_LenRem{0};
    Sync sync{Sync::UNK};
    bool r_status{false}, b_status{false};
    bool r_is_spp{false}, b_is_spp{false};
    bool b_eof{true};
    unsigned char r_buff[MAXNOVDLEN], b_buff[MAXNOVDLEN];
    Range r_range, r_temp_range;
    Range b_range, b_temp_range;
    RcvRes r_res, r_init;
    RcvRes b_res, b_init;
    SatRes b_sat_res[MAXSATNUM], r_sat_res[MAXSATNUM];
    GPSEph geph[MAXGPSSAT];
    BDSEph beph[MAXBDSSAT];
    SDEpochObs sd_obs, sd_temp_obs;
    DDObs dd_obs;
    std::vector<std::tuple<Range&, Range&, RcvRes&, RcvRes&, SatRes*, bool&>> spp_data;
    spp_data.emplace_back(b_range, b_temp_range, b_res, b_init, b_sat_res, b_is_spp);
    spp_data.emplace_back(r_range, r_temp_range, r_res, r_init, r_sat_res, r_is_spp);
    if (mode == 1) {// real-time
        SOCKET b_server, r_server;
        if (!OpenSocket(b_server, "47.114.134.129", 7190) || !OpenSocket(r_server, "8.140.46.126", 3002)) {
            std::cout << "server connect error" << std::endl;
            return;
        }
        while (true) {
            if (sync != Sync::ROV) {
                r_LenRead = recv(r_server, (char*)r_buff+r_LenRem, MAXNOVDLEN-r_LenRem, 0);
                r_LenRem=r_LenRem+r_LenRead;
            }
            if (sync != Sync::BAS) {
                b_LenRead = recv(b_server, (char*)b_buff+b_LenRem, MAXNOVDLEN-b_LenRem, 0);
                b_LenRem=b_LenRem+b_LenRead;
            }
            if (!r_LenRead || !b_LenRead) continue;
            if (sync != Sync::BAS) b_status = DecodeNovOem7(b_buff, b_LenRem, &b_range, geph, beph, 1);
            if (sync != Sync::ROV) r_status = DecodeNovOem7(r_buff, r_LenRem, &r_range, geph, beph, 1);
            sync = timeSync(b_range, b_status, r_range, r_status);
            if (sync != Sync::SYN) continue;
            printf("base: %lf rover: %lf b_status: %d r_status: %d \n",b_range.time.sow,r_range.time.sow,b_status,r_status);
            std::cout << std::format("Sync at {}\n", b_range.time.sow);
            for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                DetectOutlier(&range, &temp_range);
                is_spp = SPP(&res, sat_res, &init, &range, geph, beph);
            }
            if (b_is_spp && r_is_spp) {
                b_init = b_res;
                r_init = r_res;
                FormSDEpochObs(b_range, r_range, sd_obs);
                DetectCycleSlip(sd_temp_obs, sd_obs);
                if (DetermineRefSat(b_range, r_range, b_sat_res, r_sat_res, sd_obs, dd_obs)) {
                    int sat_num = dd_obs.dd_sat_num[0]+dd_obs.dd_sat_num[1];
                    auto a = std::make_unique<double[]>(2*sat_num);
                    auto Q_a = std::make_unique<double[]>(4*sat_num*sat_num);
                    RTKFloat(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs, a.get(), Q_a.get());
                    std::cout << "Float: ";
                    OutRTKSol(b_range.time, r_range.time, r_res, dd_obs);
                    lambda(2*sat_num, 2, a.get(), Q_a.get(), dd_obs.fixed_amb, dd_obs.res_amb);
                    dd_obs.ratio = dd_obs.res_amb[1] / dd_obs.res_amb[0];
                    RTKFixed(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs);
                    std::cout << "Fixed: ";
                    OutRTKSol(b_range.time, r_range.time, r_res, dd_obs);
                }
            }
            for (auto& sat : b_sat_res) {
                sat = SatRes();
            }
            for (auto& sat : r_sat_res) {
                sat = SatRes();
            }
            sd_temp_obs = sd_obs;
            for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                res = RcvRes();
                temp_range = range;
                range = Range();
                is_spp = false;
            }
        }
    }
    else if (mode == 0) {// file postprocess
        FILE *r_file, *b_file;
        if (
                fopen_s(&b_file, "data\\Zero-baseline\\oem719-202203031500-1.bin", "r+b") ||
            fopen_s(&r_file, "data\\Zero-baseline\\oem719-202203031500-2.bin", "r+b")) {
//            fopen_s(&b_file, "data\\Short-baseline\\oem719-202203170900-1.bin", "r+b") ||
//            fopen_s(&r_file, "data\\Short-baseline\\oem719-202203170900-2.bin", "r+b")) {
            std::cout << "file open error" << std::endl;
            return;
        }
        while (!feof(r_file) || !feof(b_file)) {
            if (sync != Sync::BAS) r_LenRead = fread_s(r_buff + r_LenRem, sizeof(unsigned char) * MAXNOVDLEN, sizeof(unsigned char), MAXNOVDLEN - r_LenRem, r_file);
            if (sync != Sync::ROV && b_eof) b_LenRead = fread_s(b_buff + b_LenRem, sizeof(unsigned char) * MAXNOVDLEN, sizeof(unsigned char), MAXNOVDLEN - b_LenRem, b_file);
            if (r_LenRead < (MAXNOVDLEN - r_LenRem)) {
                std::cout << "finished!" << std::endl;
                fclose(r_file);
                return;
            }
            if (b_LenRead < (MAXNOVDLEN - b_LenRem)) {
                fclose(b_file);
                b_eof = false;
                continue;
            }
            if (sync != Sync::BAS && b_eof) b_status = DecodeNovOem7(b_buff, b_LenRem, &b_range, geph, beph, 0);// Decode base
            if (sync != Sync::ROV) r_status = DecodeNovOem7(r_buff, r_LenRem, &r_range, geph, beph, 0);// Decode rover
            sync = timeSync(b_range, b_status, r_range, r_status);
            if (sync == Sync::SYN) {
                std::cout << "Time Synchronized in " << b_range.time.sow << std::endl;
                for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                    DetectOutlier(&range, &temp_range);
                    is_spp = SPP(&res, sat_res, &init, &range, geph, beph);
                }
                if (b_is_spp && r_is_spp) {
                    b_init = b_res;
                    r_init = r_res;
                    FormSDEpochObs(b_range, r_range, sd_obs);
                    DetectCycleSlip(sd_temp_obs, sd_obs);
                    if (DetermineRefSat(b_range, r_range, b_sat_res, r_sat_res, sd_obs, dd_obs)) {
                        int sat_num = dd_obs.dd_sat_num[0]+dd_obs.dd_sat_num[1];
                        auto a = std::make_unique<double[]>(2*sat_num);
                        auto Q_a = std::make_unique<double[]>(4*sat_num*sat_num);
                        RTKFloat(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs, a.get(), Q_a.get());
                        std::cout << "Float: ";
                        OutRTKSol(b_range.time, r_range.time, r_res, dd_obs);
                        lambda(2*sat_num, 2, a.get(), Q_a.get(), dd_obs.fixed_amb, dd_obs.res_amb);
                        dd_obs.ratio = dd_obs.res_amb[1] / dd_obs.res_amb[0];
                        RTKFixed(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs);
                        std::cout << "Fixed: ";
                        OutRTKSol(b_range.time, r_range.time, r_res, dd_obs);
                    }
                }
                for (auto& sat : b_sat_res) {
                    sat = SatRes();
                }
                for (auto& sat : r_sat_res) {
                    sat = SatRes();
                }
                sd_temp_obs = sd_obs;
                for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                    res = RcvRes();
                    temp_range = range;
                    range = Range();
                    is_spp = false;
                }
            }
        }
    }
}
int main()
{
	std::cout << "RealTime Kinematic Solution:" << std::endl;
    // sppPos(0);
    rtkPos(1);
}
