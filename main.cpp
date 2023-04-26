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

void sppPos(const Config& config, std::ostream& rcvout = std::cout, std::ostream& satout = std::cout)
{
	int LenRead = 0, LenRem = 0;
	bool Status = false, flag = false;
	unsigned char Buff[MAXNOVDLEN];
	Range range, temp_range;
	RcvRes res, init;
	SatRes satres[MAXSATNUM];
	GPSEph geph[MAXGPSSAT];
	BDSEph beph[MAXBDSSAT];
	if (config.online) {// network
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
					flag = SPP(&res, satres, &init, &range, geph, beph, config);// SPP
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
	else {// input file
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
				flag = SPP(&res, satres, &init, &range, geph, beph, config);// SPP
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
}
void rtkPos(const Config& config)
{
    int r_LenRead{0}, r_LenRem{0};
    int b_LenRead{0}, b_LenRem{0};
    int epoch_float{0}, epoch_fixed{0};
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
    std::ofstream float_file("data//Output//20230426-float.csv", std::ios::out);
    std::ofstream fixed_file("data//Output//20230426-fixed.csv", std::ios::out);
    if (config.out_mode == 1) {
        if (!float_file.is_open() || !fixed_file.is_open()) {
            std::cout << "Could not open file!" << std::endl;
            return;
        }
        if (config.out_head) {
            std::string head = "rover time(s),base time(s),";
            head += config.sol_format == 2 ? "E(m),N(m),U(m)," : (config.sol_format == 0 ? "B(rad),L(rad),H(m)," : "X(m),Y(m),Z(m),");
            head += "baseline0(m),baseline1(m),baseline2(m),";
            head += "sigma,PDOP,HDOP,VDOP,";
            head += "GPS satnum,BDS satnum,ratio\n";
            float_file << head;
            fixed_file << head;
        }
    }
    if (config.online) {// real-time
        SOCKET b_server, r_server;
        if (!OpenSocket(b_server, config.IP_b.c_str(), config.port_b) || !OpenSocket(r_server, config.IP_r.c_str(), config.port_r)) {
            std::cout << "server connect error" << std::endl;
            return;
        }
        while (true) {
            if (sync != Sync::ROV) {
                r_LenRead = recv(r_server, (char*)r_buff+r_LenRem, MAXNOVDLEN-r_LenRem, 0);
                r_LenRem=r_LenRem + r_LenRead;
            }
            if (sync != Sync::BAS) {
                b_LenRead = recv(b_server, (char*)b_buff+b_LenRem, MAXNOVDLEN-b_LenRem, 0);
                b_LenRem=b_LenRem + b_LenRead;
            }
            if (!r_LenRead || !b_LenRead) continue;
            if (sync != Sync::BAS) b_status = DecodeNovOem7(b_buff, b_LenRem, &b_range, geph, beph, 1);
            if (sync != Sync::ROV) r_status = DecodeNovOem7(r_buff, r_LenRem, &r_range, geph, beph, 1);
            sync = timeSync(b_range, b_status, r_range, r_status);
            if (sync != Sync::SYN) continue;
            std::cout << std::format("Sync at {}\n", b_range.time.sow);
            for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                DetectOutlier(&range, &temp_range);
                is_spp = SPP(&res, sat_res, &init, &range, geph, beph, config);
            }
            if (b_is_spp && r_is_spp) {
                b_init = b_res;
                r_init = r_res;
                FormSDEpochObs(b_range, r_range, sd_obs);
                DetectCycleSlip(sd_temp_obs, sd_obs);
                if (DetermineRefSat(b_range, r_range, b_sat_res, r_sat_res, sd_obs, dd_obs, config)) {
                    int sat_num = dd_obs.dd_sat_num[0]+dd_obs.dd_sat_num[1];
                    auto a = std::make_unique<double[]>(2*sat_num);
                    auto Q_a = std::make_unique<double[]>(4*sat_num*sat_num);
                    RTKFloat(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs, a.get(), Q_a.get());
                    std::cout << "Float: ";
                    OutRTKSol(b_range.time, r_range.time, r_res, dd_obs, config, float_file);
                    epoch_float++;
                    auto info = lambda(2*sat_num, 2, a.get(), Q_a.get(), dd_obs.fixed_amb, dd_obs.res_amb);
                    dd_obs.ratio = dd_obs.res_amb[1] / dd_obs.res_amb[0];
                    if (dd_obs.ratio > config.ratio_thres && info != -1) {
                        RTKFixed(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs);
                        std::cout << "Fixed: ";
                        OutRTKSol(b_range.time, r_range.time, r_res, dd_obs, config, fixed_file);
                        epoch_fixed++;
                    }
                    std::cout << std::format("fixed rate: {:7.4f}\n", static_cast<double>(epoch_fixed)/static_cast<double>(epoch_float));
                }
            }
            for (auto& sat : b_sat_res) {
                sat = SatRes();
            }
            for (auto& sat : r_sat_res) {
                sat = SatRes();
            }
            sd_temp_obs = sd_obs;
            dd_obs = DDObs();
            for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                res = RcvRes();
                temp_range = range;
                range = Range();
                is_spp = false;
            }
        }
    }
    else {// file postprocess
        FILE *r_file, *b_file;
        if (fopen_s(&b_file, config.infile_b.c_str(), "r+b") ||
        fopen_s(&r_file, config.infile_r.c_str(), "r+b")) {
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
                    is_spp = SPP(&res, sat_res, &init, &range, geph, beph, config);
                }
                if (b_is_spp && r_is_spp) {
                    b_init = b_res;
                    r_init = r_res;
                    FormSDEpochObs(b_range, r_range, sd_obs);
                    DetectCycleSlip(sd_temp_obs, sd_obs);
                    if (DetermineRefSat(b_range, r_range, b_sat_res, r_sat_res, sd_obs, dd_obs, config)) {
                        int sat_num = dd_obs.dd_sat_num[0]+dd_obs.dd_sat_num[1];
                        auto a = std::make_unique<double[]>(2*sat_num);
                        auto Q_a = std::make_unique<double[]>(4*sat_num*sat_num);
                        RTKFloat(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs, a.get(), Q_a.get());
                        std::cout << "Float: ";
                        OutRTKSol(b_range.time, r_range.time, r_res, dd_obs, config, float_file);
                        lambda(2*sat_num, 2, a.get(), Q_a.get(), dd_obs.fixed_amb, dd_obs.res_amb);
                        dd_obs.ratio = dd_obs.res_amb[1] / dd_obs.res_amb[0];
                        RTKFixed(r_res, b_sat_res, r_sat_res, sd_obs, dd_obs);
                        std::cout << "Fixed: ";
                        OutRTKSol(b_range.time, r_range.time, r_res, dd_obs, config, fixed_file);
                    }
                }
                for (auto& sat : b_sat_res) {
                    sat = SatRes();
                }
                for (auto& sat : r_sat_res) {
                    sat = SatRes();
                }
                sd_temp_obs = sd_obs;
                dd_obs = DDObs();
                for (auto& [range, temp_range, res, init, sat_res, is_spp] : spp_data) {
                    res = RcvRes();
                    temp_range = range;
                    range = Range();
                    is_spp = false;
                }
            }
        }
    }
    float_file.close();
    fixed_file.close();
}
int main()
{
    Config config;
	std::cout << "RealTime Kinematic Solution:" << std::endl;
    ReadConfigureFile("option.conf", config);
    config.pos_mode == 0 ? sppPos(config) : rtkPos(config);
    std::cout << "End" << std::endl;
    // sppPos(0);
    // rtkPos(1);
}
