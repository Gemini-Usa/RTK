#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <utility>
#include "constants.h"
#include "function.h"

/// <summary>
/// 标准单点定位
/// </summary>
/// <param name="res">待求接收机结果</param>
/// <param name="init">上一次接收机结果</param>
/// <param name="range">观测值</param>
/// <param name="geph">GPS星历</param>
/// <param name="beph">北斗星历</param>
/// <param name="satres">待求卫星结果</param>
/// <returns>true:解算成功 false:解算失败</returns>
bool SPP(RcvRes* res, SatRes satres[], const RcvRes* init, const Range* range, const GPSEph geph[], const BDSEph beph[])
{
	if (range->sat_num < 4) return false;
	double poserr = 1E-6; //定位误差
	int iterator = 0, prn, gsatnum = 0, bsatnum = 0;
	GPSTime t = range->time, ttr;
	NavSys sys;
	RcvRes resk;
	LSQInput input;
	std::pair<SatRes, XYZ> gres[MAXGPSSAT];
	std::pair<SatRes, XYZ> bres[MAXBDSSAT];
	std::vector<LSQInput> lsqinput;
	if (init->Status) *res = *init;
	
	for (resk.rcvpos.dxyz[0] = 1; GetDist(res->rcvpos, resk.rcvpos) > poserr; iterator++)
	{
		gsatnum = bsatnum = 0;
		resk = *res;
		lsqinput.clear();
		//计算或排除每颗卫星的卫星位置、钟差
		for (int n = 0; n < range->sat_num; n++)
		{
			if (!range->obs[n].Status) continue;

			prn = range->obs[n].prn;
			sys = range->obs[n].sys;
			switch (sys)
			{
			case GPS:
				if (!geph[prn - 1].Status || geph[prn - 1].Health) continue;
				//计算卫星钟差、卫星位置
				gres[prn - 1].first.dts = 0;
				ttr.week = t.week;
				ttr.sow = t.sow - (range->obs[n].combination.P_IF / CLIGHT);
				GPSGetSatState(range->obs[n], geph, gres[prn - 1].first, ttr);
				ttr.sow -= gres[prn - 1].first.dts;
				GPSGetSatState(range->obs[n], geph, gres[prn - 1].first, ttr);
				if (!gres[prn - 1].first.Status) continue;
				//将第一次计算结果存储起来
				if (iterator == 0) gres[prn - 1].second = gres[prn - 1].first.satpos;
				//地球自转改正
				EarthRotate(gres[prn - 1].second.dxyz, gres[prn - 1].first.satpos.dxyz, GPS, (GetDist(gres[prn - 1].second, res->rcvpos) / CLIGHT));
				//高度角计算,若|elevation|<15°,剔除该卫星
				res->rcvposblh = XYZ2BLH(res->rcvpos, WGS84);
				gres[prn - 1].first.elev = GetElev(res->rcvpos, gres[prn - 1].first.satpos, WGS84);
				if (fabs(gres[prn - 1].first.elev) < 15 * D2R)
				{
					gres[prn - 1].first.Status = false;
					continue;
				}
				//对流层改正
				gres[prn - 1].first.trop = Hopfield(res->rcvposblh.blh.c3, gres[prn - 1].first.elev);
				gres[prn - 1].first.tgd = 0.0;
				//中间变量存储
				input.SysType = GPS;
				input.satpos = gres[prn - 1].first.satpos;
				input.satclk = gres[prn - 1].first.dts;
				input.P_IF = range->obs[n].combination.P_IF;
				input.trop = gres[prn - 1].first.trop;
				input.tgd = 0.0;
				lsqinput.push_back(input);

				gsatnum++;
				break;
			case BDS:
				if (!beph[prn - 1].Status || beph[prn - 1].Health) continue;
				//tgd改正
				bres[prn - 1].first.tgd = (CLIGHT * beph[prn - 1].tgd[0]) / ((FREQ_B3 / FREQ_B1) * (FREQ_B3 / FREQ_B1) - 1);
				//计算卫星钟差、卫星位置
				bres[prn - 1].first.dts = 0;
				ttr.week = t.week;
				ttr.sow = t.sow - ((range->obs[n].combination.P_IF + bres[prn - 1].first.tgd) / CLIGHT);
				BDSGetSatState(range->obs[n], beph, bres[prn - 1].first, ttr);
				ttr.sow -= bres[prn - 1].first.dts;
				BDSGetSatState(range->obs[n], beph, bres[prn - 1].first, ttr);
				if (!bres[prn - 1].first.Status) continue;
				//将第一次计算结果存储起来
				if (iterator == 0) bres[prn - 1].second = bres[prn - 1].first.satpos;
				//地球自转改正
				EarthRotate(bres[prn - 1].second.dxyz, bres[prn - 1].first.satpos.dxyz, BDS, (GetDist(bres[prn - 1].second, res->rcvpos) / CLIGHT));
				res->rcvposblh = XYZ2BLH(res->rcvpos, WGS84);
				bres[prn - 1].first.elev = GetElev(res->rcvpos, bres[prn - 1].first.satpos, WGS84);
				//高度角计算,若|elevation|<15°,剔除该卫星
				if (fabs(bres[prn - 1].first.elev) < 15 * D2R)
				{
					bres[prn - 1].first.Status = false;
					continue;
				}
				//对流层改正
				bres[prn - 1].first.trop = Hopfield(res->rcvposblh.blh.c3, bres[prn - 1].first.elev);
				//中间变量存储
				input.SysType = BDS;
				input.satpos = bres[prn - 1].first.satpos;
				input.satclk = bres[prn - 1].first.dts;
				input.P_IF = range->obs[n].combination.P_IF;
				input.trop = bres[prn - 1].first.trop;
				input.tgd = bres[prn - 1].first.tgd;
				lsqinput.push_back(input);

				bsatnum++;
				break;
			default:continue;
			}
		}
		//最小二乘解算
		LSQ(lsqinput, res, gsatnum, bsatnum);
        if (fabs(Norm(res->rcvpos.dxyz, 3)) < 1.0E-10) return false;
		res->satnum[0] = gsatnum;
		res->satnum[1] = bsatnum;

		if (iterator > 10) return false;
	}
	gsatnum = 0;
	bsatnum = 0;
	//卫星结果储存
	for (int i = 0; i < MAXGPSSAT; i++)
	{
		if (!gres[i].first.Status) continue;
		satres[gsatnum] = gres[i].first;
		satres[gsatnum].sys = GPS;
		gsatnum++;
	}
	for (int i = 0; i < MAXBDSSAT; i++)
	{
		if (!bres[i].first.Status) continue;
		satres[gsatnum + bsatnum] = bres[i].first;
		satres[gsatnum + bsatnum].sys = BDS;
		bsatnum++;
	}
	return true;
}

/// <summary>
/// 最小二乘解
/// </summary>
/// <param name="input">中间变量</param>
/// <param name="rcvres">待求接收机结果</param>
/// <param name="gsatnum">GPS卫星数量</param>
/// <param name="bsatnum">北斗卫星数量</param>
void LSQ(std::vector<LSQInput>& input, RcvRes* rcvres, int gsatnum, int bsatnum)
{
	int sysnum;
	double pos[3];
	double los[3];
	if (gsatnum == 0 || bsatnum == 0)
	{
		if (gsatnum == 0 && bsatnum >= 4) sysnum = 1;
		else if (bsatnum == 0 && gsatnum >= 4) sysnum = 1;
		else
		{
			rcvres->Status = false;
			return;
		}
	}
	else if (gsatnum + bsatnum >= 5) sysnum = 2;
	else
	{
		rcvres->Status = false;
		return;
	}

	int obsnum = gsatnum + bsatnum;
	auto* B = new double[obsnum * (3 + sysnum)];
	auto* w = new double[obsnum];
	auto* BT = new double[obsnum * (3 + sysnum)];
	auto* NBB = new double[(3 + sysnum) * (3 + sysnum)];
	auto* invNBB = new double[(3 + sysnum) * (3 + sysnum)];
	auto* W = new double[(3 + sysnum)];
	auto* hatx = new double[(3 + sysnum)];
	auto* Bx = new double[obsnum];
	auto* v = new double[obsnum];
	auto* vT = new double[obsnum];
	auto* vTv = new double[1];

	for (int i = 0; i < obsnum; i++)
	{
		for (int j = 0; j < 3 + sysnum; j++) B[i * (3 + sysnum) + j] = 0.0;
		VectorSub(input[i].satpos.dxyz, rcvres->rcvpos.dxyz, pos, 3);
		VectorNormalize(pos, los, 3);
		B[i * (3 + sysnum)] = -los[0];
		B[i * (3 + sysnum) + 1] = -los[1];
		B[i * (3 + sysnum) + 2] = -los[2];
		B[i * (3 + sysnum) + 3 + input[i].SysType] = 1.0;
		w[i] = input[i].P_IF + input[i].tgd - (GetDist(input[i].satpos, rcvres->rcvpos) + rcvres->dtr[input[i].SysType] - CLIGHT * input[i].satclk + input[i].trop);
	}
	/*printf("w:\n");
	MatrixPrint(w, obsnum, 1);*/
	MatrixTranspose(B, BT, obsnum, (3 + sysnum));
	MatrixMultiply(BT, B, NBB, (3 + sysnum), obsnum, (3 + sysnum));
	MatrixInverse(NBB, invNBB, (3 + sysnum));
	MatrixMultiply(BT, w, W, (3 + sysnum), obsnum, 1);
	MatrixMultiply(invNBB, W, hatx, (3 + sysnum), (3 + sysnum), 1);
	MatrixMultiply(B, hatx, Bx, obsnum, (3 + sysnum), 1);
	MatrixSub(Bx, w, v, obsnum, 1);
	MatrixTranspose(v, vT, obsnum, 1);
	MatrixMultiply(vT, v, vTv, 1, obsnum, 1);

	for (int i = 0; i < 3; i++)
	{
		rcvres->rcvpos.dxyz[i] += hatx[i];
	}
	for (int i = 0; i < sysnum; i++)
	{
		rcvres->dtr[i] += hatx[3 + i];
	}
	rcvres->sigma = sqrt(vTv[0] / (obsnum - 3.0 - sysnum));
	rcvres->PDOP = sqrt(invNBB[0] + invNBB[3 + sysnum + 1] + invNBB[2 * (3 + sysnum) + 2]);
	rcvres->rcvposblh = XYZ2BLH(rcvres->rcvpos, WGS84);
	rcvres->Status = true;

	delete[] B;
	delete[] w;
	delete[] BT;
	delete[] NBB;
	delete[] invNBB;
	delete[] W;
	delete[] hatx;
	delete[] Bx;
	delete[] v;
	delete[] vT;
	delete[] vTv;
	B = w = BT = NBB = invNBB = W = hatx = Bx = v = vT = vTv = nullptr;
}

/// <summary>
/// 标准单点测速
/// </summary>
/// <param name="rcvres">待求接收机结果</param>
/// <param name="satres">卫星结果</param>
/// <param name="range">观测值</param>
void SPV(RcvRes* rcvres, const SatRes satres[], const Range* range)
{
	int obsnum = rcvres->satnum[0] + rcvres->satnum[1];
	double pos[3];
	double los[3];
	double vel[3], velECEF[3];
	double lambda[2] = { CLIGHT / FREQ_L1, CLIGHT / FREQ_B1 };
	double dopp = 0.0;
	double tau = 0.0;

	double* B = new double[obsnum * 4];
	double* w = new double[obsnum];
	double* BT = new double[obsnum * 4];
	double NBB[4 * 4];
	double invNBB[4 * 4];
	double W[4];
	double hatx[4];
	double* Bx = new double[obsnum];
	double* v = new double[obsnum];
	double* vT = new double[obsnum];
	double vTv[1];
	for (int i = 0; i < obsnum; i++)
	{
		if (!satres[i].Status) continue;
		for (int j = 0; j < range->sat_num; j++)
		{
			if (range->obs[j].prn == satres[i].prn && range->obs[j].sys == satres[i].sys)
			{
				dopp = range->obs[j].D[0];
				break;
			}
		}
		vel[0] = satres[i].vx;
		vel[1] = satres[i].vy;
		vel[2] = satres[i].vz;
		//地球自转改正
		EarthRotate(vel, velECEF, satres[i].sys, GetDist(satres[i].satpos, rcvres->rcvpos) / CLIGHT);
		for (int j = 0; j < 4; j++) B[i * 4 + j] = 0.0;
		VectorSub(satres[i].satpos.dxyz, rcvres->rcvpos.dxyz, pos, 3);
		VectorNormalize(pos, los, 3);
		B[i * 4] = -los[0];
		B[i * 4 + 1] = -los[1];
		B[i * 4 + 2] = -los[2];
		B[i * 4 + 3] = 1.0;
		w[i] = dopp * lambda[satres[i].sys] - (VectorInnerProduct(los, vel, 3) - CLIGHT * satres[i].dv);
	}
	//MatrixPrint(w, obsnum, 1);
	MatrixTranspose(B, BT, obsnum, 4);
	MatrixMultiply(BT, B, NBB, 4, obsnum, 4);
	MatrixInverse(NBB, invNBB, 4);
	MatrixMultiply(BT, w, W, 4, obsnum, 1);
	MatrixMultiply(invNBB, W, hatx, 4, 4, 1);
	MatrixMultiply(B, hatx, Bx, obsnum, 4, 1);
	MatrixSub(Bx, w, v, obsnum, 1);
	MatrixTranspose(v, vT, obsnum, 1);
	MatrixMultiply(vT, v, vTv, 1, obsnum, 1);
	rcvres->vx = hatx[0];
	rcvres->vy = hatx[1];
	rcvres->vz = hatx[2];
	rcvres->vdtr = hatx[3];
	rcvres->velsigma = sqrt(vTv[0]);

	delete[] B;
	delete[] w;
	delete[] BT;
	delete[] Bx;
	delete[] v;
	delete[] vT;
	B = w = BT = Bx = v = vT = nullptr;
}