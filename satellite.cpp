#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "function.h"

/// <summary>
/// 计算GPS卫星位置
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="GEph">GPS星历</param>
/// <param name="GRes">传出GPS结果</param>
/// <param name="t">信号发射时刻</param>
void GPSGetSatpos(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t)
{
	double e, n0, n, Mk, Ek, Ek0, vk, sinvk, cosvk, Phik, del_uk, del_rk, del_ik, uk, rk, ik, xkp, ykp, OMGk;
	double F = (-2.0 * sqrt(GM_GPS)) / (CLIGHT * CLIGHT);
	double itethres = 1.0E-10; //迭代阈值
	double tk;

	tk = (t.week - (GEph + prn - 1)->week) * 604800.0 + t.sow - (GEph + prn - 1)->toc;
	n0 = sqrt(GM_GPS / pow((GEph + prn - 1)->A, 3));
	n = n0 + (GEph + prn - 1)->del_n;
	Mk = (GEph + prn - 1)->M0 + n * tk;
	e = (GEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}
	sinvk = (sqrt(1 - e * e) * sin(Ek)) / (1 - e * cos(Ek));
	cosvk = (cos(Ek) - e) / (1 - e * cos(Ek));
	vk = atan2(sinvk, cosvk);
	Phik = vk + (GEph + prn - 1)->omg;
	del_uk = (GEph + prn - 1)->cus * sin(2 * Phik) + (GEph + prn - 1)->cuc * cos(2 * Phik);
	del_rk = (GEph + prn - 1)->crs * sin(2 * Phik) + (GEph + prn - 1)->crc * cos(2 * Phik);
	del_ik = (GEph + prn - 1)->cis * sin(2 * Phik) + (GEph + prn - 1)->cic * cos(2 * Phik);
	uk = Phik + del_uk;
	rk = (GEph + prn - 1)->A * (1 - e * cos(Ek)) + del_rk;
	ik = (GEph + prn - 1)->I0 + del_ik + (GEph + prn - 1)->Idot * tk;
	xkp = rk * cos(uk);
	ykp = rk * sin(uk);
	OMGk = (GEph + prn - 1)->OMG0 + ((GEph + prn - 1)->OMGdot - OMGE_GPS) * tk - OMGE_GPS * (GEph + prn - 1)->toe;
	GRes.satpos.xyz.c1 = xkp * cos(OMGk) - ykp * cos(ik) * sin(OMGk);
	GRes.satpos.xyz.c2 = xkp * sin(OMGk) + ykp * cos(ik) * cos(OMGk);
	GRes.satpos.xyz.c3 = ykp * sin(ik);
}

/// <summary>
/// 计算GPS卫星钟差
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="GEph">GPS星历</param>
/// <param name="GRes">传出GPS结果</param>
/// <param name="t">信号发射时刻</param>
void GPSGetSatclk(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t)
{
	double e, n0, n, Mk, Ek, Ek0, trel;
	double F = (-2 * sqrt(GM_GPS)) / (CLIGHT * CLIGHT);
	double itethres = 1.0E-10;//迭代阈值
	double tk;
	double del_t;

	tk = (t.week - (GEph + prn - 1)->week) * 604800.0 + t.sow - (GEph + prn - 1)->toc;
	n0 = sqrt(GM_GPS / pow((GEph + prn - 1)->A, 3));
	n = n0 + (GEph + prn - 1)->del_n;
	Mk = (GEph + prn - 1)->M0 + n * tk;
	e = (GEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}
	trel = F * e * sqrt((GEph + prn - 1)->A) * sin(Ek);
	del_t = gpst2Sec(t) - (double)(GEph + prn - 1)->week * 604800.0 - (GEph + prn - 1)->toc;
	GRes.dts = (GEph + prn - 1)->f0 + (GEph + prn - 1)->f1 * del_t + (GEph + prn - 1)->f2 * pow(del_t, 2) + trel;
}

/// <summary>
/// 计算GPS卫星速度和钟速
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="GEph">GPS星历</param>
/// <param name="GRes">传出GPS结果</param>
/// <param name="t">信号发射时刻</param>
void GPSGetSatvel(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t)
{
	double n0, n, Mk, Ek, Ek0, vk, sinvk, cosvk, Phik, del_uk, del_rk, del_ik, uk, rk, ik, xkp, ykp, OMGk;
	double e, Ekdot, Phikdot, ukdot, rkdot, ikdot, OMGkdot, xkpdot, ykpdot, del_trel;
	double F = (-2 * sqrt(GM_GPS)) / (CLIGHT * CLIGHT);
	double itethres = 1.0E-10;//迭代阈值
	double tk;
	double del_t;

	tk = (t.week - (GEph + prn - 1)->week) * 604800.0 + t.sow - (GEph + prn - 1)->toc;
	n0 = sqrt(GM_GPS / pow((GEph + prn - 1)->A, 3));
	n = n0 + (GEph + prn - 1)->del_n;
	Mk = (GEph + prn - 1)->M0 + n * tk;
	e = (GEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}
	sinvk = (sqrt(1 - e * e) * sin(Ek)) / (1 - e * cos(Ek));
	cosvk = (cos(Ek) - e) / (1 - e * cos(Ek));
	vk = atan2(sinvk, cosvk);
	Phik = vk + (GEph + prn - 1)->omg;
	del_uk = (GEph + prn - 1)->cus * sin(2 * Phik) + (GEph + prn - 1)->cuc * cos(2 * Phik);
	del_rk = (GEph + prn - 1)->crs * sin(2 * Phik) + (GEph + prn - 1)->crc * cos(2 * Phik);
	del_ik = (GEph + prn - 1)->cis * sin(2 * Phik) + (GEph + prn - 1)->cic * cos(2 * Phik);
	uk = Phik + del_uk;
	rk = (GEph + prn - 1)->A * (1 - e * cos(Ek)) + del_rk;
	ik = (GEph + prn - 1)->I0 + del_ik + (GEph + prn - 1)->Idot * tk;
	xkp = rk * cos(uk);
	ykp = rk * sin(uk);
	OMGk = (GEph + prn - 1)->OMG0 + ((GEph + prn - 1)->OMGdot - OMGE_GPS) * tk - OMGE_GPS * (GEph + prn - 1)->toe;

	Ekdot = n / (1 - e * cos(Ek));
	Phikdot = sqrt((1 + e) / (1 - e)) * pow((cos(vk / 2) / cos(Ek / 2)), 2) * Ekdot;
	ukdot = 2 * ((GEph + prn - 1)->cus * cos(2 * Phik) - (GEph + prn - 1)->cuc * sin(2 * Phik)) * Phikdot + Phikdot;
	rkdot = (GEph + prn - 1)->A * e * sin(Ek) * Ekdot + 2 * ((GEph + prn - 1)->crs * cos(2 * Phik) - (GEph + prn - 1)->crc * sin(2 * Phik)) * Phikdot;
	ikdot = (GEph + prn - 1)->Idot + 2 * ((GEph + prn - 1)->cis * cos(2 * Phik) - (GEph + prn - 1)->cic * sin(2 * Phik)) * Phikdot;
	OMGkdot = (GEph + prn - 1)->OMGdot - OMGE_GPS;
	double Rdot[12] = { cos(OMGk), -sin(OMGk) * cos(ik), -(xkp * sin(OMGk) + ykp * cos(OMGk) * cos(ik)), ykp * sin(OMGk) * sin(ik),
						sin(OMGk), cos(OMGk) * cos(ik), (xkp * cos(OMGk) - ykp * sin(OMGk) * cos(ik)), -ykp * cos(OMGk) * sin(ik),
						0, sin(ik), 0, ykp * cos(ik) };
	xkpdot = rkdot * cos(uk) - rk * ukdot * sin(uk);
	ykpdot = rkdot * sin(uk) + rk * ukdot * cos(uk);
	double c[4] = { xkpdot,
					ykpdot,
					OMGkdot,
					ikdot };
	double rdot[3] = { 0.0 };
	MatrixMultiply(Rdot, c, rdot, 3, 4, 1);
	GRes.vx = rdot[0];
	GRes.vy = rdot[1];
	GRes.vz = rdot[2];
	del_trel = F * e * sqrt((GEph + prn - 1)->A) * cos(Ek) * Ekdot;
	del_t = gpst2Sec(t) - (double)(GEph + prn - 1)->week * 604800.0 - (GEph + prn - 1)->toc;
	GRes.dv = (GEph + prn - 1)->f1 + 2 * (GEph + prn - 1)->f2 * del_t + del_trel;
}

/// <summary>
/// 计算北斗卫星位置
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="BEph">北斗星历</param>
/// <param name="BRes">传出北斗结果</param>
/// <param name="t">信号发射时刻</param>
void BDSGetSatpos(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t)
{
	double e,n0,n,Mk,Ek,Ek0,sinvk,cosvk,vk,Phik,del_uk,del_rk,del_ik,uk,rk,ik,xkp,ykp,OMGk,xgk,ygk,zgk;
	double itethres = 1.0E-10;
	double F = (-2 * sqrt(GM_BDS)) / (CLIGHT * CLIGHT);
	double tk;

	BDSTime bdt = gpst2bdst(t);
	tk = (bdt.week - (BEph + prn - 1)->week) * 604800.0 + bdt.sow - (BEph + prn - 1)->toc;
	n0 = sqrt(GM_BDS / pow((BEph + prn - 1)->A, 3));
	n = n0 + (BEph + prn - 1)->del_n;
	Mk = (BEph + prn - 1)->M0 + n * tk;
	e = (BEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}
	sinvk = (sqrt(1 - e * e) * sin(Ek) / (1 - e * cos(Ek)));
	cosvk = (cos(Ek) - e) / (1 - e * cos(Ek));
	vk = atan2(sinvk, cosvk);
	Phik = vk + (BEph + prn - 1)->omg;
	del_uk = (BEph + prn - 1)->cus * sin(2 * Phik) + (BEph + prn - 1)->cuc * cos(2 * Phik);
	del_rk = (BEph + prn - 1)->crs * sin(2 * Phik) + (BEph + prn - 1)->crc * cos(2 * Phik);
	del_ik = (BEph + prn - 1)->cis * sin(2 * Phik) + (BEph + prn - 1)->cic * cos(2 * Phik);
	uk = Phik + del_uk;
	rk = (BEph + prn - 1)->A * (1 - e * cos(Ek)) + del_rk;
	ik = (BEph + prn - 1)->I0 + del_ik + (BEph + prn - 1)->Idot * tk;
	xkp = rk * cos(uk);
	ykp = rk * sin(uk);
	//IGSO/MEO卫星
	if (prn > 5 && prn < 59 && R2D * (BEph + prn - 1)->I0 > 30)
	{
		OMGk = (BEph + prn - 1)->OMG0 + ((BEph + prn - 1)->OMGdot - OMGE_BDS) * tk - OMGE_BDS * (BEph + prn - 1)->toe;
		BRes.satpos.xyz.c1 = xkp * cos(OMGk) - ykp * cos(ik) * sin(OMGk);
		BRes.satpos.xyz.c2 = xkp * sin(OMGk) + ykp * cos(ik) * cos(OMGk);
		BRes.satpos.xyz.c3 = ykp * sin(ik);
	}
	//GEO卫星
	else
	{
		OMGk = (BEph + prn - 1)->OMG0 + (BEph + prn - 1)->OMGdot * tk - OMGE_BDS * (BEph + prn - 1)->toe;
		xgk = xkp * cos(OMGk) - ykp * cos(ik) * sin(OMGk);
		ygk = xkp * sin(OMGk) + ykp * cos(ik) * cos(OMGk);
		zgk = ykp * sin(ik);
		double Rx[9] = { 1, 0, 0,
						0, cos(-5 * D2R), sin(-5 * D2R),
						0, -sin(-5 * D2R), cos(-5 * D2R) };
		double Rz[9] = { cos(OMGE_BDS * tk), sin(OMGE_BDS * tk), 0,
						-sin(OMGE_BDS * tk), cos(OMGE_BDS * tk), 0,
						0, 0, 1 };
		double R[3] = { 0.0 };
		double RGK[3] = { xgk,
						ygk,
						zgk };
		MatrixMultiply(Rx, RGK, R, 3, 3, 1);
		double Rk[3] = { 0.0 };
		MatrixMultiply(Rz, R, Rk, 3, 3, 1);
		BRes.satpos.xyz.c1 = Rk[0];
		BRes.satpos.xyz.c2 = Rk[1];
		BRes.satpos.xyz.c3 = Rk[2];
	}
}

/// <summary>
/// 计算北斗卫星钟差
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="BEph">北斗星历</param>
/// <param name="BRes">传出北斗结果</param>
/// <param name="t">信号发射时刻</param>
void BDSGetSatclk(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t)
{
	double e, n0, n, Mk, Ek, Ek0, trel;
	double itethres = 1.0E-10;
	double F = (-2 * sqrt(GM_BDS)) / (CLIGHT * CLIGHT);
	double tk;
	double del_t;

	BDSTime bdt = gpst2bdst(t);
	tk = (bdt.week - (BEph + prn - 1)->week) * 604800.0 + bdt.sow - (BEph + prn - 1)->toc;
	n0 = sqrt(GM_BDS / pow((BEph + prn - 1)->A, 3));
	n = n0 + (BEph + prn - 1)->del_n;
	Mk = (BEph + prn - 1)->M0 + n * tk;
	e = (BEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}

	trel = F * e * sqrt((BEph + prn - 1)->A) * sin(Ek);
	del_t = (bdt.week - (BEph + prn - 1)->week) * 604800 + bdt.sow - (BEph + prn - 1)->toc;
	BRes.dts = (BEph + prn - 1)->f0 + (BEph + prn - 1)->f1 * del_t + (BEph + prn - 1)->f2 * del_t * del_t + trel;
}

/// <summary>
/// 计算北斗卫星速度和钟速
/// </summary>
/// <param name="prn">prn序列号</param>
/// <param name="BEph">北斗星历</param>
/// <param name="BRes">传出北斗结果</param>
/// <param name="t">信号发射时刻</param>
void BDSGetSatvel(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t)
{
	double e, n0, n, Mk, Ek, Ek0, sinvk, cosvk, vk, Phik, del_uk, del_rk, del_ik, uk, rk, ik, xkp, ykp, OMGk, xgk, ygk, zgk;
	double Ekdot, Phikdot, ukdot, rkdot, ikdot, OMGkdot, xkpdot, ykpdot, del_trel;
	double itethres = 1.0E-10;
	double F = (-2 * sqrt(GM_BDS)) / (CLIGHT * CLIGHT);
	double tk;
	double del_t;

	BDSTime bdt = gpst2bdst(t);
	tk = (bdt.week - (BEph + prn - 1)->week) * 604800.0 + bdt.sow - (BEph + prn - 1)->toc;
	n0 = sqrt(GM_BDS / pow((BEph + prn - 1)->A, 3));
	n = n0 + (BEph + prn - 1)->del_n;
	Mk = (BEph + prn - 1)->M0 + n * tk;
	e = (BEph + prn - 1)->ecc;
	int iterator = 0;
	for (Ek = Mk, Ek0 = 0; fabs(Ek - Ek0) > itethres;)
	{
		Ek0 = Ek;
		Ek = Mk + e * sin(Ek0);
		iterator++;

		if (iterator > 10) break;
	}
	sinvk = (sqrt(1 - e * e) * sin(Ek) / (1 - e * cos(Ek)));
	cosvk = (cos(Ek) - e) / (1 - e * cos(Ek));
	vk = atan2(sinvk, cosvk);
	Phik = vk + (BEph + prn - 1)->omg;
	del_uk = (BEph + prn - 1)->cus * sin(2 * Phik) + (BEph + prn - 1)->cuc * cos(2 * Phik);
	del_rk = (BEph + prn - 1)->crs * sin(2 * Phik) + (BEph + prn - 1)->crc * cos(2 * Phik);
	del_ik = (BEph + prn - 1)->cis * sin(2 * Phik) + (BEph + prn - 1)->cic * cos(2 * Phik);
	uk = Phik + del_uk;
	rk = (BEph + prn - 1)->A * (1 - e * cos(Ek)) + del_rk;
	ik = (BEph + prn - 1)->I0 + del_ik + (BEph + prn - 1)->Idot * tk;
	xkp = rk * cos(uk);
	ykp = rk * sin(uk);
	//IGSO/MEO卫星
	if (prn > 5 && prn < 59 && R2D * (BEph + prn - 1)->I0 > 30)
	{
		OMGk = (BEph + prn - 1)->OMG0 + ((BEph + prn - 1)->OMGdot - OMGE_BDS) * tk - OMGE_BDS * (BEph + prn - 1)->toe;
	}
	//GEO卫星
	else
	{
		OMGk = (BEph + prn - 1)->OMG0 + (BEph + prn - 1)->OMGdot * tk - OMGE_BDS * (BEph + prn - 1)->toe;
		xgk = xkp * cos(OMGk) - ykp * cos(ik) * sin(OMGk);
		ygk = xkp * sin(OMGk) + ykp * cos(ik) * cos(OMGk);
		zgk = ykp * sin(ik);
		double Rx[9] = { 1, 0, 0,
						0, cos(-5 * D2R), sin(-5 * D2R),
						0, -sin(-5 * D2R), cos(-5 * D2R) };
		double Rz[9] = { cos(OMGE_BDS * tk), sin(OMGE_BDS * tk), 0,
						-sin(OMGE_BDS * tk), cos(OMGE_BDS * tk), 0,
						0, 0, 1 };
		double R[3] = { 0.0 };
		double RGK[3] = { xgk,
						ygk,
						zgk };
		MatrixMultiply(Rx, RGK, R, 3, 3, 1);
		double Rk[3] = { 0.0 };
		MatrixMultiply(Rz, R, Rk, 3, 3, 1);
	}

	Ekdot = n / (1 - e * cos(Ek));
	Phikdot = sqrt((1 + e) / (1 - e)) * pow((cos(vk / 2) / cos(Ek / 2)), 2) * Ekdot;
	ukdot = 2 * ((BEph + prn - 1)->cus * cos(2 * Phik) - (BEph + prn - 1)->cuc * sin(2 * Phik)) * Phikdot + Phikdot;
	rkdot = (BEph + prn - 1)->A * e * sin(Ek) * Ekdot + 2 * ((BEph + prn - 1)->crs * cos(2 * Phik) - (BEph + prn - 1)->crc * sin(2 * Phik)) * Phikdot;
	ikdot = (BEph + prn - 1)->Idot + 2 * ((BEph + prn - 1)->cis * cos(2 * Phik) - (BEph + prn - 1)->cic * sin(2 * Phik)) * Phikdot;
	OMGkdot = (BEph + prn - 1)->OMGdot - OMGE_GPS;
	double Rdot[12] = { cos(OMGk), -sin(OMGk) * cos(ik), -(xkp * sin(OMGk) + ykp * cos(OMGk) * cos(ik)), ykp * sin(OMGk) * sin(ik),
						sin(OMGk), cos(OMGk) * cos(ik), (xkp * cos(OMGk) - ykp * sin(OMGk) * cos(ik)), -ykp * cos(OMGk) * sin(ik),
						0, sin(ik), 0, ykp * cos(ik) };
	xkpdot = rkdot * cos(uk) - rk * ukdot * sin(uk);
	ykpdot = rkdot * sin(uk) + rk * ukdot * cos(uk);
	double c[4] = { xkpdot,
					ykpdot,
					OMGkdot,
					ikdot };
	double rdot[3] = { 0.0 };
	MatrixMultiply(Rdot, c, rdot, 3, 4, 1);
	BRes.vx = rdot[0];
	BRes.vy = rdot[1];
	BRes.vz = rdot[2];
	del_trel = F * e * sqrt((BEph + prn - 1)->A) * cos(Ek) * Ekdot;
	del_t = (bdt.week - (BEph + prn - 1)->week) * 604800 + bdt.sow - (BEph + prn - 1)->toc;
	BRes.dv = (BEph + prn - 1)->f1 + 2 * (BEph + prn - 1)->f2 * del_t + del_trel;
}

/// <summary>
/// 获取GPS卫星状态
/// </summary>
/// <param name="obs">观测值</param>
/// <param name="GEph">GPS星历</param>
/// <param name="GRes">传出GPS结果</param>
/// <param name="t">信号发射时刻</param>
void GPSGetSatState(const Obs& obs, const GPSEph GEph[], SatRes& GRes, GPSTime t)
{
	int prn = obs.prn;
	if (prn < 0 || prn > MAXGPSSAT || !obs.Status) return;
	if (!GEph[prn - 1].Status) return;
	if (fabs(gpst2Sec(t) - GEph[prn - 1].week * 604800.0 - GEph[prn - 1].toe) > 7220) return;
	GRes.sys = GPS;
	GRes.prn = prn;

	GPSGetSatpos(prn, GEph, GRes, t);
	GPSGetSatclk(prn, GEph, GRes, t);
	GPSGetSatvel(prn, GEph, GRes, t);
	GRes.Status = true;
}

/// <summary>
/// 获取北斗卫星状态
/// </summary>
/// <param name="obs">观测值</param>
/// <param name="BEph">北斗星历</param>
/// <param name="BRes">传出北斗结果</param>
/// <param name="t">信号发射时刻</param>
void BDSGetSatState(const Obs& obs, const BDSEph BEph[], SatRes& BRes, GPSTime t)
{
	int prn = obs.prn;
	if (prn < 0 || prn > MAXBDSSAT || !obs.Status) return;
	if (!BEph[prn - 1].Status) return;
	if (fabs(bdst2Sec(gpst2bdst(t)) - BEph[prn - 1].week * 604800.0 - BEph[prn - 1].toe) > 3620) return;
	BRes.sys = GPS;
	BRes.prn = prn;

	BDSGetSatpos(prn, BEph, BRes, t);
	BDSGetSatclk(prn, BEph, BRes, t);
	BDSGetSatvel(prn, BEph, BRes, t);
	BRes.Status = true;
}

/// <summary>
/// 地球自转改正
/// </summary>
/// <param name="cECI">改正前坐标</param>
/// <param name="cECEF">待求改正后坐标</param>
/// <param name="sys">导航系统</param>
/// <param name="tau">信号传播时间</param>
void EarthRotate(const double cECI[3], double cECEF[3], NavSys sys, double tau)
{
	double OMGE;
	switch (sys)
	{
	case GPS:
		OMGE = OMGE_GPS;
		break;
	case BDS:
		OMGE = OMGE_BDS;
		break;
	default:
		OMGE = OMGE_GPS;
		break;
	}

	double cosOMGT = cos(OMGE * tau);
	double sinOMGT = sin(OMGE * tau);
	double Rotate[9] = { cosOMGT, sinOMGT, 0,
						-sinOMGT, cosOMGT, 0,
						0, 0, 1 };
	MatrixMultiply(Rotate, cECI, cECEF, 3, 3, 1);
}