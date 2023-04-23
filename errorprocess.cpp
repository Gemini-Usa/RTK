#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "function.h"

/// <summary>
/// 粗差探测
/// </summary>
/// <param name="rCurr">当前历元观测值</param>
/// <param name="rPrev">上一历元观测值</param>
void DetectOutlier(Range* rCurr, const Range* rPrev)
{
	int num = rCurr->sat_num;
	int prn, n, epc;
	double MW_k1, MW_k2, MW_k3, MW_k4, lambda1, lambda2, dMW, dGF, IF_k1, IF_k2;
	double timedur = gpst2Sec(rCurr->time) - gpst2Sec(rPrev->time);
	/* 对Range中的每个Obs进行循环,完成粗差探测并赋Status值 */
	for (int i = 0; i < num; i++)
	{
		/* 需要直接跳过循环的情况:
			当前历元没有双频数据
		*/
		bool dualfreq1 = fabs(rCurr->obs[i].P[0]) < 1.0E-10 || fabs(rCurr->obs[i].P[1]) < 1.0E-10 || fabs(rCurr->obs[i].L[0]) < 1.0E-10 || fabs(rCurr->obs[i].L[1]) < 1.0E-10;
		if (dualfreq1 == true)
		{
			rCurr->obs[i].Status = false;
			rCurr->obs[i].combination.epoch = 0;
			continue;
		}
		
		//储存线性组合观测值
		switch (rCurr->obs[i].sys)
		{
			case GPS:
			{
				MW_k1 = CLIGHT / (FREQ_L1 - FREQ_L2);
				MW_k2 = -CLIGHT / (FREQ_L1 - FREQ_L2);
				MW_k3 = -FREQ_L1 / (FREQ_L1 + FREQ_L2);
				MW_k4 = -FREQ_L2 / (FREQ_L1 + FREQ_L2);
				lambda1 = CLIGHT / FREQ_L1;
				lambda2 = CLIGHT / FREQ_L2;
				IF_k1 = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
				IF_k2 = -(FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
				break;
			}
			case BDS:
			{
				MW_k1 = CLIGHT / (FREQ_B1 - FREQ_B3);
				MW_k2 = -CLIGHT / (FREQ_B1 - FREQ_B3);
				MW_k3 = -FREQ_B1 / (FREQ_B1 + FREQ_B3);
				MW_k4 = -FREQ_B3 / (FREQ_B1 + FREQ_B3);
				lambda1 = CLIGHT / FREQ_B1;
				lambda2 = CLIGHT / FREQ_B3;
				IF_k1 = (FREQ_B1 * FREQ_B1) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
				IF_k2 = -(FREQ_B3 * FREQ_B3) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
				break;
			}
			default:
			{
				continue;
			}
		}
		rCurr->obs[i].combination.MW = MW_k1 * rCurr->obs[i].L[0] + MW_k2 * rCurr->obs[i].L[1] + MW_k3 * rCurr->obs[i].P[0] + MW_k4 * rCurr->obs[i].P[1];
		rCurr->obs[i].combination.L_GF = lambda1 * rCurr->obs[i].L[0] - lambda2 * rCurr->obs[i].L[1];
		rCurr->obs[i].combination.P_GF = rCurr->obs[i].P[0] - rCurr->obs[i].P[1];

		/* 进行粗差探测,MW平滑值计算,返回对应的Status,需要跳过循环的情况:
			1.当前历元卫星未出现在上一历元
			2.前一历元没有双频数据
			3.当前历元观测值和前一历元观测值的时间超过1s
		*/
		prn = 0;
		for (int j = 0; j < rPrev->sat_num; j++)
		{
			if (rPrev->obs[j].prn == rCurr->obs[i].prn && rPrev->obs[j].sys == rCurr->obs[i].sys)
			{
				prn = rCurr->obs[i].prn;
				n = j;
				break;
			}
		}
		if (prn == 0)
		{
			rCurr->obs[i].combination.epoch = 0;
			rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
			rCurr->obs[i].Status = false;
			continue;
		}

		if (rPrev->obs[n].combination.epoch != 0)
		{
			if (timedur > 1.1 || timedur < 0)
			{
				rCurr->obs[i].combination.epoch = 1;
				rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
				rCurr->obs[i].Status = false;
				continue;
			}

			//进行历元间差分和粗差探测
			dMW = rCurr->obs[i].combination.MW - rPrev->obs[n].combination.MW_smooth;
			dGF = rCurr->obs[i].combination.L_GF - rPrev->obs[n].combination.L_GF;
			if (fabs(dMW) > 1.5 || fabs(dGF) > 0.5)
			{
				rCurr->obs[i].combination.epoch = 0;
				rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
				rCurr->obs[i].Status = false;
				continue;
			}
			//进行MW平滑值的计算
			epc = rPrev->obs[n].combination.epoch + 1;
			rCurr->obs[i].combination.MW_smooth = (1.0 / epc) * rCurr->obs[i].combination.MW + ((epc - 1.0) / epc) * rPrev->obs[n].combination.MW_smooth;
			rCurr->obs[i].combination.epoch = epc;
			//进行IF组合的计算
			rCurr->obs[i].combination.P_IF = IF_k1 * rCurr->obs[i].P[0] + IF_k2 * rCurr->obs[i].P[1];
			rCurr->obs[i].combination.L_IF = IF_k1 * rCurr->obs[i].L[0] + IF_k2 * rCurr->obs[i].L[1];
			rCurr->obs[i].Status = true;
		}
		else
		{
			//进行MW平滑值的计算
			rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
			rCurr->obs[i].combination.epoch = 1;
			//进行IF组合的计算
			rCurr->obs[i].combination.P_IF = IF_k1 * rCurr->obs[i].P[0] + IF_k2 * rCurr->obs[i].P[1];
			rCurr->obs[i].combination.L_IF = IF_k1 * rCurr->obs[i].L[0] + IF_k2 * rCurr->obs[i].L[1];
			rCurr->obs[i].Status = false;
		}
	}
}

/// <summary>
/// 计算对流层延迟
/// </summary>
/// <param name="h">高程</param>
/// <param name="elev">高度角</param>
/// <returns>对流层延迟改正T</returns>
double Hopfield(double h, double elev)
{
	if (fabs(h) < 10.0 * 1E3)
	{
		double hd = 40136.0 + 148.72 * (T0 - 273.16);
		double hw = 11000.0;
		double T = T0 - 0.0065 * (h - H0);
		double p = p0 * pow((1.0 - 0.0000226 * (h - H0)), 5.225);
		double RH = RH0 * exp(-0.0006396 * (h - H0));
		double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
		double Kd = (155.2 * 1.0E-7 * p * (hd - h)) / T;
		double Kw = (155.2 * 1.0E-7 * 4810.0 * e * (hw - h)) / (T * T);
		double del_d = Kd / (sin(D2R * sqrt(pow(elev * R2D, 2) + 6.25 * 6.25)));
		double del_w = Kw / (sin(D2R * sqrt(pow(elev * R2D, 2) + 2.25 * 2.25)));
		return (del_d + del_w);
	}
	else return 0;
}