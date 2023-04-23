#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "function.h"

/// <summary>
/// �ֲ�̽��
/// </summary>
/// <param name="rCurr">��ǰ��Ԫ�۲�ֵ</param>
/// <param name="rPrev">��һ��Ԫ�۲�ֵ</param>
void DetectOutlier(Range* rCurr, const Range* rPrev)
{
	int num = rCurr->sat_num;
	int prn, n, epc;
	double MW_k1, MW_k2, MW_k3, MW_k4, lambda1, lambda2, dMW, dGF, IF_k1, IF_k2;
	double timedur = gpst2Sec(rCurr->time) - gpst2Sec(rPrev->time);
	/* ��Range�е�ÿ��Obs����ѭ��,��ɴֲ�̽�Ⲣ��Statusֵ */
	for (int i = 0; i < num; i++)
	{
		/* ��Ҫֱ������ѭ�������:
			��ǰ��Ԫû��˫Ƶ����
		*/
		bool dualfreq1 = fabs(rCurr->obs[i].P[0]) < 1.0E-10 || fabs(rCurr->obs[i].P[1]) < 1.0E-10 || fabs(rCurr->obs[i].L[0]) < 1.0E-10 || fabs(rCurr->obs[i].L[1]) < 1.0E-10;
		if (dualfreq1 == true)
		{
			rCurr->obs[i].Status = false;
			rCurr->obs[i].combination.epoch = 0;
			continue;
		}
		
		//����������Ϲ۲�ֵ
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

		/* ���дֲ�̽��,MWƽ��ֵ����,���ض�Ӧ��Status,��Ҫ����ѭ�������:
			1.��ǰ��Ԫ����δ��������һ��Ԫ
			2.ǰһ��Ԫû��˫Ƶ����
			3.��ǰ��Ԫ�۲�ֵ��ǰһ��Ԫ�۲�ֵ��ʱ�䳬��1s
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

			//������Ԫ���ֺʹֲ�̽��
			dMW = rCurr->obs[i].combination.MW - rPrev->obs[n].combination.MW_smooth;
			dGF = rCurr->obs[i].combination.L_GF - rPrev->obs[n].combination.L_GF;
			if (fabs(dMW) > 1.5 || fabs(dGF) > 0.5)
			{
				rCurr->obs[i].combination.epoch = 0;
				rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
				rCurr->obs[i].Status = false;
				continue;
			}
			//����MWƽ��ֵ�ļ���
			epc = rPrev->obs[n].combination.epoch + 1;
			rCurr->obs[i].combination.MW_smooth = (1.0 / epc) * rCurr->obs[i].combination.MW + ((epc - 1.0) / epc) * rPrev->obs[n].combination.MW_smooth;
			rCurr->obs[i].combination.epoch = epc;
			//����IF��ϵļ���
			rCurr->obs[i].combination.P_IF = IF_k1 * rCurr->obs[i].P[0] + IF_k2 * rCurr->obs[i].P[1];
			rCurr->obs[i].combination.L_IF = IF_k1 * rCurr->obs[i].L[0] + IF_k2 * rCurr->obs[i].L[1];
			rCurr->obs[i].Status = true;
		}
		else
		{
			//����MWƽ��ֵ�ļ���
			rCurr->obs[i].combination.MW_smooth = rCurr->obs[i].combination.MW;
			rCurr->obs[i].combination.epoch = 1;
			//����IF��ϵļ���
			rCurr->obs[i].combination.P_IF = IF_k1 * rCurr->obs[i].P[0] + IF_k2 * rCurr->obs[i].P[1];
			rCurr->obs[i].combination.L_IF = IF_k1 * rCurr->obs[i].L[0] + IF_k2 * rCurr->obs[i].L[1];
			rCurr->obs[i].Status = false;
		}
	}
}

/// <summary>
/// ����������ӳ�
/// </summary>
/// <param name="h">�߳�</param>
/// <param name="elev">�߶Ƚ�</param>
/// <returns>�������ӳٸ���T</returns>
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