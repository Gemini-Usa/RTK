#include <stdio.h>
#include <iostream>
#include <string.h>
#include "constants.h"
#include "function.h"

/// <summary>
/// 无符号字符解码
/// </summary>
/// <typeparam name="T">待转数据类型</typeparam>
/// <param name="buff">字符指针</param>
/// <returns>解码结果</returns>
template <typename T>
T DecodeUchar(unsigned char* buff)
{
	T res = 0;
	int bytes = sizeof(T);
	memcpy(&res, buff, bytes);

	return res;
}

/// <summary>
/// crc校验
/// </summary>
/// <param name="buff">字符指针</param>
/// <param name="len">长度</param>
/// <returns>crc结果</returns>
unsigned int crc32(const unsigned char* buff, int len)
{
	int i, j;
	unsigned int crc = 0;
	for (i = 0; i < len; i++)
	{
		crc ^= buff[i];
		for (j = 0; j < 8; j++)
		{
			if (crc & 1) crc = (crc >> 1) ^ POLYCRC32;
			else crc >>= 1;
		}
	}
	return crc;
}

/// <summary>
/// 观测值解码
/// </summary>
/// <param name="buff">字符指针</param>
/// <param name="range">传出观测值</param>
/// <returns>true:含有观测值 false:不含观测值</returns>
bool DecodeRange(unsigned char* buff, Range* range)
{
	int n, i, j, f, prn;
	unsigned char* p = buff + 28;
	unsigned int track;
	NavSys sys = NavSys::UNK;
	bool Status = false;

	range->sat_num = DecodeUchar<unsigned int>(p);
	if (range->sat_num <= 0) return false;
	memset(range->obs, 0, MAXFREQ * sizeof(MAXSATNUM));
	n = 0;
	for (p += 4, j = 0; j < range->sat_num; j++, p += 44)
	{
		//卫星系统
		track = DecodeUchar<unsigned int>(p + 40);
		auto a = (track >> 16) & 7;
		switch (a)
		{
			case 0:
				sys = GPS;
				break;
			case 4:
				sys = BDS;
				break;
			default:
				sys = NavSys::UNK;
				break;
		}

		if (sys == NavSys::UNK) continue;

		//信号类型
		if (sys == GPS)
		{
			switch ((track >> 21) & 0x1F)
			{
				case 0:f = 0; break;
				case 9:f = 1; break;
				default:f = 2; break;
			}
		}
		if (sys == BDS)
		{
			switch ((track >> 21) & 0x1F)
			{
				case 0:f = 0; break;
				case 2:f = 1; break;
				case 4:f = 0; break;
				case 6:f = 1; break;
				default:f = 2; break;
			}
		}
		if (f == 2) continue;

		//观测值赋值
		prn = DecodeUchar<unsigned short>(p);
		for (i = 0; i < n; i++)
		{
			if (range->obs[i].sys == sys && range->obs[i].prn == prn)
			{
				n = i;
				break;
			}
		}
		range->obs[n].prn = prn;
		range->obs[n].sys = sys;
        range->obs[n].parity[f] = static_cast<short>((track >> 11) & 1);
		range->obs[n].P[f] = DecodeUchar<double>(p + 4);
		range->obs[n].Pstd[f] = DecodeUchar<float>(p + 12);
		range->obs[n].L[f] = -DecodeUchar<double>(p + 16);
		range->obs[n].Lstd[f] = DecodeUchar<float>(p + 24);
		range->obs[n].D[f] = -DecodeUchar<float>(p + 28);
		range->obs[n].SNR[f] = DecodeUchar<float>(p + 32);
		n++;
	}
	range->sat_num = n;
	if (range->sat_num >= 0) Status = true;
	return Status;
}

/// <summary>
/// GPS星历解码
/// </summary>
/// <param name="buff">字符指针</param>
/// <param name="GEph">传出GPS星历指针</param>
/// <returns>true:有GPS星历 false:无GPS星历</returns>
bool DecodeGPSEph(unsigned char buff[], GPSEph* GEph)
{
	int prn;
	unsigned char* p = buff + 28;
	prn = DecodeUchar<unsigned int>(p);
	if (prn < 1 || prn > 32) return false;

	(GEph + prn - 1)->PRN = DecodeUchar<unsigned int>(p);
	(GEph + prn - 1)->tow = DecodeUchar<double>(p + 4);
	(GEph + prn - 1)->Health = DecodeUchar<unsigned int>(p + 12);
	(GEph + prn - 1)->IODE[0] = DecodeUchar<unsigned int>(p + 16);
	(GEph + prn - 1)->IODE[1] = DecodeUchar<unsigned int>(p + 20);
	(GEph + prn - 1)->week = DecodeUchar<unsigned int>(p + 24);
	(GEph + prn - 1)->z_week = DecodeUchar<unsigned int>(p + 28);
	(GEph + prn - 1)->toe = DecodeUchar<double>(p + 32);
	(GEph + prn - 1)->A = DecodeUchar<double>(p + 40);
	(GEph + prn - 1)->del_n = DecodeUchar<double>(p + 48);
	(GEph + prn - 1)->M0 = DecodeUchar<double>(p + 56);
	(GEph + prn - 1)->ecc = DecodeUchar<double>(p + 64);
	(GEph + prn - 1)->omg = DecodeUchar<double>(p + 72);
	(GEph + prn - 1)->cuc = DecodeUchar<double>(p + 80);
	(GEph + prn - 1)->cus = DecodeUchar<double>(p + 88);
	(GEph + prn - 1)->crc = DecodeUchar<double>(p + 96);
	(GEph + prn - 1)->crs = DecodeUchar<double>(p + 104);
	(GEph + prn - 1)->cic = DecodeUchar<double>(p + 112);
	(GEph + prn - 1)->cis = DecodeUchar<double>(p + 120);
	(GEph + prn - 1)->I0 = DecodeUchar<double>(p + 128);
	(GEph + prn - 1)->Idot = DecodeUchar<double>(p + 136);
	(GEph + prn - 1)->OMG0 = DecodeUchar<double>(p + 144);
	(GEph + prn - 1)->OMGdot = DecodeUchar<double>(p + 152);
	(GEph + prn - 1)->iodc = DecodeUchar<unsigned int>(p + 160);
	(GEph + prn - 1)->toc = DecodeUchar<double>(p + 164);
	(GEph + prn - 1)->tgd = DecodeUchar<double>(p + 172);
	(GEph + prn - 1)->f0 = DecodeUchar<double>(p + 180);
	(GEph + prn - 1)->f1 = DecodeUchar<double>(p + 188);
	(GEph + prn - 1)->f2 = DecodeUchar<double>(p + 196);
	(GEph + prn - 1)->AS = (bool)DecodeUchar<unsigned int>(p + 204);
	(GEph + prn - 1)->URA = DecodeUchar<double>(p + 216);
	(GEph + prn - 1)->Status = true;
	return true;
}


/// <summary>
/// 北斗星历解码
/// </summary>
/// <param name="buff">字符指针</param>
/// <param name="CEph">传出北斗星历指针</param>
/// <returns>true:有北斗星历 false:无北斗星历</returns>
bool DecodeBDSEph(unsigned char buff[], BDSEph* CEph)
{
	int prn;
	unsigned char* p = buff + 28;
	prn = DecodeUchar<unsigned int>(p);
	if (prn < 1 || prn > 63) return false;

	(CEph + prn - 1)->PRN = DecodeUchar<unsigned int>(p);
	(CEph + prn - 1)->week = DecodeUchar<unsigned int>(p + 4);
	(CEph + prn - 1)->URA = DecodeUchar<double>(p + 8);
	(CEph + prn - 1)->Health = DecodeUchar<unsigned int>(p + 16);
	(CEph + prn - 1)->tgd[0] = DecodeUchar<double>(p + 20);
	(CEph + prn - 1)->tgd[1] = DecodeUchar<double>(p + 28);
	(CEph + prn - 1)->AODC = DecodeUchar<unsigned int>(p + 36);
	(CEph + prn - 1)->toc = DecodeUchar<unsigned int>(p + 40);
	(CEph + prn - 1)->f0 = DecodeUchar<double>(p + 44);
	(CEph + prn - 1)->f1 = DecodeUchar<double>(p + 52);
	(CEph + prn - 1)->f2 = DecodeUchar<double>(p + 60);
	(CEph + prn - 1)->AODE = DecodeUchar<unsigned int>(p + 68);
	(CEph + prn - 1)->toe = DecodeUchar<unsigned int>(p + 72);
	(CEph + prn - 1)->A = pow(DecodeUchar<double>(p + 76), 2);
	(CEph + prn - 1)->ecc = DecodeUchar<double>(p + 84);
	(CEph + prn - 1)->omg = DecodeUchar<double>(p + 92);
	(CEph + prn - 1)->del_n = DecodeUchar<double>(p + 100);
	(CEph + prn - 1)->M0 = DecodeUchar<double>(p + 108);
	(CEph + prn - 1)->OMG0 = DecodeUchar<double>(p + 116);
	(CEph + prn - 1)->OMGdot = DecodeUchar<double>(p + 124);
	(CEph + prn - 1)->I0 = DecodeUchar<double>(p + 132);
	(CEph + prn - 1)->Idot = DecodeUchar<double>(p + 140);
	(CEph + prn - 1)->cuc = DecodeUchar<double>(p + 148);
	(CEph + prn - 1)->cus = DecodeUchar<double>(p + 156);
	(CEph + prn - 1)->crc = DecodeUchar<double>(p + 164);
	(CEph + prn - 1)->crs = DecodeUchar<double>(p + 172);
	(CEph + prn - 1)->cic = DecodeUchar<double>(p + 180);
	(CEph + prn - 1)->cis = DecodeUchar<double>(p + 188);
	(CEph + prn - 1)->Status = true;
	return true;
}

/// <summary>
/// 解码OEM7格式数据
/// </summary>
/// <param name="buff">二进制数据</param>
/// <param name="LenRem">传出数组剩余长度</param>
/// <param name="range">传出观测值指针</param>
/// <param name="GEph">传出GPS星历指针</param>
/// <param name="CEph">传出北斗星历指针</param>
/// <returns>true:含有观测值 false:不含观测值</returns>
bool DecodeNovOem7(unsigned char buff[], int &LenRem, Range *range, GPSEph *GEph, BDSEph *CEph, int mode)
{
	int i, MsgLen, MsgID;
	GPSTime gt;
	bool Status = false;
    int byte_num = mode == 1 ? LenRem : MAXNOVDLEN;

	for (i = 0; i < byte_num - 2; i++)
	{
		//同步字
		if ((buff[i] == 0xAA && buff[i + 1] == 0x44 && buff[i + 2] == 0x12) == false) continue;
		//文件头
		if ((i + 28) >= byte_num) break;
		MsgID = DecodeUchar<unsigned short>(buff + i + 4);
		MsgLen = DecodeUchar<unsigned short>(buff + i + 8) + 32;
		gt.week = DecodeUchar<unsigned short>(buff + i + 14);
		gt.sow = DecodeUchar<unsigned int>(buff + i + 16) * 1.0E-3;
		if ((i + MsgLen) >= byte_num) break;
		//CRC校验
		if (crc32(buff + i, MsgLen - 4) != DecodeUchar<unsigned int>(buff + i + MsgLen - 4))
		{
			i += MsgLen;
			continue;
		}
		range->time.week = gt.week;
		range->time.sow = gt.sow;
		//数据解码
		switch (MsgID)
		{
			case 43:
				Status = DecodeRange(buff + i, range);
				LenRem = 0;
				for (i = i + MsgLen; i < byte_num; i++) {
					buff[LenRem] = buff[i];
					LenRem++;
				}
				return Status;
			case 7:
				DecodeGPSEph(buff + i, GEph);
				break;
			case 1696:
				DecodeBDSEph(buff + i, CEph);
				break;
			default:
				break;
		}
		i += MsgLen - 1;
	}
	LenRem = 0;
	for (; i < byte_num; i++) {
		buff[LenRem] = buff[i];
		LenRem++;
	}

	return Status;
}