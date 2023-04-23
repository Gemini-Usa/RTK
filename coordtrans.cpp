#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "function.h"
#include "time.h"

/// <summary>
/// ������ת��������
/// </summary>
/// <param name="t">������</param>
/// <returns>��������</returns>
ModJulianDay JD2MJD(JulianDay t)
{
	ModJulianDay time;
	time.days = t.days - 2400000;
	time.fracday = t.fracday - 0.5;
	return time;
}

/// <summary>
/// ��������ת������
/// </summary>
/// <param name="t">��������</param>
/// <returns>������</returns>
JulianDay MJD2JD(ModJulianDay t)
{
	JulianDay time;
	time.days = t.days + 2400000;
	time.fracday = t.fracday + 0.5;
	return time;
}

/// <summary>
/// GPSʱת����ʱ
/// </summary>
/// <param name="t">GPSʱ</param>
/// <returns>����ʱ</returns>
BDSTime gpst2bdst(GPSTime t)
{
	BDSTime time;
	time.week = t.week - 1356;
	time.sow = t.sow - 14.0;
	return time;
}

/// <summary>
/// ����ʱתGPSʱ
/// </summary>
/// <param name="t">����ʱ</param>
/// <returns>GPSʱ</returns>
GPSTime bdst2gpst(BDSTime t)
{
	GPSTime time;
	time.week = t.week + 1356;
	time.sow = t.sow + 14.0;
	return time;
}

/// <summary>
/// GPSʱת��
/// </summary>
/// <param name="t">GPSʱ</param>
/// <returns>����</returns>
double gpst2Sec(GPSTime t)
{
	return (604800.0 * t.week + t.sow);
}

/// <summary>
/// ����ʱת��
/// </summary>
/// <param name="t">����ʱ</param>
/// <returns>����</returns>
double bdst2Sec(BDSTime t)
{
	return (604800.0 * t.week + t.sow);
}

/// <summary>
/// ����ʱת������
/// </summary>
/// <param name="t">����ʱ</param>
/// <returns>������</returns>
JulianDay Common2Julian(CommonDay t)
{
	JulianDay time;
	short Y = t.year;
	unsigned short M = t.month;
	unsigned short D = t.date;
	double UT = (double)(t.hour) + (double)(t.minute) / 60.0 + ((double)(t.second) + t.fracsec) / 3600.0;
	double jd;

	if (M <= 2)
	{
		Y -= 1;
		M += 12;
	}
	jd = (double)((int)(365.25 * (double)(Y)) + (int)(30.6001 * (double)(M + 1)) + D + UT / 24.0 + 1720981.5);
	time.days = (int)jd;
	time.fracday = jd - time.days;
	return time;
}

/// <summary>
/// ������ת����ʱ
/// </summary>
/// <param name="t">������</param>
/// <returns>����ʱ</returns>
CommonDay Julian2Common(JulianDay t)
{
	double hh, mm, ss = 0.0;
	CommonDay time;
	int a = t.days + (int)(t.fracday + 0.5);
	int b = a + 1537;
	int c = (int)((b - 122.1) / 365.25);
	int d = (int)(365.25 * c);
	int e = (int)((double)(b - d) / 30.6001);
	double days = (double)(b - d - (int)(30.6001 * e)) + t.fracday + 0.5;
	time.date = (int)days;
	time.month = e - 1 - 12 * (int)(e / 14);
	time.year = c - 4715 - (int)((7 + time.month) / 10);
	hh = (t.fracday + 0.5 - (int)(t.fracday + 0.5)) * 24;
	mm = (hh - (int)hh) * 60;
	ss = (mm - (int)mm) * 60;
	time.hour = (int)hh;
	time.minute = (int)mm;
	time.second = (int)ss;
	time.fracsec = ss - time.second;
	return time;
}

/// <summary>
/// ��������תGPSʱ
/// </summary>
/// <param name="t">��������</param>
/// <returns>GPSʱ</returns>
GPSTime MJD2gpst(ModJulianDay t)
{
	GPSTime time;
	time.week = (int)(((double)(t.days - 44244) + t.fracday) / 7.0);
	time.sow = ((double)(t.days - 44244) + t.fracday - (double)time.week * 7.0) * 86400.0;
	return time;
}

/// <summary>
/// GPSʱת��������
/// </summary>
/// <param name="t">GPSʱ</param>
/// <returns>��������</returns>
ModJulianDay gpst2MJD(GPSTime t)
{
	ModJulianDay time;
	double mjd = 44244.0 + (double)(t.week) * 7.0 + t.sow / 86400.0;
	time.days = (int)mjd;
	time.fracday = mjd - time.days;
	return time;
}

/// <summary>
/// XYZ����תBLH����
/// </summary>
/// <param name="c">XYZ����</param>
/// <param name="EllipseType">�ο�����</param>
/// <returns>BLH����</returns>
BLH XYZ2BLH(XYZ c, int EllipseType)
{
	BLH coord;
	double RE, FE = 0.0;
	if (fabs(c.xyz.c1) < 1E-10 && fabs(c.xyz.c2) < 1E-10 && fabs(c.xyz.c3) < 1E-10) return coord;

	switch (EllipseType)
	{
	case WGS84:
		RE = RE_WGS84;
		FE = FE_WGS84;
		break;
	case CGCS2000:
		RE = RE_CGCS2000;
		FE = FE_CGCS2000;
		break;
	default:
		RE = 0.0;
		FE = 0.0;
		break;
	}
	double e2 = 2.0 * FE - FE * FE;//��һƫ���ʵ�ƽ��
	double threshold = 1E-10;//������ֵ
	double r = sqrt(pow(c.xyz.c1, 2) + pow(c.xyz.c2, 2) + pow(c.xyz.c3, 2));
	double phi = atan2(c.xyz.c3, sqrt(c.xyz.c1 * c.xyz.c1 + c.xyz.c2 * c.xyz.c2));
	double b, bk, W = 0.0;
	int iterator = 0;
	for (b = phi, bk = 0.0;fabs(b - bk) > threshold;)
	{
		bk = b;
		W = sqrt(1 - e2 * sin(b) * sin(b));
		b = atan2((c.xyz.c3 * W + RE * e2 * sin(b)) * tan(phi), c.xyz.c3 * W);
		iterator++;

		if (iterator > 10) break;
	}
	coord.blh.c1 = r > 1E-12 ? b : (c.xyz.c3 > 0 ? PI / 2 : -PI / 2);
	coord.blh.c2 = r > 1E-12 ? atan2(c.xyz.c2, c.xyz.c1) : 0.0;
	coord.blh.c3 = r > 1E-12 ? (r * cos(phi)) / cos(b) - RE / W : 0;
	return coord;
}

/// <summary>
/// BLH����תXYZ����
/// </summary>
/// <param name="c">BLH����</param>
/// <param name="EllipseType">�ο�����</param>
/// <returns>XYZ����</returns>
XYZ BLH2XYZ(BLH c, int EllipseType)
{
	XYZ coord;
	double RE, FE = 0.0;
	switch (EllipseType)
	{
	case WGS84:
		RE = RE_WGS84;
		FE = FE_WGS84;
		break;
	case CGCS2000:
		RE = RE_CGCS2000;
		FE = FE_CGCS2000;
		break;
	default:
		RE = 0.0;
		FE = 0.0;
		break;
	}
	double e2 = 2 * FE - FE * FE;//��һƫ���ʵ�ƽ��
	double radius = RE / (sqrt(1 - e2 * pow(sin(c.blh.c1), 2)));//î��Ȧ�뾶
	coord.xyz.c1 = (radius + c.blh.c3) * cos(c.blh.c1) * cos(c.blh.c2);
	coord.xyz.c2 = (radius + c.blh.c3) * cos(c.blh.c1) * sin(c.blh.c2);
	coord.xyz.c3 = (radius * (1 - e2) + c.blh.c3) * sin(c.blh.c1);
	return coord;
}

/// <summary>
/// XYZ����תվ������
/// </summary>
/// <param name="refpoint">վ��XYZ����</param>
/// <param name="point">��תXYZ����</param>
/// <param name="EllipseType">�ο�����</param>
/// <returns>վ������</returns>
ENU XYZ2ENU(XYZ refpoint, XYZ point, int EllipseType)
{
	ENU coord;
	BLH refblh = XYZ2BLH(refpoint, EllipseType);
	double sinl = sin(refblh.blh.c2), cosl = cos(refblh.blh.c2), sinb = sin(refblh.blh.c1), cosb = cos(refblh.blh.c1);
	double E[9] = { -sinl, cosl, 0,
					-sinb * cosl, -sinb * sinl, cosb,
					 cosb * cosl,  cosb * sinl, sinb };
	double del_r[3] = { point.xyz.c1 - refpoint.xyz.c1,
						point.xyz.c2 - refpoint.xyz.c2,
						point.xyz.c3 - refpoint.xyz.c3 };
	double enu[3] = { 0 };
	MatrixMultiply(E, del_r, enu, 3, 3, 1);
	coord.enu.c1 = enu[0];
	coord.enu.c2 = enu[1];
	coord.enu.c3 = enu[2];
	return coord;
}

void GetENUTransMat(const XYZ& refpoint, int EllipseType, double Mat[9])
{
	BLH refblh = XYZ2BLH(refpoint, EllipseType);
	double sinl = sin(refblh.blh.c2), cosl = cos(refblh.blh.c2), sinb = sin(refblh.blh.c1), cosb = cos(refblh.blh.c1);
	double E[9] = { -sinl, cosl, 0,
					-sinb * cosl, -sinb * sinl, cosb,
					 cosb * cosl,  cosb * sinl, sinb };
	for (int i = 0; i < 9; i++)
	{
		Mat[i] = E[i];
	}
}

/// <summary>
/// ����߶Ƚ�
/// </summary>
/// <param name="refpoint">վ��XYZ����</param>
/// <param name="point">����XYZ����</param>
/// <param name="EllipseType">�ο�����</param>
/// <returns>�߶Ƚ�</returns>
double GetElev(XYZ refpoint, XYZ point, int EllipseType)
{
	if (fabs(refpoint.xyz.c1) < 1E-10 && fabs(refpoint.xyz.c2) < 1E-10 && fabs(refpoint.xyz.c3) < 1E-10) return 90 * D2R;
	ENU coord = XYZ2ENU(refpoint, point, EllipseType);

	double sinE = coord.enu.c3 / Norm(coord.denu, 3);
	return asin(sinE);
}

/// <summary>
/// ����XYZ��������
/// </summary>
/// <param name="R1">XYZ����1</param>
/// <param name="R2">XYZ����2</param>
/// <returns>����</returns>
double GetDist(XYZ R1, XYZ R2)
{
	double dR[3] = { R1.dxyz[0] - R2.dxyz[0], R1.dxyz[1] - R2.dxyz[1], R1.dxyz[2] - R2.dxyz[2] };
	return Norm(dR, 3);
}

double GetAzim(const XYZ& refpoint, const XYZ& point, int EllipseType)
{
	ENU coord = XYZ2ENU(refpoint, point, EllipseType);
	return atan2(coord.enu.c1, coord.enu.c2);
}

void GetLOSVector(const XYZ& sat_pos, const XYZ& rcv_pos, double* los)
{
    double dr[3];
    VectorSub(sat_pos.dxyz, rcv_pos.dxyz, dr, 3);
    VectorNormalize(dr, los, 3);
}