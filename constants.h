#ifndef constants_H
#define constants_H

#define PI (4.0 * atan(1.0))		//圆周率
#define R2D (180.0 / PI)			//弧度转角度
#define D2R (PI / 180.0)			//角度转弧度
#define CLIGHT (2.99792458 * 1.0E8) //光速

#define WGS84 1						//WGS84椭球
#define CGCS2000 2					//CGCS2000椭球
#define RE_WGS84 6378137.0			//WGS84地球长半轴
#define FE_WGS84 (1.0/298.257223563)//WGS84地球扁率
#define RE_CGCS2000 6378137.0		//CGCS2000地球长半轴
#define FE_CGCS2000 (1.0/298.257222101) //CGCS2000地球扁率

#define MAXNOVDLEN 40480			//缓冲区长度
#define MAXFREQ 2					//频点数
#define POLYCRC32 0xEDB88320u		//32位CRC校验码
#define MAXGPSSAT 32				//GPS最多32颗卫星
#define MAXBDSSAT 63				//BDS最多63颗卫星
#define MAXSATNUM (32 + 63)			//最多卫星数, GPS最多32颗, BDS最多63颗

#define GM_GPS (3.986005 * 1.0E14)	//GPS万有引力常量
#define GM_BDS (3.986004418 * 1.0E14)	//BDS万有引力常量
#define OMGE_GPS (7.2921151467 * 1.0E-5)//GPS地球自转角速度
#define OMGE_BDS (7.2921150 * 1.0E-5)	//BDS地球自转角速度
#define FREQ_L1 1.57542E9			//GPS L1频率(Hz)
#define FREQ_L2 1.22760E9			//GPS L2频率(Hz)
#define FREQ_L5 1.17645E9           //GPS L5频率(Hz)
#define FREQ_B1 1.561098E9			//BDS B1频率(Hz)
#define FREQ_B2 1.20714E9           //BDS B2频率(Hz)
#define FREQ_B3 1.26852E9           //BDS B3频率(Hz)

#define H0 0.0						//海平面高度
#define T0 (15 + 273.16)			//海平面温度
#define p0 (1013.25)				//海平面气压
#define RH0 (0.5)					//相对湿度

constexpr double maxtimediff = 0.1; // 基站和流动站间最大容许时间差

#endif