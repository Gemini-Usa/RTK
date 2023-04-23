#ifndef constants_H
#define constants_H

#define PI (4.0 * atan(1.0))		//Բ����
#define R2D (180.0 / PI)			//����ת�Ƕ�
#define D2R (PI / 180.0)			//�Ƕ�ת����
#define CLIGHT (2.99792458 * 1.0E8) //����

#define WGS84 1						//WGS84����
#define CGCS2000 2					//CGCS2000����
#define RE_WGS84 6378137.0			//WGS84���򳤰���
#define FE_WGS84 (1.0/298.257223563)//WGS84�������
#define RE_CGCS2000 6378137.0		//CGCS2000���򳤰���
#define FE_CGCS2000 (1.0/298.257222101) //CGCS2000�������

#define MAXNOVDLEN 40480			//����������
#define MAXFREQ 2					//Ƶ����
#define POLYCRC32 0xEDB88320u		//32λCRCУ����
#define MAXGPSSAT 32				//GPS���32������
#define MAXBDSSAT 63				//BDS���63������
#define MAXSATNUM (32 + 63)			//���������, GPS���32��, BDS���63��

#define GM_GPS (3.986005 * 1.0E14)	//GPS������������
#define GM_BDS (3.986004418 * 1.0E14)	//BDS������������
#define OMGE_GPS (7.2921151467 * 1.0E-5)//GPS������ת���ٶ�
#define OMGE_BDS (7.2921150 * 1.0E-5)	//BDS������ת���ٶ�
#define FREQ_L1 1.57542E9			//GPS L1Ƶ��(Hz)
#define FREQ_L2 1.22760E9			//GPS L2Ƶ��(Hz)
#define FREQ_L5 1.17645E9           //GPS L5Ƶ��(Hz)
#define FREQ_B1 1.561098E9			//BDS B1Ƶ��(Hz)
#define FREQ_B2 1.20714E9           //BDS B2Ƶ��(Hz)
#define FREQ_B3 1.26852E9           //BDS B3Ƶ��(Hz)

#define H0 0.0						//��ƽ��߶�
#define T0 (15 + 273.16)			//��ƽ���¶�
#define p0 (1013.25)				//��ƽ����ѹ
#define RH0 (0.5)					//���ʪ��

constexpr double maxtimediff = 0.1; // ��վ������վ���������ʱ���

#endif