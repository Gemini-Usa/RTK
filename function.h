#ifndef Function_H
#define Function_H

#include "constants.h"
#include "time.h"
#include "function.h"
#include <vector>
#include <iostream>
#include <queue>

/*-----------------------------enumerate-------------------------------*/
enum NavSys
{
	GPS = 0, BDS, GLO, GAL, QZS, SBS, UNK
};
enum class Sync
{
    ROV = 0, BAS, SYN, UNK
};
/*--------------------------struct---------------------------*/
struct CommonDay
{
	short year;
	unsigned short month;
	unsigned short date;
	unsigned short hour;
	unsigned short minute;
	unsigned short second;
	double fracsec;

	CommonDay()
	{
		year = 0;
		month = date = hour = minute = second = 0;
		fracsec = 0.0;
	}
};

struct JulianDay
{
	int days;
	double fracday;

	JulianDay() { days = 0; fracday = 0.0; }
};

struct ModJulianDay
{
	int days;
	double fracday;

	ModJulianDay() { days = 0; fracday = 0.0; }
};

struct GPSTime
{
	int week;
	double sow;

	GPSTime() { week = 0; sow = 0.0; }
};

struct BDSTime
{
	int week;
	double sow;

	BDSTime() { week = 0; sow = 0.0; }
};

struct Coord
{
	double c1;
	double c2;
	double c3;

	Coord() { c1 = 0.0; c2 = 0.0; c3 = 0.0; }
};

union XYZ
{
	Coord xyz;
	double dxyz[3];

	XYZ() { xyz = Coord(); }
};

union BLH
{
	Coord blh;
	double dblh[3];

	BLH() { blh = Coord(); }
};

union ENU
{
	Coord enu;
	double denu[3];

	ENU() { enu = Coord(); }
};

struct LinearCom
{
	double MW;					//MW组和观测值
	double P_GF;				//伪距电离层残差组合观测值
	double L_GF;				//相位电离层残差组合观测值
	double P_IF;				//伪距无电离层组和观测值
	double L_IF;				//相位无电离层组合观测值
	double MW_smooth;			//MW组合平滑值
	unsigned int epoch;			//当前进行MW平滑中的历元序数

	LinearCom()
	{
		MW = 0.0; P_GF = 0.0; L_GF = 0.0; P_IF = 0.0; L_IF = 0.0; MW_smooth = 0.0; 
		epoch = 0;
	}
};

struct Obs
{
	unsigned short prn;			//卫星序列号
	NavSys sys;				    //导航系统
    short parity[MAXFREQ];      //parity known flag
	double P[MAXFREQ];			//伪距观测值,pseudorange
	double L[MAXFREQ];			//载波相位观测值,carrier-phase
	float D[MAXFREQ];			//多普勒观测值,doppler
	float SNR[MAXFREQ];			//载噪比
	float Pstd[MAXFREQ];		//伪距精度
	float Lstd[MAXFREQ];		//载波精度
	LinearCom combination;		//观测值线性组合
	bool Status;				//粗差检测通过否,不含粗差=true,含粗差=false

	Obs()
	{
        prn = 0;
        sys = NavSys::UNK;
		for (int i = 0; i < MAXFREQ; i++)
		{
            parity[i] = 0;
			P[i] = L[i] = D[i] = SNR[i] = Pstd[i] = Lstd[i] = 0;
		}
		combination = LinearCom();
		Status = false;
	}
    [[nodiscard]] bool isDoubleFreq() const {
        if (fabs(P[0]) < 1.0E-10 || fabs(P[1]) < 1.0E-10 || fabs(L[0]) < 1.0E-10 || fabs(L[1]) < 1.0E-10) {
            return false;
        }
        else return true;
    }
    [[nodiscard]] double formMW() const {
        double MW_k1, MW_k2, MW_k3, MW_k4;
        switch (sys)
        {
            case GPS:
            {
                MW_k1 = CLIGHT / (FREQ_L1 - FREQ_L2);
                MW_k2 = -CLIGHT / (FREQ_L1 - FREQ_L2);
                MW_k3 = -FREQ_L1 / (FREQ_L1 + FREQ_L2);
                MW_k4 = -FREQ_L2 / (FREQ_L1 + FREQ_L2);
                break;
            }
            case BDS:
            {
                MW_k1 = CLIGHT / (FREQ_B1 - FREQ_B3);
                MW_k2 = -CLIGHT / (FREQ_B1 - FREQ_B3);
                MW_k3 = -FREQ_B1 / (FREQ_B1 + FREQ_B3);
                MW_k4 = -FREQ_B3 / (FREQ_B1 + FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return MW_k1 * L[0] + MW_k2 * L[1] + MW_k3 * P[0] + MW_k4 * P[1];
    }
    [[nodiscard]] double formPGF() const {
        if (sys == NavSys::UNK) return 0.0;
        else return P[0] - P[1];
    }
    [[nodiscard]] double formLGF() const {
        double lambda1, lambda2;
        switch (sys)
        {
            case GPS:
            {
                lambda1 = CLIGHT / FREQ_L1;
                lambda2 = CLIGHT / FREQ_L2;
                break;
            }
            case BDS:
            {
                lambda1 = CLIGHT / FREQ_B1;
                lambda2 = CLIGHT / FREQ_B3;
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return lambda1 * L[0] - lambda2 * L[1];
    }
    [[nodiscard]] double formPIF() const {
        double IF_k1, IF_k2;
        switch (sys)
        {
            case GPS:
            {
                IF_k1 = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                IF_k2 = -(FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                break;
            }
            case BDS:
            {
                IF_k1 = (FREQ_B1 * FREQ_B1) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                IF_k2 = -(FREQ_B3 * FREQ_B3) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return IF_k1 * P[0] + IF_k2 * P[1];
    }
    [[nodiscard]] double formLIF() const {
        double IF_k1, IF_k2;
        switch (sys)
        {
            case GPS:
            {
                IF_k1 = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                IF_k2 = -(FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                break;
            }
            case BDS:
            {
                IF_k1 = (FREQ_B1 * FREQ_B1) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                IF_k2 = -(FREQ_B3 * FREQ_B3) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return IF_k1 * L[0] + IF_k2 * L[1];
    }
};

struct Range
{
	GPSTime time;				//当前历元
	Obs obs[MAXSATNUM];			//观测值
	unsigned int sat_num;		//卫星数量

	Range() 
	{ 
		time = GPSTime();
		for (int i = 0; i < MAXSATNUM; i++)
		{
			obs[i] = Obs();
		}
        sat_num = 0;
	}
};

struct GPSEph
{
	unsigned long PRN;										//卫星PRN号
	double tow;												//子帧1的时间戳
	unsigned long Health;									//卫星健康状况
	unsigned long IODE[2];									//星历, 数据的信息
	unsigned long week;										//由z计数推得的toe周
	unsigned long z_week;									//z计数
	double toe;												//卫星星历参考时间
	double A, del_n, M0, ecc, omg, I0, Idot, OMG0, OMGdot;	//开普勒轨道根数
	double cuc, cus, crc, crs, cic, cis;					//摄动改正项
	unsigned long iodc;										//数据钟信息
	double toc;												//卫星钟差改正项
	double tgd;												//群延迟参数
	double f0, f1, f2;										//钟偏, 钟速, 钟漂
	bool AS;												//防欺骗标识
	double URA;												//用户测量精度
	bool Status;											//星历是否可用

	GPSEph()
	{
		PRN = 0;
		tow = toe = toc = 0.0;
		Health = 0;
		IODE[0] = IODE[1] = 0;
		week = z_week = 0;
		A = del_n = M0 = ecc = omg = I0 = Idot = OMG0 = OMGdot = 0.0;
		cuc = cus = crc = crs = cic = cis = 0.0;
		iodc = 0;
		tgd = 0;
		f0 = f1 = f2 = 0.0;
		AS = false;
		URA = 0.0;
		Status = false;
	}
};

struct BDSEph
{
	unsigned long PRN;
	unsigned long week;
	double URA;
	unsigned long Health;
	double tgd[2];
	unsigned long AODC;
	unsigned long toc;
	double f0, f1, f2;
	unsigned long AODE;
	unsigned long toe;
	double A, ecc, omg, del_n, M0, OMG0, OMGdot, I0, Idot;
	double cuc, cus, crc, crs, cic, cis;
	bool Status;

	BDSEph()
	{
		PRN = 0;
		week = 0;
		URA = 0.0;
		Health = 0;
		tgd[0] = tgd[1] = 0.0;
		AODC = 0;
		toc = 0;
		f0 = f1 = f2 = 0.0;
		AODE = 0;
		toe = 0;
		A = del_n = M0 = ecc = omg = I0 = Idot = OMG0 = OMGdot = 0.0;
		cuc = cus = crc = crs = cic = cis = 0.0;
		Status = false;
	}
};

struct SatRes
{
	unsigned long prn;										//卫星PRN号
	XYZ satpos;												//卫星位置
	NavSys sys;											    //卫星系统
	double elev;											//高度角
	double vx, vy, vz;										//卫星速度
	double dts;												//卫星钟差
	double dv;												//卫星钟速
	double iono;											//电离层改正
	double trop;											//对流层改正
	double tgd;												//卫星设备延迟
	bool Status;											//计算结果有效性

	SatRes()
	{
        prn = 0;
		satpos = XYZ();
        sys = NavSys::UNK;
		elev = 90 * D2R;
		vx = vy = vz = 0.0;
		dts = 0.0;
		dv = 0.0;
		iono = trop = 0.0;
		tgd = 0.0;
		Status = false;
	}
};

struct RcvRes
{
	XYZ rcvpos;
	BLH rcvposblh;
	double dtr[4];
	double vx, vy, vz;
	double vdtr;
	bool Status;
	double sigma;
	double GDOP;
	double PDOP;
	double HDOP;
	double VDOP;
	double velsigma;
	int satnum[4];

	RcvRes()
	{
		rcvpos = XYZ();
		dtr[0] = dtr[1] = dtr[2] = dtr[3] = 0.0;
		vx = vy = vz = 0.0;
		vdtr = 0.0;
		Status = false;
		sigma = velsigma = 999.9;
		GDOP = PDOP = HDOP = VDOP = 999.9;
		satnum[0] = satnum[1] = satnum[2] = satnum[3] = 0;
	}
};

struct LSQInput
{
	NavSys SysType;
	XYZ satpos;
	double P_IF;
	double tgd;
	double satclk;
	double trop;
};

struct SDSatObs
{
    unsigned short prn{0};
    NavSys sys{NavSys::UNK};
    short valid{-1};
    double dP[2]{ 0.0, 0.0 }, dL[2]{ 0.0, 0.0 };
    short n_bas{0}, n_rov{0};
    
    [[nodiscard]] double formMW() const {
        double MW_k1, MW_k2, MW_k3, MW_k4;
        switch (sys)
        {
            case GPS:
            {
                MW_k1 = CLIGHT / (FREQ_L1 - FREQ_L2);
                MW_k2 = -CLIGHT / (FREQ_L1 - FREQ_L2);
                MW_k3 = -FREQ_L1 / (FREQ_L1 + FREQ_L2);
                MW_k4 = -FREQ_L2 / (FREQ_L1 + FREQ_L2);
                break;
            }
            case BDS:
            {
                MW_k1 = CLIGHT / (FREQ_B1 - FREQ_B3);
                MW_k2 = -CLIGHT / (FREQ_B1 - FREQ_B3);
                MW_k3 = -FREQ_B1 / (FREQ_B1 + FREQ_B3);
                MW_k4 = -FREQ_B3 / (FREQ_B1 + FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return MW_k1 * dL[0] + MW_k2 * dL[1] + MW_k3 * dP[0] + MW_k4 * dP[1];
    }
    [[nodiscard]] double formPGF() const {
        if (sys == NavSys::UNK) return 0.0;
        else return dP[0] - dP[1];
    }
    [[nodiscard]] double formLGF() const {
        double lambda1, lambda2;
        switch (sys)
        {
            case GPS:
            {
                lambda1 = CLIGHT / FREQ_L1;
                lambda2 = CLIGHT / FREQ_L2;
                break;
            }
            case BDS:
            {
                lambda1 = CLIGHT / FREQ_B1;
                lambda2 = CLIGHT / FREQ_B3;
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return lambda1 * dL[0] - lambda2 * dL[1];
    }
    [[nodiscard]] double formPIF() const {
        double IF_k1, IF_k2;
        switch (sys)
        {
            case GPS:
            {
                IF_k1 = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                IF_k2 = -(FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                break;
            }
            case BDS:
            {
                IF_k1 = (FREQ_B1 * FREQ_B1) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                IF_k2 = -(FREQ_B3 * FREQ_B3) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return IF_k1 * dP[0] + IF_k2 * dP[1];
    }
    [[nodiscard]] double formLIF() const {
        double IF_k1, IF_k2;
        switch (sys)
        {
            case GPS:
            {
                IF_k1 = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                IF_k2 = -(FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
                break;
            }
            case BDS:
            {
                IF_k1 = (FREQ_B1 * FREQ_B1) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                IF_k2 = -(FREQ_B3 * FREQ_B3) / (FREQ_B1 * FREQ_B1 - FREQ_B3 * FREQ_B3);
                break;
            }
            default:
            {
                return 0.0;
            }
        }
        return IF_k1 * dL[0] + IF_k2 * dL[1];
    }
};

struct SDEpochObs
{
    GPSTime time;
    unsigned int sat_num{0};
    SDSatObs sat_obs[MAXSATNUM];
    LinearCom com_obs[MAXSATNUM];
};

struct DDObs {
    int ref_prn[2]{0,0};
    int b_ref_idx[2]{0,0};
    int r_ref_idx[2]{0,0};
    int dd_sat_num[2]{0,0};
    double fixed_amb[MAXSATNUM * 4]{0.0};
    double res_amb[2]{0.0,0.0}, ratio{0.0};
    float fix_rms[2]{0.0,0.0};
    double dpos[3]{0,0,0};
    bool fixed{false};
};

struct Config {
    int pos_mode{0};        //0-spp, 1-rtk
    bool online{false};     //file or network
    double el_mask{0};      //elevation mask(deg)
    double snr_mask_r{0.0}; //snr mask of rover
    double snr_mask_b{0.0}; //snr mask of base
    int iono_opt{0};        //0-none, 1-dual freq
    int trop_opt{0};        //0-none, 1-hopfield
    double ratio_thres{0.0};//AR ratio threshold
    int sol_format{0};      //solution output format
    bool out_head{false};   //output head of file or not
    bool out_vel{false};    //output velocity or not
    std::string infile_r;   //input rover file
    std::string infile_b;   //input base file
    std::string IP_r;       //IP of rover
    std::string IP_b;       //IP of base
    int port_r{0};          //port of rover
    int port_b{0};          //port of base
    int out_mode{0};        //way to show result
    int basepos_type{0};    //0-llh, 1-xyz
    double basepos[3]{0,0,0};//base position
};
/*-------------------------matrix function-------------------------*/
void VectorAdd(const double* a, const double* b, double* c, int m);
void VectorSub(const double* a, const double* b, double* c, int m);
void VectorScalarMult(const double* a, double* b, double scalar, int m);
double VectorInnerProduct(const double* a, const double* b, int m);
double Norm(const double* a, int m);
void VectorNormalize(const double* a, double* b, int m);
void VectorOuterProduct(const double* a, const double* b, double* c);
void MatrixAdd(const double* a, const double* b, double* c, int m, int n);
void MatrixSub(const double* a, const double* b, double* c, int m, int n);
void MatrixMultiply(const double* a, const double* b, double* c, int m, int n, int s);
void MatrixTranspose(const double* a, double* b, int m, int n);
void MatrixExtend(const double* a, double* b, int m);
void MatrixRowExchange(double* a, int m, int n, int row1, int row2);
void MatrixInverse(const double* a, double* b, int m);
void MatrixPrint(const double* a, int m, int n);
/*----------------------------time function------------------------------*/
ModJulianDay JD2MJD(JulianDay t);
JulianDay MJD2JD(ModJulianDay t);
BDSTime gpst2bdst(GPSTime t);
GPSTime bdst2gpst(BDSTime t);
double gpst2Sec(GPSTime t);
double bdst2Sec(BDSTime t);
JulianDay Common2Julian(CommonDay t);
CommonDay Julian2Common(JulianDay t);
GPSTime MJD2gpst(ModJulianDay t);
ModJulianDay gpst2MJD(GPSTime t);
/*-----------------------------coord function-------------------------------*/
BLH XYZ2BLH(XYZ c, int EllipseType);
XYZ BLH2XYZ(BLH c, int EllipseType);
ENU XYZ2ENU(XYZ refpoint, XYZ point, int EllipseType);
double GetElev(XYZ refpoint, XYZ point, int EllipseType);
double GetDist(XYZ R1, XYZ R2);
double GetAzim(const XYZ& refpoint, const XYZ& point, int EllipseType);
void GetENUTransMat(const XYZ& refpoint, int EllipseType, double Mat[9]);
void GetLOSVector(const XYZ& sat_pos, const XYZ& rcv_pos, double* los);
/*-----------------------------------decoding function-------------------------------------*/
bool DecodeNovOem7(unsigned char buff[], int &LenRem, Range *obs, GPSEph *GEph, BDSEph *CEph, int mode);
bool DecodeRange(unsigned char* buff, Range* obs);
unsigned int crc32(const unsigned char* buff, int len);
bool DecodeGPSEph(unsigned char buff[], GPSEph* GEph);
bool DecodeBDSEph(unsigned char buff[], BDSEph* CEph);
/*----------------------------------------satellite position and velocity------------------------------------------*/
void GPSGetSatpos(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t);
void GPSGetSatclk(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t);
void GPSGetSatvel(int prn, const GPSEph GEph[], SatRes& GRes, GPSTime t);
void BDSGetSatpos(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t);
void BDSGetSatclk(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t);
void BDSGetSatvel(int prn, const BDSEph BEph[], SatRes& BRes, GPSTime t);
void GPSGetSatState(const Obs& obs, const GPSEph GEph[], SatRes& GRes, GPSTime t);
void BDSGetSatState(const Obs& obs, const BDSEph BEph[], SatRes& BRes, GPSTime t);
void EarthRotate(const double cECI[3], double cECEF[3], NavSys sys, double tau);
/*---------------------------------Error Process-------------------------------------*/
void DetectOutlier(Range* rCurr, const Range* rPrev);
double Hopfield(double h, double elev);
/*------------------------------------SPP and SPV----------------------------------------*/
bool SPP(RcvRes *res, SatRes satres[], const RcvRes *init, const Range *range, const GPSEph geph[], const BDSEph beph[],
         const Config& config);
void LSQ(std::vector<LSQInput>& input, RcvRes* rcvres, int gsatnum, int bsatnum);
void SPV(RcvRes* rcvres, const SatRes satres[], const Range* range);
/*-------------------------------------Time Sync--------------------------------------*/
Sync timeSync(const Range& b_range, bool b_status, const Range& r_range, bool r_status);
/*--------------------------------Single Differencing----------------------------------*/
void FormSDEpochObs(const Range& b_obs, const Range& r_obs, SDEpochObs& sd_obs);
void DetectCycleSlip(SDEpochObs& prev_obs, SDEpochObs& curr_obs);
/*--------------------------------Double Differencing----------------------------------*/
bool
DetermineRefSat(const Range &b_range, const Range &r_range, SatRes b_sat_res[], SatRes r_sat_res[],
                const SDEpochObs &sd_obs, DDObs &dd_obs, const Config &config);
/*---------------------------------------RTK-----------------------------------------*/
bool
RTKFloat(RcvRes &r_res, const SatRes b_sat_res[], const SatRes r_sat_res[], const SDEpochObs &sd_obs,
         DDObs &dd_obs, double *a, double *Q_a);
int lambda(int n, int m, const double* a, const double* Q, double* F, double* s);
bool RTKFixed(RcvRes& r_res, const SatRes b_sat_res[], const SatRes r_sat_res[], const SDEpochObs &sd_obs,
              DDObs &dd_obs);
/*------------------------------------Configure---------------------------------------*/
bool ReadConfigureFile(const char* filename, Config& config);
/*----------------------------------Format Control-----------------------------------*/
void OutRcvSol(const RcvRes& rcvres, const Range& range, std::ostream& os = std::cout);
void OutSatSol(const RcvRes& rcvres, const SatRes& satres, const Range& range, std::ostream& os = std::cout);
void OutRTKSol(const GPSTime &b_time, const GPSTime &r_time, const RcvRes &rcv_res, const DDObs &dd_obs,
               const Config &config, std::ostream &ofs = std::cout);
#endif