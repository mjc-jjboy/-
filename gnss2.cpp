#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip> // 用于 std::setprecision 和 std::fixed
#include <stdexcept> // 用于 std::stod

// --- 常量定义 ---

// WGS-84 / GRS-80 常量 (用于 GPS)
const double GM_GPS = 3.986005E14; // m^3/s^2
// CGCS2000 / BDS 常量 (与 GRS-80 非常接近)
const double GM_BDS = 3.986004418E14; // m^3/s^2
// WGS-84 地球自转角速度
const double OMEGA_E_DOT = 7.2921151467E-5; // rad/s
// 圆周率
const double PI = 3.141592653589793;

/**
 * @brief 辅助函数：将 'D' 替换为 'E' 以便 std::stod 转换
 */
double d_to_e(std::string s) {
    size_t pos = s.find('D');
    if (pos != std::string::npos) {
        s.replace(pos, 1, "E");
    }
    try {
        return std::stod(s);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error converting string to double: " << s << std::endl;
        return 0.0;
    }
}

/**
 * @brief 存储广播星历参数的结构体
 */
struct Ephemeris {
    std::string satID;
    double gm; // 使用的引力常数

    // 星历参数
    double af0, af1, af2;
    double iode, crs, delta_n, m0;
    double cuc, e, cus, sqrtA;
    double toe, cic, omega0, cis;
    double i0, crc, omega, omega_dot;
    double i_dot;
};

/**
 * @brief 存储 ECEF 坐标
 */
struct ECEF {
    double x, y, z;
};

/**
 * @brief 从图片中加载星历数据
 * @return 包含 G03, C01, C19 星历的向量
 */
std::vector<Ephemeris> populateData() {
    std::vector<Ephemeris> all_eph;
    Ephemeris eph;

    // --- G03 (GPS) ---
    eph.satID = "G03";
    eph.gm = GM_GPS;
    eph.af0 = d_to_e("-4.852190613747D-06");
    eph.af1 = d_to_e("9.348487101919D-12");
    eph.af2 = d_to_e("0.000000000000D+00");
    eph.iode = d_to_e("1.000000000000D+00");
    eph.crs = d_to_e("-6.481250000000D-01");
    eph.delta_n = d_to_e("4.315902391137D-09");
    eph.m0 = d_to_e("1.595863730907D+00");
    eph.cuc = d_to_e("-5.039379461839D-06");
    eph.e = d_to_e("1.691120279374D-09");
    eph.cus = d_to_e("8.972361683846D-06");
    eph.sqrtA = d_to_e("5.153753458808D+03");
    eph.toe = d_to_e("2.016000000000D+05");
    eph.cic = d_to_e("4.312235746159D-09");
    eph.omega0 = d_to_e("5.491223920760D-01");
    eph.cis = d_to_e("-2.793677238460D-08");
    eph.i0 = d_to_e("8.607925160420D-01");
    eph.crc = d_to_e("2.071562500000D-02");
    eph.omega = d_to_e("6.463909090989D-01");
    eph.omega_dot = d_to_e("-8.043217119060D-09");
    eph.i_dot = d_to_e("-2.178914673211D-10");
    all_eph.push_back(eph);

    // --- C01 (BDS) ---
    // 数据来自图片1和图片2的顶部
    eph.satID = "C01";
    eph.gm = GM_BDS;
    eph.af0 = d_to_e("7.481554761119D-04");
    eph.af1 = d_to_e("4.765964001111D-11");
    eph.af2 = d_to_e("0.000000000000D+00");
    eph.iode = d_to_e("1.000000000000D+00");
    eph.crs = d_to_e("6.764062500000D+02");
    eph.delta_n = d_to_e("-5.053781938875D-10");
    eph.m0 = d_to_e("2.808428941360D+00");
    eph.cuc = d_to_e("1.958291977644D-05");
    eph.e = d_to_e("8.439762536444D-04");
    eph.cus = d_to_e("-2.025163984299D-06");
    eph.sqrtA = d_to_e("6.493545307159D+03");
    // C01 数据在第二张图的顶部继续
    eph.toe = d_to_e("2.016000000000D+05");
    eph.cic = d_to_e("1.117587089539D-08");
    eph.omega0 = d_to_e("1.190128731033D+00");
    eph.cis = d_to_e("-3.725290298462D-09");
    eph.i0 = d_to_e("1.067111437234D-01");
    eph.crc = d_to_e("7.660937500000D+01");
    eph.omega = d_to_e("2.963024331937D+00");
    eph.omega_dot = d_to_e("1.433631145063D-09");
    eph.i_dot = d_to_e("1.146176326769D-10");
    all_eph.push_back(eph);

    // --- C19 (BDS) ---
    // 数据来自图片2
    eph.satID = "C19";
    eph.gm = GM_BDS;
    eph.af0 = d_to_e("5.070250794593D-09");
    eph.af1 = d_to_e("5.0191906293D-10");
    eph.af2 = d_to_e("0.000000000000D+00");
    eph.iode = d_to_e("1.000000000000D+00");
    eph.crs = d_to_e("1.890625000000D+01");
    eph.delta_n = d_to_e("3.805158501250D-09");
    eph.m0 = d_to_e("1.744387019202D+00");
    eph.cuc = d_to_e("-3.191176801920D-06");
    eph.e = d_to_e("5.805316613987D-04");
    eph.cus = d_to_e("1.124059781432D-05");
    eph.sqrtA = d_to_e("5.282623334853D+03");
    eph.toe = d_to_e("2.016000000000D+05");
    eph.cic = d_to_e("3.881903171539D-09");
    eph.omega0 = d_to_e("7.913500807130D-01");
    eph.cis = d_to_e("1.024154832077D-08");
    eph.i0 = d_to_e("9.603097794083D-01");
    eph.crc = d_to_e("1.304375000000D+02");
    eph.omega = d_to_e("1.447370829605D+00");
    eph.omega_dot = d_to_e("-6.769210536181D-09");
    eph.i_dot = d_to_e("-2.712971398626D-10");
    all_eph.push_back(eph);

    return all_eph;
}

/**
 * @brief 计算卫星在指定时刻的 ECEF 坐标
 * @param eph 卫星的星历
 * @param t_target 目标计算时间 (GPS 周内秒)
 * @return ECEF 坐标
 */
ECEF calculatePosition(const Ephemeris& eph, double t_target) {
    // 1. 基本参数
    double A = eph.sqrtA * eph.sqrtA; // 轨道长半轴
    double n0 = std::sqrt(eph.gm / (A * A * A)); // 计算的平均运动
    double tk = t_target - eph.toe; // 距离星历参考时刻的时间
    
    // --- 在本题中, t_target == eph.toe, 所以 tk = 0 ---
    // --- 但我们仍按通用公式计算, 0 会自动简化各项 ---

    // 2. 改正后的平均运动
    double n = n0 + eph.delta_n;

    // 3. 归化平近点角
    double Mk = eph.m0 + n * tk;

    // 