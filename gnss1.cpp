#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip> // 用于格式化输出

// 为 Windows / Linux 兼容性定义 M_PI
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// WGS-84 椭球参数
const double WGS84_A = 6378137.0; // 椭球长半轴
const double WGS84_F = 1.0 / 298.257223563; // 扁率
const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F; // 第一偏心率的平方
const double WGS84_E_PRIME_SQ = WGS84_E2 / (1.0 - WGS84_E2); // 第二偏心率的平方

// 定义结构体
struct ECEF {
    std::string name;
    double X, Y, Z;
};

struct Geodetic {
    std::string name;
    double B_deg, L_deg, H; // 纬度 (B), 经度 (L) (单位: 度), 大地高 (H)
};

struct GaussKruger {
    std::string name;
    double x, y; // 高斯平面坐标
    int band;      // 带号
};

struct ENU {
    std::string name;
    double E, N, U; // 站心坐标 (East, North, Up)
};

// 辅助函数：弧度转角度
double radToDeg(double rad) {
    return rad * 180.0 / M_PI;
}

// 辅助函数：角度转弧度
double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

/**
 * @brief 任务 1 & 2: WGS-84 空间直角坐标 (X, Y, Z) 转换为 大地坐标 (B, L, H)
 * @param xyz ECEF 坐标
 * @return Geodetic 坐标
 */
Geodetic xyzToBLH(const ECEF& xyz) {
    double L = atan2(xyz.Y, xyz.X);
    double p = sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y);

    double B_i = atan2(xyz.Z, p * (1.0 - WGS84_E2)); // 纬度 B 的初始值
    double H = 0;
    double N = 0;
    
    // 迭代计算 B 和 H
    // 设定一个很小的阈值，用于判断迭代是否收敛
    double tolerance = 1e-12;
    double B_new;

    do {
        double sinB = sin(B_i);
        N = WGS84_A / sqrt(1.0 - WGS84_E2 * sinB * sinB);
        H = p / cos(B_i) - N;
        
        B_new = atan2(xyz.Z, p * (1.0 - WGS84_E2 * N / (N + H)));

        if (std::abs(B_new - B_i) < tolerance) {
            break; // 收敛
        }
        B_i = B_new;

    } while (true);

    // 循环结束后，使用最终的 B_i 重新计算 H
    double sinB_final = sin(B_i);
    N = WGS84_A / sqrt(1.0 - WGS84_E2 * sinB_final * sinB_final);
    H = p / cos(B_i) - N;

    return {xyz.name, radToDeg(B_i), radToDeg(L), H};
}

/**
 * @brief 任务 2: 大地坐标 (B, L, H) 转换为 3度带高斯-克吕格平面坐标 (x, y)
 * @param blh Geodetic 坐标
 * @return GaussKruger 坐标
 */
GaussKruger blhToGaussKruger(const Geodetic& blh) {
    double B_rad = degToRad(blh.B_deg);
    double L_deg = blh.L_deg;

    // 1. 确定带号和中央子午线
    int band = static_cast<int>(round(L_deg / 3.0));
    double L0_rad = degToRad(band * 3.0);
    double L_rad = degToRad(L_deg);
    double l = L_rad - L0_rad; // 经差 (弧度)

    // 2. 计算辅助量
    double sinB = sin(B_rad);
    double cosB = cos(B_rad);
    double tanB = tan(B_rad);
    double t = tanB;
    double t2 = t * t;
    double t4 = t2 * t2;

    double eta_sq = WGS84_E_PRIME_SQ * cosB * cosB;
    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * sinB * sinB);

    // 3. 计算子午线弧长 (X_arc)
    // 使用级数展开
    double e2 = WGS84_E2;
    double m0 = WGS84_A * (1 - e2);
    double m2 = (3.0/2.0) * e2 * m0;
    double m4 = (5.0/4.0) * e2 * m2;
    double m6 = (7.0/6.0) * e2 * m4;
    double m8 = (9.0/8.0) * e2 * m6;

    double a0 = m0 + m2/2.0 + 3.0*m4/8.0 + 5.0*m6/16.0 + 35.0*m8/128.0;
    double a2 = m2/2.0 + m4/2.0 + 15.0*m6/32.0 + 7.0*m8/16.0;
    double a4 = m4/8.0 + 3.0*m6/16.0 + 7.0*m8/32.0;
    double a6 = m6/32.0 + m8/16.0;
    double a8 = m8/128.0;

    double X_arc = a0*B_rad - (a2/2.0)*sin(2*B_rad) + (a4/4.0)*sin(4*B_rad) - (a6/6.0)*sin(6*B_rad) + (a8/8.0)*sin(8*B_rad);

    // 4. 计算高斯平面坐标 (x, y)
    double l2 = l * l;
    double l3 = l2 * l;
    double l4 = l3 * l;
    double l5 = l4 * l;
    double l6 = l5 * l;

    double cosB2 = cosB * cosB;
    double cosB3 = cosB2 * cosB;
    double cosB4 = cosB3 * cosB;
    double cosB5 = cosB4 * cosB;
    double cosB6 = cosB5 * cosB;

    // 计算 x (Northing)
    double x = X_arc;
    x += (N * t * cosB2 * l2 / 2.0);
    x += (N * t * (5.0 - t2 + 9.0 * eta_sq + 4.0 * eta_sq * eta_sq) * cosB4 * l4 / 24.0);
    x += (N * t * (61.0 - 58.0*t2 + t4 + 600.0*eta_sq - 330.0*WGS84_E_PRIME_SQ) * cosB6 * l6 / 720.0);

    // 计算 y (Easting)，不含带号和假东移
    double y_calc = N * cosB * l;
    y_calc += (N * (1.0 - t2 + eta_sq) * cosB3 * l3 / 6.0);
    y_calc += (N * (5.0 - 18.0*t2 + t4 + 14.0*eta_sq - 58.0*eta_sq*t2) * cosB5 * l5 / 120.0);

    // 5. 应用假东移 (500,000 m) 和带号
    // 标准 3 度带 y = (带号 * 1,000,000) + 500,000 + y_calc
    double y = y_calc + 500000.0 + band * 1000000.0;

    return {blh.name, x, y, band};
}

/**
 * @brief 任务 3: WGS-84 空间直角坐标 (X, Y, Z) 转换为 站心坐标 (E, N, U)
 * @param point 要转换的点
 * @param origin 坐标原点 (DPMC)
 * @param origin_blh 坐标原点的大地坐标
 * @return ENU 坐标
 */
ENU xyzToENU(const ECEF& point, const ECEF& origin, const Geodetic& origin_blh) {
    double B0_rad = degToRad(origin_blh.B_deg);
    double L0_rad = degToRad(origin_blh.L_deg);

    double sinB0 = sin(B0_rad);
    double cosB0 = cos(B0_rad);
    double sinL0 = sin(L0_rad);
    double cosL0 = cos(L0_rad);

    // 1. 计算 ECEF 坐标差
    double dX = point.X - origin.X;
    double dY = point.Y - origin.Y;
    double dZ = point.Z - origin.Z;

    // 2. 旋转
    // [ E ]   [ -sinL0         cosL0          0      ] [ dX ]
    // [ N ] = [ -sinB0*cosL0   -sinB0*sinL0   cosB0  ] [ dY ]
    // [ U ]   [  cosB0*cosL0    cosB0*sinL0   sinB0  ] [ dZ ]

    double E = -sinL0 * dX + cosL0 * dY;
    double N = -sinB0 * cosL0 * dX - sinB0 * sinL0 * dY + cosB0 * dZ;
    double U =  cosB0 * cosL0 * dX + cosB0 * sinL0 * dY + sinB0 * dZ;

    return {point.name, E, N, U};
}


int main() {
    // 1. 初始化输入数据 (来自表 1)
    std::vector<ECEF> points = {
        {"DPMC", -2366462.0712, 4813719.8322, 3439397.7856},
        {"MPN1", -2368711.2115, 4812974.2208, 3438897.5484},
        {"MPN2", -2365546.0571, 4815124.9575, 3438094.6390},
        {"MPN3", -2364618.8793, 4814092.9370, 3440139.5922},
        {"MPN4", -2366798.1206, 4814098.4032, 3438639.8475},
        {"MPN5", -2366082.5847, 4812909.4119, 3440782.6667}
    };

    // 存储计算结果
    std::vector<Geodetic> geodetic_coords;
    std::vector<GaussKruger> gauss_coords;
    std::vector<ENU> enu_coords;

    // 2. 执行任务 1 和 2
    for (const auto& point : points) {
        // 任务 1 & 2.1: XYZ -> BLH
        Geodetic blh = xyzToBLH(point);
        geodetic_coords.push_back(blh);

        // 任务 2.2: BLH -> Gauss-Kruger
        GaussKruger gk = blhToGaussKruger(blh);
        gauss_coords.push_back(gk);
    }

    // 3. 执行任务 3
    ECEF origin_xyz = points[0];
    Geodetic origin_blh = geodetic_coords[0];

    // 从 i = 1 开始，因为 DPMC 是原点
    for (size_t i = 1; i < points.size(); ++i) {
        ENU enu = xyzToENU(points[i], origin_xyz, origin_blh);
        enu_coords.push_back(enu);
    }

    // 4. 打印所有结果
    std::cout << "--- 任务 1 & 2: 大地坐标 (B, L, H) 和 3度带高斯平面坐标 (x, y) ---" << std::endl;
    std::cout << std::left << std::setw(6) << "点名"
              << std::setw(16) << "纬度 B (deg)"
              << std::setw(16) << "经度 L (deg)"
              << std::setw(14) << "大地高 H (m)"
              << std::setw(6)  << "带号"
              << std::setw(20) << "x (m)"
              << std::setw(20) << "y (m)" << std::endl;
    std::cout << std::string(98, '-') << std::endl;
    
    std::cout << std::fixed;
    for (size_t i = 0; i < geodetic_coords.size(); ++i) {
        std::cout << std::left << std::setw(6) << geodetic_coords[i].name
                  << std::setprecision(8)
                  << std::setw(16) << geodetic_coords[i].B_deg
                  << std::setw(16) << geodetic_coords[i].L_deg
                  << std::setprecision(4)
                  << std::setw(14) << geodetic_coords[i].H
                  << std::setw(6)  << gauss_coords[i].band
                  << std::setprecision(4)
                  << std::setw(20) << gauss_coords[i].x
                  << std::setw(20) << gauss_coords[i].y << std::endl;
    }

    std::cout << "\n--- 任务 3: 站心坐标 (E, N, U) ---" << std::endl;
    std::cout << "原点: " << origin_xyz.name << std::endl;
    std::cout << std::left << std::setw(6) << "点名"
              << std::setw(18) << "E (m)"
              << std::setw(18) << "N (m)"
              << std::setw(18) << "U (m)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    std::cout << std::fixed << std::setprecision(4);
    for (const auto& enu : enu_coords) {
        std::cout << std::left << std::setw(6) << enu.name
                  << std::setw(18) << enu.E
                  << std::setw(18) << enu.N
                  << std::setw(18) << enu.U << std::endl;
    }

    return 0;
}
