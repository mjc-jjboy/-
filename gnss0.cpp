#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
using namespace std;

// ===== 时间系统常量 =====
const time_t GPS_EPOCH = 315964800;   // 1980-01-06 00:00:00 UTC
const time_t BDS_EPOCH = 1136073600;  // 2006-01-01 00:00:00 UTC
const int WEEK_SEC = 604800;
const int BDS_GPS_DIFF = 14;          // BDS 比 GPS 慢14秒

// ===== 工具函数：计算年积日（DOY） =====
int day_of_year(int year, int month, int day) {
    static const int mdays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int doy = day;
    for(int i=0; i<month-1; i++) doy += mdays[i];
    // 闰年修正
    if(month > 2 && ((year%4==0 && year%100!=0) || (year%400==0))) doy++;
    return doy;
}

// ===== 将年月日时分秒转换为 time_t =====
time_t ymdhms_to_time_t(int y, int m, int d, int h, int mi, int s) {
    tm t{};
    t.tm_year = y - 1900;
    t.tm_mon = m - 1;
    t.tm_mday = d;
    t.tm_hour = h;
    t.tm_min = mi;
    t.tm_sec = s;
    t.tm_isdst = 0;
#ifdef _WIN32
    return _mkgmtime(&t);
#else
    return timegm(&t);
#endif
}

// ===== 正算：公历 → GPS/BDS 周秒 =====
void calendar_to_weeksec(int y, int m, int d, int h, int mi, int s) {
    time_t t = ymdhms_to_time_t(y,m,d,h,mi,s);
    int doy = day_of_year(y,m,d);

    // 儒略日计算
    double JD = (t / 86400.0) + 2440587.5;
    double MJD = JD - 2400000.5;

    // GPS 时间
    double gps_sec = difftime(t, GPS_EPOCH);
    int gps_week = gps_sec / WEEK_SEC;
    double gps_tow = fmod(gps_sec, WEEK_SEC);

    // BDS 时间
    double bds_sec = difftime(t, BDS_EPOCH);
    int bds_week = bds_sec / WEEK_SEC;
    double bds_tow = fmod(bds_sec, WEEK_SEC);

    cout << fixed << setprecision(3);
    cout << "======== 正算结果 ========\n";
    cout << "年积日 DOY = " << doy << endl;
    cout << "儒略日 JD = " << JD << endl;
    cout << "约化儒略日 MJD = " << MJD << endl;
    cout << "GPS:  周=" << gps_week << "  秒=" << gps_tow << endl;
    cout << "BDS:  周=" << bds_week << "  秒=" << bds_tow << endl;
}

// ===== 反算：GPS/BDS 周秒 → 公历 =====
void weeksec_to_calendar(int week, double sow, string system) {
    time_t epoch = (system=="BDS") ? BDS_EPOCH : GPS_EPOCH;
    double total_sec = week * WEEK_SEC + sow;
    time_t t = epoch + (time_t)total_sec;

    tm *utc = gmtime(&t);
    int doy = day_of_year(utc->tm_year+1900, utc->tm_mon+1, utc->tm_mday);

    cout << "======== 反算结果 (" << system << ") ========\n";
    cout << "年: " << utc->tm_year+1900 
         << " 月: " << utc->tm_mon+1 
         << " 日: " << utc->tm_mday << endl;
    cout << "时间: " << setfill('0') << setw(2) << utc->tm_hour << ":"
         << setw(2) << utc->tm_min << ":" << setw(2) << utc->tm_sec << endl;
    cout << "年积日 DOY = " << doy << endl;
}

int main() {
    cout << "时间系统转换程序 (C++)\n";
    cout << "选择模式: 1=正算  2=反算\n> ";
    int mode; cin >> mode;

    if(mode == 1){
        int y,m,d,h,mi,s;
        cout << "输入年月日时分秒 (例如 2019 9 10 10 13 37):\n> ";
        cin >> y >> m >> d >> h >> mi >> s;
        calendar_to_weeksec(y,m,d,h,mi,s);
    } 
    else if(mode == 2){
        int week; double sow; string sys;
        cout << "输入系统(GPS/BDS) 周 秒:\n> ";
        cin >> sys >> week >> sow;
        weeksec_to_calendar(week,sow,sys);
    }
    else cout << "输入错误！\n";

    return 0;
}