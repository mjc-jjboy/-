// gnss_opt.cpp
// Course design: Optimal Estimation / GNSS SPP (LS + EKF static + EKF dynamic CV/CA)
// Author: Jincheng Ma
// Supports:
//  - Experiment 1: Least-squares pseudorange Single Point Positioning (SPP)
//  - Experiment 2: Static EKF SPP (state: x y z cb)
//  - Experiment 3: Dynamic EKF SPP (CV: x y z vx vy vz cb cd ; CA: +ax ay az)
//
// Reads:
//  - RINEX observation files (RINEX 2.xx and 3.xx, basic parsing for pseudorange)
//  - RINEX navigation (broadcast ephemeris) RINEX 3 multi-GNSS (GPS, BDS supported; limited Galileo/others)
//
// Notes / Simplifications:
//  - No ionosphere / troposphere correction (can be added if needed).
//  - Uses broadcast ephemeris with standard Kepler propagation, Earth rotation correction.
//  - BDS: supports MEO/IGSO (Kepler) and GEO (with BDS GEO transformation).
//  - Time systems: assumes observation epoch time system from RINEX header if present; supports GPS/BDT.
//    For BDS calculations, ephemeris times are in BDT; conversion uses constant BDT = GPST - 14s.
//
// Build: g++ -O2 -std=c++17 gnss_opt.cpp -o gnss_opt
// Usage examples:
//   ./gnss_opt --exp ls --obs RF5Y0470.25o --nav BRDC00IGS_R_20250470000_01D_MN.rnx --sys C --max-epochs 200 --out exp1_ls_report.txt
//   ./gnss_opt --exp kf_static --obs RF5Y0470.25o --nav BRDC...rnx --sys C --out exp2_static_kf_report.txt
//   ./gnss_opt --exp kf_dyn_cv --obs 9uy70470.25o --nav BRDC...rnx --sys C --out exp3_dyn_cv_report.txt
//   ./gnss_opt --exp kf_dyn_ca --obs 9uy70470.25o --nav BRDC...rnx --sys C --out exp3_dyn_ca_report.txt
//
// Output report: one epoch per line:
//   idx week sow used X Y Z cb_m [sigma0]  (sigma0 only for LS)

#include <bits/stdc++.h>
using namespace std;

// ------------------------ Constants ------------------------
static constexpr double C_LIGHT = 299792458.0;          // m/s
static constexpr double OMEGA_E = 7.2921151467e-5;      // rad/s (WGS-84)
static constexpr double MU_GPS  = 3.986005e14;          // m^3/s^2 (GPS)
static constexpr double MU_BDS  = 3.986004418e14;       // m^3/s^2 (BDS, close to WGS-84)
static constexpr double F_REL   = -4.442807633e-10;     // relativistic correction constant
static constexpr double BDT_GPS_OFFSET = 14.0;          // seconds: BDT = GPST - 14s (so GPST = BDT + 14)
static constexpr int    GPS_WEEK_AT_BDT_EPOCH = 1356;    // GPS week number at 2006-01-01 00:00:00 BDT (BDT week 0)
static constexpr double PI = 3.14159265358979323846;

// ------------------------ Utility ------------------------
static inline double sqr(double x){ return x*x; }

static inline double wrapToPi(double x){
    while(x >  PI) x -= 2*PI;
    while(x < -PI) x += 2*PI;
    return x;
}
static inline int str2int(const string& s){
    try{ return stoi(s); } catch(...) { return 0; }
}
static inline double str2dbl(const string& s){
    try{ return stod(s); } catch(...) { return 0.0; }
}

struct GpsTime {
    int week = 0;
    double sow = 0.0;
};

// Calendar date-time
struct DateTime {
    int y=0, mo=0, d=0, hh=0, mm=0;
    double ss=0;
};

static bool isLeap(int y){
    return (y%4==0 && y%100!=0) || (y%400==0);
}
static int daysInMonth(int y,int m){
    static int dm[] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
    if(m==2) return dm[m] + (isLeap(y)?1:0);
    return dm[m];
}

// Convert calendar date-time (assumed GPST-like continuous) to GPS week & seconds-of-week.
// This uses a simple Julian day conversion.
static GpsTime dateTimeToGps(const DateTime& t){
    // Julian Day for Gregorian calendar
    int Y=t.y, M=t.mo, D=t.d;
    int A = (14 - M)/12;
    int y = Y + 4800 - A;
    int m = M + 12*A - 3;
    long long JDN = D + (153*m + 2)/5 + 365LL*y + y/4 - y/100 + y/400 - 32045;
    double dayFrac = (t.hh - 12)/24.0 + t.mm/1440.0 + t.ss/86400.0; // JD starts at noon
    double JD = (double)JDN + dayFrac;
    // GPS epoch 1980-01-06 00:00:00
    // Compute JD of GPS epoch
    DateTime ge; ge.y=1980; ge.mo=1; ge.d=6; ge.hh=0; ge.mm=0; ge.ss=0;
    int Y2=ge.y, M2=ge.mo, D2=ge.d;
    int A2=(14-M2)/12;
    int y2=Y2+4800-A2;
    int m2=M2+12*A2-3;
    long long JDN0 = D2 + (153*m2 + 2)/5 + 365LL*y2 + y2/4 - y2/100 + y2/400 - 32045;
    double JD0 = (double)JDN0 - 0.5; // 1980-01-06 00:00:00 corresponds to JD with dayFrac -0.5
    // Actually JD at 1980-01-06 00:00:00 is 2444244.5
    // Our JDN0 is 2444245 for that date; subtract 0.5 gives 2444244.5. Good.

    double days = JD - JD0;
    long long w = (long long)floor(days/7.0);
    double sow = (days - 7.0*w) * 86400.0;
    GpsTime gt; gt.week = (int)w; gt.sow = sow;
    return gt;
}

// Adjust sow to be in [-302400, 302400] or [0,604800) etc.
static inline double adjustTime(double t){
    // wrap to +/-302400
    while(t > 302400.0) t -= 604800.0;
    while(t < -302400.0) t += 604800.0;
    return t;
}

// Vector3
struct Vec3 {
    double x=0,y=0,z=0;
    Vec3()=default;
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
};
static inline double dot(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline double norm(const Vec3&a){ return sqrt(dot(a,a)); }

// Rotate vector about Z axis by angle (rad)
static Vec3 rotZ(const Vec3& v, double ang){
    double ca=cos(ang), sa=sin(ang);
    return { ca*v.x + sa*v.y, -sa*v.x + ca*v.y, v.z };
}

// ------------------------ RINEX Navigation (broadcast) ------------------------
struct Eph {
    char sys='G'; // 'G' GPS, 'C' BDS, etc
    int prn=0;
    int week=0; // reference week of toe (GPS week or BDS week)
    double toe=0; // seconds of week
    double toc=0; // seconds of week
    double af0=0, af1=0, af2=0;
    double iode=0, iodc=0;
    double crs=0, crc=0, cus=0, cuc=0, cis=0, cic=0;
    double deltan=0;
    double M0=0;
    double e=0;
    double sqrtA=0;
    double OMG0=0; // Omega0
    double omg=0;  // argument of perigee
    double i0=0;
    double OMGdot=0;
    double idot=0;
    double tgd=0; // group delay (GPS: TGD; BDS: TGD1 maybe)
    double health=0;

    // BDS additional? (ignored): AODE, etc.
    // Flags:
    bool valid=false;
};

struct NavStore {
    // key: (sys, prn) -> list of ephemeris (sorted by toc/toe)
    map<pair<char,int>, vector<Eph>> ephs;
};

// Parse a double from fixed-width substring safely
static double parseDoubleField(const string& line, int start, int len){
    if(start < 0 || start >= (int)line.size()) return 0.0;
    string s = line.substr(start, min(len, (int)line.size()-start));
    for(char& c: s) if(c=='D' || c=='d') c='E';
    // trim
    auto l = s.find_first_not_of(' ');
    if(l==string::npos) return 0.0;
    auto r = s.find_last_not_of(' ');
    s = s.substr(l, r-l+1);
    return str2dbl(s);
}

// Parse int from substring
static int parseIntField(const string& line, int start, int len){
    if(start < 0 || start >= (int)line.size()) return 0;
    string s = line.substr(start, min(len, (int)line.size()-start));
    auto l = s.find_first_not_of(' ');
    if(l==string::npos) return 0;
    auto r = s.find_last_not_of(' ');
    s = s.substr(l, r-l+1);
    return str2int(s);
}

// For RINEX 3 NAV: epoch fields like "G01 2025  2  ..."
static DateTime parseNavEpochRinex3(const string& line){
    DateTime t;
    // sys+prn occupy [0..2], then year at 4..7 etc, but spacing can vary.
    // We'll use stringstream on substring after PRN.
    string rest = line.substr(3);
    stringstream ss(rest);
    ss >> t.y >> t.mo >> t.d >> t.hh >> t.mm >> t.ss;
    return t;
}

// For RINEX 2 NAV: " 1 25  2  ..." (PRN first col)
static DateTime parseNavEpochRinex2(const string& line){
    DateTime t;
    // PRN (2 chars), then year (2), month, day, hour, min, sec
    int yy = parseIntField(line, 3, 2);
    if(yy < 80) t.y = 2000 + yy; else t.y = 1900 + yy;
    t.mo = parseIntField(line, 6, 2);
    t.d  = parseIntField(line, 9, 2);
    t.hh = parseIntField(line, 12, 2);
    t.mm = parseIntField(line, 15, 2);
    t.ss = parseDoubleField(line, 18, 5);
    return t;
}

static bool loadRinexNav(const string& navPath, NavStore& store){
    ifstream in(navPath);
    if(!in) return false;
    string line;
    bool isRinex3 = false;
    // Read header
    while(getline(in,line)){
        if(line.find("RINEX VERSION / TYPE")!=string::npos){
            double ver = parseDoubleField(line,0,9);
            if(ver >= 3.0) isRinex3 = true;
        }
        if(line.find("END OF HEADER")!=string::npos) break;
    }
    // Read body
    while(true){
        if(!getline(in,line)) break;
        if(line.size()<1) continue;
        Eph e;
        DateTime toc_dt;
        if(isRinex3){
            // first char system
            e.sys = line[0];
            e.prn = parseIntField(line,1,2);
            toc_dt = parseNavEpochRinex3(line);
            e.af0 = parseDoubleField(line, 23, 19);
            e.af1 = parseDoubleField(line, 42, 19);
            e.af2 = parseDoubleField(line, 61, 19);
        }else{
            e.sys = 'G'; // RINEX2 nav usually GPS
            e.prn = parseIntField(line,0,2);
            toc_dt = parseNavEpochRinex2(line);
            e.af0 = parseDoubleField(line, 23, 19);
            e.af1 = parseDoubleField(line, 42, 19);
            e.af2 = parseDoubleField(line, 61, 19);
        }
        // Read next 7 lines for GPS-like ephemeris in RINEX3 (8 lines total).
        // For BDS in RINEX3, format is also 8 lines with same main fields.
        string l2,l3,l4,l5,l6,l7,l8;
        if(!getline(in,l2)) break;
        if(!getline(in,l3)) break;
        if(!getline(in,l4)) break;
        if(!getline(in,l5)) break;
        if(!getline(in,l6)) break;
        if(!getline(in,l7)) break;
        if(!getline(in,l8)) break;

        e.iode = parseDoubleField(l2, 4, 19);
        e.crs  = parseDoubleField(l2,23, 19);
        e.deltan = parseDoubleField(l2,42,19);
        e.M0   = parseDoubleField(l2,61,19);

        e.cuc  = parseDoubleField(l3, 4,19);
        e.e    = parseDoubleField(l3,23,19);
        e.cus  = parseDoubleField(l3,42,19);
        e.sqrtA= parseDoubleField(l3,61,19);

        e.toe  = parseDoubleField(l4, 4,19);
        e.cic  = parseDoubleField(l4,23,19);
        e.OMG0 = parseDoubleField(l4,42,19);
        e.cis  = parseDoubleField(l4,61,19);

        e.i0   = parseDoubleField(l5, 4,19);
        e.crc  = parseDoubleField(l5,23,19);
        e.omg  = parseDoubleField(l5,42,19);
        e.OMGdot = parseDoubleField(l5,61,19);

        e.idot = parseDoubleField(l6, 4,19);
        e.health = parseDoubleField(l6,23,19);
        // l6 contains other fields like codes, week, etc, but positions depend by constellation.
        // We'll try to parse week from l6 third/4th field if present.
        int week_guess = (int)parseDoubleField(l6,42,19);
        if(week_guess!=0) e.week = week_guess;

        // l7 has tgd and others; for GPS TGD in first field; for BDS TGD1 may be in first field.
        e.tgd  = parseDoubleField(l7, 4,19);

        // toc, toe in seconds-of-week should be tied to week.
        // Determine toc (seconds-of-week) from toc_dt.
        GpsTime gt_toc = dateTimeToGps(toc_dt);
        e.toc = gt_toc.sow;
        if(e.week==0) e.week = gt_toc.week; // fallback
        e.valid = true;

        // Push to store
        store.ephs[{e.sys,e.prn}].push_back(e);
    }

    // Sort each list by toc
    for(auto &kv: store.ephs){
        auto &v = kv.second;
        sort(v.begin(), v.end(), [](const Eph&a,const Eph&b){
            if(a.week!=b.week) return a.week<b.week;
            return a.toc<b.toc;
        });
    }
    return true;
}

// Select best ephemeris for given time (week,sow) for specific satellite.
static bool selectEph(const NavStore& store, char sys, int prn, int week, double sow, Eph& out){
    auto it = store.ephs.find({sys,prn});
    if(it==store.ephs.end() || it->second.empty()) return false;
    const auto &vec = it->second;
    // choose ephemeris with minimal |t - toc| within 4 hours
    double best = 1e100;
    int bestIdx = -1;
    for(int i=0;i<(int)vec.size();++i){
        const Eph& e = vec[i];
        if(e.week!=week){
            // allow week rollover near boundary
            int dw = week - e.week;
            if(abs(dw)>1) continue;
        }
        double dt = (week - e.week)*604800.0 + (sow - e.toc);
        dt = adjustTime(dt);
        double adt = fabs(dt);
        if(adt < best){
            best=adt; bestIdx=i;
        }
    }
    if(bestIdx<0 || best>4*3600.0) return false;
    out = vec[bestIdx];
    return true;
}

// ------------------------ Satellite position from broadcast ephemeris ------------------------
struct SatState {
    Vec3 pos;        // ECEF (m)
    double clk = 0;  // satellite clock bias (s)
    bool ok=false;
};

// Solve Kepler equation E - e sin E = M
static double solveKepler(double M, double e){
    double E = M;
    for(int k=0;k<20;k++){
        double f = E - e*sin(E) - M;
        double fp = 1 - e*cos(E);
        double dE = -f/fp;
        E += dE;
        if(fabs(dE) < 1e-12) break;
    }
    return E;
}

// Compute satellite clock (s) with broadcast parameters + relativistic correction
static double satClockBias(const Eph& eph, double dt){
    // dt = t - toc (seconds)
    double clk = eph.af0 + eph.af1*dt + eph.af2*dt*dt;
    return clk;
}

// Compute satellite position and clock for GPS-like constellation (GPS, BDS MEO/IGSO).
// Input time: week, sow in the appropriate constellation time system (GPST for GPS, BDT for BDS).
static SatState satPosGpsLike(const Eph& eph, int week, double sow, bool bdsGEO){
    SatState st;
    if(!eph.valid) return st;

    // Time from ephemeris reference epoch
    double tk = (week - eph.week)*604800.0 + (sow - eph.toe);
    tk = adjustTime(tk);

    double mu = (eph.sys=='C') ? MU_BDS : MU_GPS;
    double n0 = sqrt(mu) / pow(eph.sqrtA,3);
    double n = n0 + eph.deltan;

    double M = eph.M0 + n*tk;
    M = wrapToPi(M);

    double E = solveKepler(M, eph.e);
    double sinE = sin(E), cosE = cos(E);

    double v = atan2(sqrt(1 - eph.e*eph.e)*sinE, cosE - eph.e);
    double phi = v + eph.omg;

    // Corrections
    double du = eph.cus*sin(2*phi) + eph.cuc*cos(2*phi);
    double dr = eph.crs*sin(2*phi) + eph.crc*cos(2*phi);
    double di = eph.cis*sin(2*phi) + eph.cic*cos(2*phi);

    double u = phi + du;
    double r = eph.sqrtA*eph.sqrtA*(1 - eph.e*cosE) + dr;
    double i = eph.i0 + di + eph.idot*tk;

    double x_op = r*cos(u);
    double y_op = r*sin(u);

    double OMG = eph.OMG0 + (eph.OMGdot - OMEGA_E)*tk - OMEGA_E*eph.toe;

    Vec3 pos;
    if(eph.sys=='C' && bdsGEO){
        // BDS GEO satellites need special handling (see BDS ICD / IS-GPS-200 adaptation).
        // A common method (RTKLIB-like):
        // 1) compute in Earth-centered inertial-ish with Omega = OMG0 + OMGdot*tk - OMEGA_E*toe
        // 2) then rotate by Earth rotation and inclination offset of 5 degrees.
        // Here we use an approximate transform:
        double sinOMG = sin(OMG), cosOMG = cos(OMG);
        double sinI = sin(i), cosI = cos(i);

        double x = x_op*cosOMG - y_op*cosI*sinOMG;
        double y = x_op*sinOMG + y_op*cosI*cosOMG;
        double z = y_op*sinI;

        // rotate for earth rotation during tk? In GEO, there is additional rotation by OMEGA_E*tk.
        // Apply rotation by OMEGA_E*tk to convert from inertial to ECEF.
        Vec3 p1 = rotZ({x,y,z}, OMEGA_E*tk);

        // Apply inclination offset of 5 deg about X axis (BDS GEO nominal inclination)
        double ang = 5.0*PI/180.0;
        double ca=cos(ang), sa=sin(ang);
        Vec3 p2;
        p2.x = p1.x;
        p2.y = ca*p1.y + sa*p1.z;
        p2.z = -sa*p1.y + ca*p1.z;
        pos = p2;
    }else{
        double cosOMG = cos(OMG), sinOMG = sin(OMG);
        double cosI = cos(i), sinI = sin(i);
        pos.x = x_op*cosOMG - y_op*cosI*sinOMG;
        pos.y = x_op*sinOMG + y_op*cosI*cosOMG;
        pos.z = y_op*sinI;
    }

    // Satellite clock
    double dt = (week - eph.week)*604800.0 + (sow - eph.toc);
    dt = adjustTime(dt);
    double clk = satClockBias(eph, dt);
    // Relativistic correction
    double dtr = F_REL * eph.e * eph.sqrtA * sinE;
    clk += dtr;
    // Group delay (TGD) correction in seconds (meters/c). eph.tgd stored in seconds.
    // Usually applied to pseudorange model as -c*TGD; so we add to clock as +TGD.
    clk -= eph.tgd; // subtract because satellite clock correction includes -TGD for iono-free? We'll keep simple.

    st.pos = pos;
    st.clk = clk;
    st.ok = true;
    return st;
}

// Wrapper for satellite position: handles BDS time system conversions and GEO detection.
static bool computeSatState(const NavStore& nav, char sys, int prn, const GpsTime& obsTimeGPST,
                           const string& obsTimeSys, SatState& st_out){
    // Determine the time in constellation system for ephemeris propagation.
    // GPS modeled in GPST; BDS broadcast ephemeris uses BDT week/sow (week count since 2006-01-01).
    int week = obsTimeGPST.week;
    double sow = obsTimeGPST.sow;

    string ts = obsTimeSys;
    for(char& c: ts) c = toupper((unsigned char)c);

    if(sys=='C'){ // BDS
        if(ts=="BDT"){
            // Observation epoch already in BDT. Convert week count (since 1980) to BDT week (since 2006).
            week = week - GPS_WEEK_AT_BDT_EPOCH;
            // sow already BDT sow
        }else{
            // Assume observation epoch is GPST (common in mixed RINEX). Convert: BDT = GPST - 14s.
            int week_gps = week;
            double sow_bdt = sow - BDT_GPS_OFFSET;
            if(sow_bdt < 0){ sow_bdt += 604800.0; week_gps -= 1; }
            week = week_gps - GPS_WEEK_AT_BDT_EPOCH;
            sow = sow_bdt;
        }
    }else{
        // Non-BDS constellations are modeled in GPST here.
        if(ts=="BDT"){
            // Convert BDT -> GPST
            sow += BDT_GPS_OFFSET;
            if(sow >= 604800.0){ sow -= 604800.0; week += 1; }
        }
    }

    Eph eph;
    if(!selectEph(nav, sys, prn, week, sow, eph)) return false;

    bool isGEO = false;
    if(sys=='C'){
        // BDS GEO PRN: C01-C05, C59-C63 (depending on numbering). In RINEX: C01..C63.
        if(prn<=5 || prn>=59) isGEO = true;
    }

    SatState st = satPosGpsLike(eph, week, sow, isGEO);
    if(!st.ok) return false;
    st_out = st;
    return true;
}

// ------------------------ RINEX Observation parsing ------------------------
struct ObsTypeInfo {
    // For each system, list of observation type strings (e.g., "C1C", "C2I")
    map<char, vector<string>> sysObsTypes;
    // time system (e.g., GPS, BDT)
    string timeSystem; // may be empty
    bool rinex3=false;
    // Optional approximate receiver position from header
    Vec3 approxPos;
    bool hasApprox=false;
};

struct EpochObs {
    GpsTime t_gps;
    DateTime t_cal; // original calendar time
    // map sat -> pseudorange (m)
    // sat key: sys + prn
    vector<pair<pair<char,int>, double>> pranges;
};

// trim right
static inline string rtrim(string s){
    while(!s.empty() && (s.back()=='\r' || s.back()=='\n' || s.back()==' ' || s.back()=='\t')) s.pop_back();
    return s;
}

// Parse RINEX 3 "SYS / # / OBS TYPES" line(s)
static void parseSysObsTypesLineR3(const string& line, ObsTypeInfo& info){
    // Format: [0]=sys, [3..5]=#types, then types in blocks of 4 chars starting at 7.
    if(line.size()<7) return;
    char sys = line[0];
    int n = parseIntField(line,3,3);
    vector<string>& v = info.sysObsTypes[sys];
    for(int i=0;i<13;i++){
        int pos = 7 + i*4;
        if(pos+3 <= (int)line.size()){
            string t = line.substr(pos,3);
            if(t.find_first_not_of(' ')!=string::npos){
                t.erase(remove(t.begin(), t.end(), ' '), t.end());
                v.push_back(t);
            }
        }
    }
    // If more than 13 types, continuation lines will appear; we handle by calling this on each such line.
    // We'll later truncate to declared n.
    if((int)v.size() > n) v.resize(n);
}

static bool loadRinexObs(const string& obsPath, vector<EpochObs>& epochs, ObsTypeInfo& info,
                         const set<char>& sysFilter, int maxEpochs=INT_MAX){
    ifstream in(obsPath);
    if(!in) return false;

    string line;
    bool rinex3=false;
    // header parsing
    while(getline(in,line)){
        line = rtrim(line);
        if(line.find("RINEX VERSION / TYPE")!=string::npos){
            double ver = parseDoubleField(line,0,9);
            if(ver >= 3.0) rinex3=true;
            info.rinex3 = rinex3;
        }
        if(rinex3 && line.find("SYS / # / OBS TYPES")!=string::npos){
            parseSysObsTypesLineR3(line, info);
        }
        if(!rinex3 && line.find("TYPES OF OBSERV")!=string::npos){
            // RINEX2: # / TYPES OF OBSERV
            int n = parseIntField(line,0,6);
            vector<string> types;
            for(int i=0;i<9;i++){
                int pos = 10 + i*6;
                if(pos+2 <= (int)line.size()){
                    string t = line.substr(pos,2);
                    if(t.find_first_not_of(' ')!=string::npos){
                        t.erase(remove(t.begin(), t.end(), ' '), t.end());
                        types.push_back(t);
                    }
                }
            }
            // continuation lines if n>9
            while((int)types.size() < n){
                streampos oldpos = in.tellg();
                string l2;
                if(!getline(in,l2)) break;
                if(l2.find("TYPES OF OBSERV")!=string::npos || l2.size()>=0){
                    for(int i=0;i<9;i++){
                        int pos = 10 + i*6;
                        if(pos+2 <= (int)l2.size()){
                            string t = l2.substr(pos,2);
                            if(t.find_first_not_of(' ')!=string::npos){
                                t.erase(remove(t.begin(), t.end(), ' '), t.end());
                                types.push_back(t);
                            }
                        }
                    }
                }else{
                    in.seekg(oldpos);
                    break;
                }
            }
            // RINEX2 types apply to GPS by default; we'll map to 'G' and 'C' (best effort)
            info.sysObsTypes['G'] = types;
            info.sysObsTypes['C'] = types;
        }
        if(line.find("TIME SYSTEM ID")!=string::npos){
            // RINEX3 header line: positions 0..2
            string ts = line.substr(0,3);
            ts.erase(remove(ts.begin(), ts.end(), ' '), ts.end());
            info.timeSystem = ts;
        }
                if(line.find("APPROX POSITION XYZ")!=string::npos){
            // First 3 values are X Y Z (m)
            string part = line.substr(0, 60);
            stringstream ss(part);
            double X=0,Y=0,Z=0;
            if(ss >> X >> Y >> Z){
                info.approxPos = Vec3(X,Y,Z);
                info.hasApprox = true;
            }
        }
if(line.find("TIME OF FIRST OBS")!=string::npos){
            // Often last 3 chars are time system (GPS, BDT)
            if(line.size()>=51){
                string ts = line.substr(48,3);
                ts.erase(remove(ts.begin(), ts.end(), ' '), ts.end());
                if(!ts.empty()) info.timeSystem = ts;
            }
        }
        if(line.find("END OF HEADER")!=string::npos) break;
    }

    // Choose default pseudorange type per system.
    // For RINEX3: prefer any type starting with 'C' (code pseudorange), in priority list.
    auto choosePRTypeR3 = [&](char sys)->int{
        auto it = info.sysObsTypes.find(sys);
        if(it==info.sysObsTypes.end()) return -1;
        const auto& t = it->second;
        // prefer C1C, C1I, C2I, C2X, C6I etc
        vector<string> pref = {"C1C","C1I","C1X","C2I","C2X","C2C","C6I","C7I","C5Q","C5X","C1P","C2P"};
        for(const auto& p: pref){
            for(int i=0;i<(int)t.size();++i){
                if(t[i]==p) return i;
            }
        }
        // fallback: first type that starts with C
        for(int i=0;i<(int)t.size();++i){
            if(!t[i].empty() && t[i][0]=='C') return i;
        }
        return -1;
    };
    auto choosePRTypeR2 = [&](char sys)->int{
        auto it = info.sysObsTypes.find(sys);
        if(it==info.sysObsTypes.end()) return -1;
        const auto& t = it->second;
        // prefer C1, P1, C2, P2
        vector<string> pref = {"C1","P1","C2","P2"};
        for(const auto& p: pref){
            for(int i=0;i<(int)t.size();++i){
                if(t[i]==p) return i;
            }
        }
        return (t.empty()? -1:0);
    };

    // Observation body
    int epochCount = 0;
    while(epochCount < maxEpochs && getline(in,line)){
        line = rtrim(line);
        if(line.empty()) continue;

        EpochObs ep;
        vector<pair<pair<char,int>, double>> pr;
        if(rinex3){
            if(line[0] != '>') continue;
            // > yyyy mm dd hh mm ss.sssss  flag  nsat
            DateTime dt;
            stringstream ss(line.substr(1));
            int flag=0, nsat=0;
            ss >> dt.y >> dt.mo >> dt.d >> dt.hh >> dt.mm >> dt.ss >> flag >> nsat;
            ep.t_cal = dt;
            ep.t_gps = dateTimeToGps(dt);

            // read nsat lines of observation
            for(int i=0;i<nsat;i++){
                string ol;
                if(!getline(in,ol)) break;
                if(ol.size()<3) continue;
                char sys = ol[0];
                int prn = parseIntField(ol,1,2);
                if(!sysFilter.empty() && sysFilter.find(sys)==sysFilter.end()) continue;

                int idx = choosePRTypeR3(sys);
                if(idx<0) continue;
                // Each obs field is 16 chars, starting at pos 3
                int pos = 3 + idx*16;
                if(pos+14 <= (int)ol.size()){
                    string field = ol.substr(pos,14);
                    double val = str2dbl(field);
                    if(val > 1e3){ // in meters
                        pr.push_back({{sys,prn}, val});
                    }
                }
            }
        }else{
            // RINEX2 epoch line: yy mm dd hh mm ss flag nsat [sat list ...]
            // Example: 25  2 15  0  0  0.0000000  0 12G01G03...
            if(line.size()<32) continue;
            DateTime dt;
            int yy = parseIntField(line,0,2);
            dt.y = (yy<80)?(2000+yy):(1900+yy);
            dt.mo = parseIntField(line,3,2);
            dt.d  = parseIntField(line,6,2);
            dt.hh = parseIntField(line,9,2);
            dt.mm = parseIntField(line,12,2);
            dt.ss = parseDoubleField(line,15,11);
            int flag = parseIntField(line,28,1);
            int nsat = parseIntField(line,29,3);
            ep.t_cal = dt;
            ep.t_gps = dateTimeToGps(dt);

            // Satellite list may continue
            vector<pair<char,int>> sats;
            int satsOnLine = min(nsat, (int)((line.size()-32)/3));
            for(int i=0;i<satsOnLine;i++){
                int pos = 32 + i*3;
                char sys = 'G';
                int prn = parseIntField(line,pos,2);
                // In RINEX2, sys not explicit. We'll assume GPS for now; if user filters 'C', still accept.
                if(!sysFilter.empty()){
                    if(sysFilter.find(sys)==sysFilter.end() && sysFilter.find('C')==sysFilter.end()) continue;
                }
                sats.push_back({sys,prn});
            }
            while((int)sats.size() < nsat){
                string l2;
                streampos oldpos = in.tellg();
                if(!getline(in,l2)) break;
                if(l2.size()<1){ in.seekg(oldpos); break; }
                // continuation has up to 12 sats
                int count = min(nsat - (int)sats.size(), (int)(l2.size()/3));
                for(int i=0;i<count;i++){
                    int pos = i*3;
                    char sys='G';
                    int prn = parseIntField(l2,pos,2);
                    sats.push_back({sys,prn});
                }
            }

            // For each sat, read observation record: #types fields, 80 chars per line, 16 char field
            int ntypes = 0;
            if(info.sysObsTypes.count('G')) ntypes = (int)info.sysObsTypes['G'].size();
            int idx = choosePRTypeR2('G');
            if(idx<0) idx = 0;

            int linesPerSat = (ntypes + 4) / 5; // 5 obs per line (each 16 chars)
            for(auto sp: sats){
                char sys = sp.first;
                int prn = sp.second;
                string rec;
                for(int k=0;k<linesPerSat;k++){
                    string rline;
                    if(!getline(in,rline)) break;
                    if((int)rline.size()<80) rline.resize(80,' ');
                    rec += rline;
                }
                int pos = idx*16;
                if(pos+14 <= (int)rec.size()){
                    double val = str2dbl(rec.substr(pos,14));
                    if(val > 1e3){
                        // Map sys as 'G' in RINEX2; allow user to request 'G' or 'C' - we still keep 'G'.
                        pr.push_back({{sys,prn}, val});
                    }
                }
            }
        }

        if(pr.size()>=4){
            ep.pranges = std::move(pr);
            epochs.push_back(std::move(ep));
            epochCount++;
        }
    }

    return !epochs.empty();
}

// ------------------------ Least Squares SPP ------------------------
struct LSResult {
    Vec3 rx;
    double cb_m = 0;   // receiver clock bias in meters
    double sigma0 = 0;
    int used = 0;
    bool ok=false;
};

struct LSConfig {
    int maxIter = 10;
    double tol = 1e-4; // meters
};

// Solve 4x4 linear system by Gaussian elimination
static bool solve4(double A[4][4], double b[4], double x[4]){
    // Augmented matrix
    double M[4][5];
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++) M[i][j]=A[i][j];
        M[i][4]=b[i];
    }
    for(int col=0;col<4;col++){
        // pivot
        int piv = col;
        double best = fabs(M[col][col]);
        for(int r=col+1;r<4;r++){
            if(fabs(M[r][col])>best){
                best=fabs(M[r][col]); piv=r;
            }
        }
        if(best < 1e-20) return false;
        if(piv!=col){
            for(int j=col;j<5;j++) swap(M[piv][j], M[col][j]);
        }
        // normalize
        double diag = M[col][col];
        for(int j=col;j<5;j++) M[col][j] /= diag;
        // eliminate
        for(int r=0;r<4;r++){
            if(r==col) continue;
            double f = M[r][col];
            for(int j=col;j<5;j++) M[r][j] -= f*M[col][j];
        }
    }
    for(int i=0;i<4;i++) x[i]=M[i][4];
    return true;
}

static LSResult sppLeastSquares(const NavStore& nav, const EpochObs& ep, const string& timeSys,
                                char sysWanted, Vec3 x0, double cb0_m, const LSConfig& cfg){
    LSResult res;
    Vec3 x = x0;
    double cb = cb0_m;

    for(int iter=0; iter<cfg.maxIter; ++iter){
        // Build normal equations: (H^T W H) dx = H^T W v
        double N[4][4] = {{0}};
        double u[4] = {0,0,0,0};
        double vv_sum = 0;
        int m_used=0;

        for(const auto& ob: ep.pranges){
            char sys = ob.first.first;
            int prn  = ob.first.second;
            double P = ob.second;

            if(sysWanted!='A' && sys!=sysWanted) continue;

            SatState st;
            if(!computeSatState(nav, sys, prn, ep.t_gps, timeSys, st)) continue;

            // Approx geometric range
            Vec3 rs = st.pos;
            Vec3 dr = rs - x;
            double rho = norm(dr);
            if(rho < 1.0) continue;

            // Earth rotation correction using signal transit time
            double tau = rho / C_LIGHT;
            Vec3 rs_corr = rotZ(rs, OMEGA_E * tau); // rotate sat position to account for Earth rotation
            dr = rs_corr - x;
            rho = norm(dr);

            // Predicted pseudorange
            double pred = rho + cb - C_LIGHT*st.clk;

            double v = P - pred; // residual (m)

            // Design matrix row (1x4): partials wrt x,y,z,cb
            double hx = -(dr.x)/rho;
            double hy = -(dr.y)/rho;
            double hz = -(dr.z)/rho;
            double hcb = 1.0;

            // Weight: simple elevation-based weight could be added; here W=1
            double w = 1.0;

            double h[4] = {hx,hy,hz,hcb};
            for(int i=0;i<4;i++){
                for(int j=0;j<4;j++){
                    N[i][j] += w * h[i]*h[j];
                }
                u[i] += w * h[i]*v;
            }
            vv_sum += w*v*v;
            m_used++;
        }

        if(m_used < 4) { res.ok=false; return res; }

        double dx[4];
        if(!solve4(N,u,dx)){ res.ok=false; return res; }

        x.x += dx[0];
        x.y += dx[1];
        x.z += dx[2];
        cb  += dx[3];

        double step = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
        if(step < cfg.tol && fabs(dx[3]) < cfg.tol){
            // compute sigma0
            int dof = m_used - 4;
            double sigma0 = (dof>0)? sqrt(vv_sum/dof) : 0;
            res.rx = x; res.cb_m = cb; res.sigma0 = sigma0; res.used = m_used; res.ok=true;
            return res;
        }
        if(iter==cfg.maxIter-1){
            int dof = m_used - 4;
            double sigma0 = (dof>0)? sqrt(vv_sum/dof) : 0;
            res.rx = x; res.cb_m = cb; res.sigma0 = sigma0; res.used = m_used; res.ok=true;
            return res;
        }
    }
    res.ok=false;
    return res;
}

// ------------------------ EKF ------------------------
struct KFResult {
    Vec3 rx;
    double cb_m=0;
    double cd_mps=0;
    int used=0;
    bool ok=false;
};

static void matMul(const vector<vector<double>>&A, const vector<vector<double>>&B, vector<vector<double>>&C){
    int n=A.size(), m=B[0].size(), k=B.size();
    C.assign(n, vector<double>(m,0.0));
    for(int i=0;i<n;i++){
        for(int kk=0;kk<k;kk++){
            double a=A[i][kk];
            for(int j=0;j<m;j++) C[i][j]+=a*B[kk][j];
        }
    }
}
static void matAddInPlace(vector<vector<double>>&A, const vector<vector<double>>&B){
    for(size_t i=0;i<A.size();i++) for(size_t j=0;j<A[i].size();j++) A[i][j]+=B[i][j];
}
static vector<vector<double>> matTranspose(const vector<vector<double>>&A){
    int n=A.size(), m=A[0].size();
    vector<vector<double>>T(m, vector<double>(n,0.0));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) T[j][i]=A[i][j];
    return T;
}
static bool matInv(vector<vector<double>> A, vector<vector<double>>& inv){
    // Gauss-Jordan for small matrices
    int n=A.size();
    inv.assign(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) inv[i][i]=1.0;
    for(int col=0;col<n;col++){
        int piv=col;
        double best=fabs(A[col][col]);
        for(int r=col+1;r<n;r++){
            if(fabs(A[r][col])>best){ best=fabs(A[r][col]); piv=r; }
        }
        if(best<1e-20) return false;
        if(piv!=col){
            swap(A[piv], A[col]);
            swap(inv[piv], inv[col]);
        }
        double diag=A[col][col];
        for(int j=0;j<n;j++){ A[col][j]/=diag; inv[col][j]/=diag; }
        for(int r=0;r<n;r++){
            if(r==col) continue;
            double f=A[r][col];
            for(int j=0;j<n;j++){
                A[r][j]-=f*A[col][j];
                inv[r][j]-=f*inv[col][j];
            }
        }
    }
    return true;
}
static vector<double> matVecMul(const vector<vector<double>>&A, const vector<double>&x){
    int n=A.size(), m=A[0].size();
    vector<double> y(n,0.0);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++) y[i]+=A[i][j]*x[j];
    }
    return y;
}
static vector<double> vecAdd(const vector<double>&a,const vector<double>&b){
    vector<double> c=a;
    for(size_t i=0;i<c.size();i++) c[i]+=b[i];
    return c;
}
static vector<double> vecSub(const vector<double>&a,const vector<double>&b){
    vector<double> c=a;
    for(size_t i=0;i<c.size();i++) c[i]-=b[i];
    return c;
}

// EKF update step for pseudorange measurements
static void ekfUpdate(vector<double>&x, vector<vector<double>>&P,
                      const vector<double>&z, const vector<double>&z_pred,
                      const vector<vector<double>>&H,
                      const vector<double>&Rdiag){
    int n = (int)x.size();
    int m = (int)z.size();

    // y = z - z_pred
    vector<double> y(m,0.0);
    for(int i=0;i<m;i++) y[i] = z[i] - z_pred[i];

    // S = H P H^T + R
    vector<vector<double>> HP, S;
    matMul(H, P, HP);                 // m x n
    vector<vector<double>> Ht = matTranspose(H); // n x m
    matMul(HP, Ht, S);                // m x m
    for(int i=0;i<m;i++) S[i][i] += Rdiag[i];

    // K = P H^T S^{-1}
    vector<vector<double>> S_inv;
    if(!matInv(S, S_inv)) return;
    vector<vector<double>> PHt;
    matMul(P, Ht, PHt);               // n x m
    vector<vector<double>> K;
    matMul(PHt, S_inv, K);            // n x m

    // x = x + K y
    vector<double> Ky(n,0.0);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++) Ky[i] += K[i][j]*y[j];
    }
    for(int i=0;i<n;i++) x[i] += Ky[i];

    // P = (I - K H) P
    vector<vector<double>> KH, IminusKH;
    matMul(K, H, KH);                 // n x n
    IminusKH.assign(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            double v = (i==j?1.0:0.0) - KH[i][j];
            IminusKH[i][j]=v;
        }
    }
    vector<vector<double>> newP;
    matMul(IminusKH, P, newP);
    P.swap(newP);
}

// EKF prediction: x = F x ; P = F P F^T + Q
static void ekfPredict(vector<double>&x, vector<vector<double>>&P,
                       const vector<vector<double>>&F,
                       const vector<vector<double>>&Q){
    vector<double> x_new = matVecMul(F, x);
    x.swap(x_new);

    vector<vector<double>> FP, FPFt;
    matMul(F, P, FP);
    vector<vector<double>> Ft = matTranspose(F);
    matMul(FP, Ft, FPFt);
    matAddInPlace(FPFt, Q);
    P.swap(FPFt);
}

// Build pseudorange measurement model for given state (pos + clock bias).
static void buildMeasModel(const NavStore& nav, const EpochObs& ep, const string& timeSys, char sysWanted,
                           const Vec3& rx, double cb_m,
                           vector<double>&z, vector<double>&z_pred, vector<vector<double>>&H,
                           vector<double>&Rdiag, int stateDim, int idx_cb, int idx_cd){
    z.clear(); z_pred.clear(); H.clear(); Rdiag.clear();

    for(const auto& ob: ep.pranges){
        char sys = ob.first.first;
        int prn  = ob.first.second;
        double Pm = ob.second;

        if(sysWanted!='A' && sys!=sysWanted) continue;

        SatState st;
        if(!computeSatState(nav, sys, prn, ep.t_gps, timeSys, st)) continue;

        Vec3 rs = st.pos;
        Vec3 dr = rs - rx;
        double rho = norm(dr);
        if(rho < 1.0) continue;

        // Sagnac/Earth rotation correction
        double tau = rho / C_LIGHT;
        Vec3 rs_corr = rotZ(rs, OMEGA_E * tau);
        dr = rs_corr - rx;
        rho = norm(dr);

        double pred = rho + cb_m - C_LIGHT*st.clk;

        // H row
        vector<double> h(stateDim, 0.0);
        h[0] = -(dr.x)/rho;
        h[1] = -(dr.y)/rho;
        h[2] = -(dr.z)/rho;
        h[idx_cb] = 1.0;
        // pseudorange doesn't directly observe clock drift; leave 0 (idx_cd) if present

        z.push_back(Pm);
        z_pred.push_back(pred);
        H.push_back(std::move(h));

        // Measurement noise variance: use 3m^2 default
        double sig = 3.0;
        Rdiag.push_back(sig*sig);
    }
}

// ------------------------ Experiments ------------------------
struct Options {
    string exp="kf_dyn_cv"; // Experiment 3 default: CV
    string obsPath;
    string navPath;
    char sys='A'; // 'A' all, 'G', 'C'
    int maxEpochs=200;
    string outPath="exp3_dyn_cv_report.txt";
};
;

static void printUsage(){
    cerr << "Usage: exp3_dyn_kf --obs <obs.rnx> --nav <nav.rnx>\n"
         << "                 [--exp <kf_dyn_cv|kf_dyn_ca>] [--sys <A|G|C>] [--max-epochs N] [--out report.txt]\n";
}


// Initial receiver position guess: Earth center offset
static Vec3 defaultInitPos(){
    return Vec3(1113194.9, -4841695.5, 3985354.7); // rough anywhere; user can tune
}

// Run Experiment 1: LS SPP for each epoch (independent)
static bool runExpLS(const Options& opt, const vector<EpochObs>& epochs, const ObsTypeInfo& info, const NavStore& nav, const Vec3& initPos){
    ofstream out(opt.outPath);
    if(!out) return false;

    Vec3 x = initPos;
    double cb=0.0;

    LSConfig cfg; cfg.maxIter=10; cfg.tol=1e-4;

    for(size_t i=0;i<epochs.size();i++){
        const auto& ep = epochs[i];
        LSResult r = sppLeastSquares(nav, ep, info.timeSystem, opt.sys, x, cb, cfg);
        if(r.ok){
            x = r.rx; cb = r.cb_m;
            out << (i+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << r.used << " "
                << setprecision(4) << r.rx.x << " " << r.rx.y << " " << r.rx.z << " "
                << r.cb_m << " " << r.sigma0 << "\n";
        }else{
            out << (i+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << 0 << " "
                << "nan nan nan nan nan\n";
        }
    }
    return true;
}

// Run Experiment 2: Static EKF (x y z cb)
static bool runExpKFStatic(const Options& opt, const vector<EpochObs>& epochs, const ObsTypeInfo& info, const NavStore& nav, const Vec3& initPos){
    ofstream out(opt.outPath);
    if(!out) return false;

    int n=4;
    vector<double> x(n,0.0);
    Vec3 x0 = initPos;
    x[0]=x0.x; x[1]=x0.y; x[2]=x0.z; x[3]=0.0;

    vector<vector<double>> P(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) P[i][i] = 100.0*100.0; // 100m std

    // Static: F=I, Q small (random walk)
    vector<vector<double>> F(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) F[i][i]=1.0;

    for(size_t k=0;k<epochs.size();k++){
        const auto& ep = epochs[k];

        // Q: tiny
        vector<vector<double>> Q(n, vector<double>(n,0.0));
        double qpos = 0.01*0.01; // m^2 per epoch
        double qcb  = 10.0*10.0; // (m)^2 per epoch
        Q[0][0]=Q[1][1]=Q[2][2]=qpos;
        Q[3][3]=qcb;

        ekfPredict(x,P,F,Q);

        Vec3 rx(x[0],x[1],x[2]);
        double cb = x[3];

        vector<double> z,zp,Rdiag;
        vector<vector<double>> H;
        buildMeasModel(nav, ep, info.timeSystem, opt.sys, rx, cb, z, zp, H, Rdiag, n, 3, -1);

        int used = (int)z.size();
        if(used>=4){
            ekfUpdate(x,P,z,zp,H,Rdiag);

            out << (k+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << used << " "
                << setprecision(4) << x[0] << " " << x[1] << " " << x[2] << " "
                << x[3] << "\n";
        }else{
            out << (k+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << 0 << " nan nan nan nan\n";
        }
    }
    return true;
}

// Build CV model matrices for dt (seconds). State: [x y z vx vy vz cb cd]
static void buildCVModel(double dt, vector<vector<double>>&F, vector<vector<double>>&Q){
    int n=8;
    F.assign(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) F[i][i]=1.0;
    for(int i=0;i<3;i++){
        F[i][i+3] = dt;
    }
    // clock: cb += cd*dt
    F[6][7] = dt;

    // Process noise: continuous white accel with spectral density sa^2, and clock random walk
    double sa = 0.5; // m/s^2 (tune)
    double q = sa*sa;

    Q.assign(n, vector<double>(n,0.0));
    for(int axis=0;axis<3;axis++){
        int p=axis, v=axis+3;
        Q[p][p] = q*dt*dt*dt/3.0;
        Q[p][v] = q*dt*dt/2.0;
        Q[v][p] = q*dt*dt/2.0;
        Q[v][v] = q*dt;
    }
    // clock bias/drift noise
    double s_cb = 10.0;   // m
    double s_cd = 1.0;    // m/s
    Q[6][6] = s_cb*s_cb*dt;
    Q[7][7] = s_cd*s_cd*dt;
}

// Build CA model matrices for dt. State: [x y z vx vy vz ax ay az cb cd] (11)
static void buildCAModel(double dt, vector<vector<double>>&F, vector<vector<double>>&Q){
    int n=11;
    F.assign(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++) F[i][i]=1.0;

    // position updates
    for(int i=0;i<3;i++){
        F[i][i+3] = dt;
        F[i][i+6] = 0.5*dt*dt;
        F[i+3][i+6] = dt;
    }
    // clock
    F[9][10] = dt;

    // Process noise: white jerk for acceleration
    double sj = 0.2; // m/s^3
    double q = sj*sj;

    Q.assign(n, vector<double>(n,0.0));
    for(int axis=0;axis<3;axis++){
        int p=axis, v=axis+3, a=axis+6;
        double dt2=dt*dt, dt3=dt2*dt, dt4=dt3*dt, dt5=dt4*dt;
        Q[p][p] = q*dt5/20.0;
        Q[p][v] = q*dt4/8.0;
        Q[p][a] = q*dt3/6.0;

        Q[v][p] = Q[p][v];
        Q[v][v] = q*dt3/3.0;
        Q[v][a] = q*dt2/2.0;

        Q[a][p] = Q[p][a];
        Q[a][v] = Q[v][a];
        Q[a][a] = q*dt;
    }
    double s_cb = 10.0, s_cd = 1.0;
    Q[9][9] = s_cb*s_cb*dt;
    Q[10][10] = s_cd*s_cd*dt;
}

// Run Experiment 3: Dynamic EKF (CV or CA)
static bool runExpKFDynamic(const Options& opt, const vector<EpochObs>& epochs, const ObsTypeInfo& info,
                            const NavStore& nav, bool useCA, const Vec3& initPos){
    ofstream out(opt.outPath);
    if(!out) return false;

    int n = useCA?11:8;
    vector<double> x(n,0.0);
    Vec3 x0 = initPos;
    x[0]=x0.x; x[1]=x0.y; x[2]=x0.z;
    if(!useCA){
        // vx vy vz
        x[3]=x[4]=x[5]=0.0;
        x[6]=0.0; // cb
        x[7]=0.0; // cd
    }else{
        x[3]=x[4]=x[5]=0.0;
        x[6]=x[7]=x[8]=0.0; // ax ay az
        x[9]=0.0; x[10]=0.0;
    }

    vector<vector<double>> P(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++){
        double s = 100.0;
        if(!useCA){
            if(i>=3 && i<=5) s=10.0; // velocity
            if(i==6) s=100.0; // cb
            if(i==7) s=10.0;  // cd
        }else{
            if(i>=3 && i<=5) s=10.0;
            if(i>=6 && i<=8) s=1.0;
            if(i==9) s=100.0;
            if(i==10) s=10.0;
        }
        P[i][i]=s*s;
    }

    double prevSow = epochs.front().t_gps.sow;
    int prevWeek = epochs.front().t_gps.week;

    for(size_t k=0;k<epochs.size();k++){
        const auto& ep = epochs[k];
        double dt = 1.0;
        // Compute dt from previous epoch (handle week rollover)
        if(k>0){
            double t1 = prevWeek*604800.0 + prevSow;
            double t2 = ep.t_gps.week*604800.0 + ep.t_gps.sow;
            dt = t2 - t1;
            if(dt<=0 || dt>30) dt = 1.0; // fallback
        }
        prevSow = ep.t_gps.sow;
        prevWeek = ep.t_gps.week;

        vector<vector<double>> F,Q;
        if(useCA) buildCAModel(dt,F,Q);
        else buildCVModel(dt,F,Q);

        ekfPredict(x,P,F,Q);

        Vec3 rx(x[0],x[1],x[2]);
        int idx_cb = useCA?9:6;
        int idx_cd = useCA?10:7;
        double cb = x[idx_cb];

        vector<double> z,zp,Rdiag;
        vector<vector<double>> H;
        buildMeasModel(nav, ep, info.timeSystem, opt.sys, rx, cb, z, zp, H, Rdiag, n, idx_cb, idx_cd);

        int used = (int)z.size();
        if(used>=4){
            ekfUpdate(x,P,z,zp,H,Rdiag);

            out << (k+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << used << " "
                << setprecision(4) << x[0] << " " << x[1] << " " << x[2] << " "
                << x[idx_cb] << " " << x[idx_cd] << "\n";
        }else{
            out << (k+1) << " " << ep.t_gps.week << " " << fixed << setprecision(3) << ep.t_gps.sow
                << " " << 0 << " nan nan nan nan nan\n";
        }
    }
    return true;
}

// ------------------------ CLI ------------------------
static bool parseArgs(int argc, char**argv, Options& opt){
    for(int i=1;i<argc;i++){
        string a=argv[i];
        if(a=="--exp" && i+1<argc){ opt.exp=argv[++i]; }
        else if(a=="--obs" && i+1<argc){ opt.obsPath=argv[++i]; }
        else if(a=="--nav" && i+1<argc){ opt.navPath=argv[++i]; }
        else if(a=="--sys" && i+1<argc){
            string s=argv[++i];
            if(!s.empty()) opt.sys = (char)toupper((unsigned char)s[0]);
        }
        else if(a=="--max-epochs" && i+1<argc){ opt.maxEpochs=stoi(argv[++i]); }
        else if(a=="--out" && i+1<argc){ opt.outPath=argv[++i]; }
        else if(a=="-h"||a=="--help"){ return false; }
        else{
            cerr<<"Unknown arg: "<<a<<"\n";
            return false;
        }
    }
    if(opt.obsPath.empty() || opt.navPath.empty()) return false;
    return true;
}

int main(int argc, char**argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Options opt;
    if(!parseArgs(argc, argv, opt)){
        printUsage();
        return 1;
    }

    NavStore nav;
    if(!loadRinexNav(opt.navPath, nav)){
        cerr << "Failed to load NAV: " << opt.navPath << "\n";
        return 2;
    }

    set<char> sysFilter;
    if(opt.sys!='A') sysFilter.insert(opt.sys);
    else{
        // Allow all systems (G,C), but keep only those present in obs parsing
    }

    vector<EpochObs> epochs;
    ObsTypeInfo info;
    if(!loadRinexObs(opt.obsPath, epochs, info, sysFilter, opt.maxEpochs)){
        cerr << "Failed to load OBS or no usable epochs: " << opt.obsPath << "\n";
        return 3;
    }
    if(info.timeSystem.empty()) info.timeSystem = "GPS"; // default
// Experiment fixed: Dynamic EKF SPP (CV/CA)
// Use: --exp kf_dyn_cv (default) or --exp kf_dyn_ca
bool isCA = (opt.exp=="kf_dyn_ca");
if(opt.exp!="kf_dyn_cv" && opt.exp!="kf_dyn_ca"){
    cerr << "Note: --exp ignored for this program. Using kf_dyn_cv by default.\n";
    isCA = false;
}
if(opt.outPath=="exp3_dyn_cv_report.txt" && isCA){
    opt.outPath = "exp3_dyn_ca_report.txt";
}

bool ok=false;
Vec3 initPos = info.hasApprox ? info.approxPos : defaultInitPos();
ok = runExpKFDynamic(opt, epochs, info, nav, isCA, initPos);

if(!ok){
        cerr << "Failed to run experiment or write output.\n";
        return 4;
    }

    cerr << "Done. Wrote report: " << opt.outPath << "\n";
    return 0;
}
