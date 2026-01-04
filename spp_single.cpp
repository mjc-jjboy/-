#include <bits/stdc++.h>
using namespace std;

// ---------------- 常量 ----------------
static constexpr double C = 299792458.0;
static constexpr double MU_GPS = 3.986005e14;
static constexpr double OMEGA_E = 7.2921151467e-5;
static constexpr double F_REL = -4.442807633e-10;

static constexpr double A_WGS84 = 6378137.0;
static constexpr double F_WGS84 = 1.0 / 298.257223563;
static constexpr double E2_WGS84 = F_WGS84 * (2.0 - F_WGS84);

// ---------------- 小工具 ----------------
static inline string rtrim(string s){
    while(!s.empty() && (s.back()=='\n' || s.back()=='\r' || s.back()==' ' || s.back()=='\t')) s.pop_back();
    return s;
}
static inline string trim(const string& s){
    size_t a = 0, b = s.size();
    while(a<b && isspace((unsigned char)s[a])) ++a;
    while(b>a && isspace((unsigned char)s[b-1])) --b;
    return s.substr(a,b-a);
}
static inline string dexp_to_exp(string s){
    for(char &ch: s){
        if(ch=='D' || ch=='d') ch='E';
    }
    return s;
}
static inline double parse_double_field(const string& s){
    string t = trim(dexp_to_exp(s));
    if(t.empty()) return numeric_limits<double>::quiet_NaN();
    try { return stod(t); } catch(...) { return numeric_limits<double>::quiet_NaN(); }
}

// ---------------- 民用日期到天数 (Howard Hinnant) ----------------
// 返回自1970-01-01起的天数
static inline int64_t days_from_civil(int y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y-399) / 400;
    const unsigned yoe = (unsigned)(y - era * 400);      // [0, 399]
    const unsigned doy = (153*(m + (m > 2 ? -3 : 9)) + 2)/5 + d-1;  // [0, 365]
    const unsigned doe = yoe * 365 + yoe/4 - yoe/100 + doy;         // [0, 146096]
    return (int64_t)era * 146097 + (int64_t)doe - 719468; // 719468 = days to 1970-01-01
}

struct GpsTime { int week{}; double tow{}; };

static inline int64_t unix_seconds_from_civil(int y,int m,int d,int hh,int mm,double ss){
    int64_t days = days_from_civil(y,(unsigned)m,(unsigned)d);
    int64_t sec = days*86400LL + hh*3600LL + mm*60LL + (int64_t)floor(ss);
    return sec;
}

static inline GpsTime gps_week_tow(int y,int m,int d,int hh,int mm,double ss){
    // GPS纪元: 1980-01-06 00:00:00
    static const int64_t gps_epoch_unix = unix_seconds_from_civil(1980,1,6,0,0,0.0);
    int64_t t_unix = unix_seconds_from_civil(y,m,d,hh,mm,ss);
    double frac = ss - floor(ss);
    double diff = (double)(t_unix - gps_epoch_unix) + frac;
    int week = (int)floor(diff / 604800.0);
    double tow = diff - week*604800.0;
    tow = fmod(tow, 604800.0);
    if(tow < 0) tow += 604800.0;
    return {week, tow};
}

static inline double adjust_week_time(double dt){
    if(dt > 302400.0) dt -= 604800.0;
    else if(dt < -302400.0) dt += 604800.0;
    return dt;
}

// ---------------- ECEF/LLH/ENU ----------------
struct Vec3 { double x{},y{},z{}; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,double k){ return {a.x*k,a.y*k,a.z*k}; }
static inline double dot(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline double norm(const Vec3&a){ return sqrt(dot(a,a)); }

static inline void ecef_to_llh(const Vec3& p, double &lat, double &lon, double &h){
    lon = atan2(p.y, p.x);
    double r = hypot(p.x, p.y);
    lat = atan2(p.z, r * (1.0 - E2_WGS84));
    for(int it=0; it<10; ++it){
        double sinlat = sin(lat);
        double N = A_WGS84 / sqrt(1.0 - E2_WGS84*sinlat*sinlat);
        h = r / cos(lat) - N;
        double lat_new = atan2(p.z, r * (1.0 - E2_WGS84 * N / (N + h)));
        if(fabs(lat_new - lat) < 1e-12) { lat = lat_new; break; }
        lat = lat_new;
    }
    double sinlat = sin(lat);
    double N = A_WGS84 / sqrt(1.0 - E2_WGS84*sinlat*sinlat);
    h = r / cos(lat) - N;
}

static inline array<array<double,3>,3> ecef_to_enu_matrix(double lat, double lon){
    double slat = sin(lat), clat = cos(lat);
    double slon = sin(lon), clon = cos(lon);
    return {{
        {{-slon,        clon,       0.0}},
        {{-slat*clon,  -slat*slon,  clat}},
        {{ clat*clon,   clat*slon,  slat}}
    }};
}

static inline Vec3 mat_mul(const array<array<double,3>,3>&R, const Vec3&v){
    return {
        R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z,
        R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z,
        R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z
    };
}

static inline Vec3 matT_mul(const array<array<double,3>,3>&R, const Vec3&v){
    return {
        R[0][0]*v.x + R[1][0]*v.y + R[2][0]*v.z,
        R[0][1]*v.x + R[1][1]*v.y + R[2][1]*v.z,
        R[0][2]*v.x + R[1][2]*v.y + R[2][2]*v.z
    };
}

static inline double elevation_angle(const Vec3& rx, const Vec3& sat){
    double lat, lon, h;
    ecef_to_llh(rx, lat, lon, h);
    auto R = ecef_to_enu_matrix(lat, lon);
    Vec3 enu = mat_mul(R, sat - rx);
    double horiz = hypot(enu.x, enu.y);
    return atan2(enu.z, horiz);
}

static inline double simple_tropo_delay(double el){
    double sinel = max(sin(el), 0.1);
    return 2.3 / sinel;
}

// ---------------- RINEX OBS ----------------
struct ObsHeader {
    vector<string> gps_types;
    Vec3 approx_xyz{};
    bool has_approx=false;
    // antenna delta H/E/N (marker->APC)
    double antH=0.0, antE=0.0, antN=0.0;
};

static inline vector<string> split_ws(const string& s){
    istringstream iss(s);
    vector<string> out;
    string t;
    while(iss >> t) out.push_back(t);
    return out;
}

static ObsHeader parse_rinex3_obs_header(ifstream &in){
    ObsHeader h;
    string line;
    while(getline(in, line)){
        line = rtrim(line);
        if(line.find("APPROX POSITION XYZ") != string::npos){
            h.approx_xyz.x = parse_double_field(line.substr(0,14));
            h.approx_xyz.y = parse_double_field(line.substr(14,14));
            h.approx_xyz.z = parse_double_field(line.substr(28,14));
            h.has_approx = true;
        } else if(line.find("ANTENNA: DELTA H/E/N") != string::npos){
            h.antH = parse_double_field(line.substr(0,14));
            h.antE = parse_double_field(line.substr(14,14));
            h.antN = parse_double_field(line.substr(28,14));
        } else if(line.find("SYS / # / OBS TYPES") != string::npos){
            char sys = line.empty()?'\0':line[0];
            if(sys != 'G') continue;
            int n = 0;
            try { n = stoi(line.substr(3,3)); } catch(...) { n = 0; }
            vector<string> types = split_ws(line.substr(7, 53));
            while((int)types.size() < n){
                string l2;
                if(!getline(in, l2)) break;
                l2 = rtrim(l2);
                auto more = split_ws(l2.substr(7, 53));
                types.insert(types.end(), more.begin(), more.end());
            }
            if((int)types.size() > n) types.resize(n);
            h.gps_types = types;
        } else if(line.find("END OF HEADER") != string::npos){
            break;
        }
    }
    return h;
}

struct EpochTime { int y,m,d,hh,mm; double ss; };

struct EpochObs {
    EpochTime t{};
    unordered_map<string, vector<double>> satObs;
};

static bool read_next_epoch(ifstream &in, const ObsHeader& hdr, EpochObs &out){
    string line;
    while(getline(in, line)){
        if(line.empty()) continue;
        if(line[0] != '>') continue;
        line = rtrim(line);
        char gt;
        istringstream iss(line);
        iss >> gt;
        if(!iss) return false;
        iss >> out.t.y >> out.t.m >> out.t.d >> out.t.hh >> out.t.mm >> out.t.ss;
        int flag=0, nsat=0;
        iss >> flag >> nsat;
        out.satObs.clear();

        int ntypes = (int)hdr.gps_types.size();
        for(int i=0;i<nsat;i++){
            string sline;
            if(!getline(in, sline)) return true;
            sline = rtrim(sline);
            if(sline.size() < 3) continue;

            string sat = trim(sline.substr(0,3));
            if(sat.empty() || sat[0] != 'G') {
                continue; // mixed systems: we only parse GPS
            }

            string fields = (sline.size() > 3) ? sline.substr(3) : "";

            // If a file is strictly RINEX fixed-width, a sat record may wrap to multiple lines.
            // Continuation lines start with 3 spaces.
            while((int)fields.size() < 16*ntypes){
                streampos pos = in.tellg();
                string cont;
                if(!getline(in, cont)) break;
                cont = rtrim(cont);
                if(cont.size()>=3 && cont[0]==' ' && cont[1]==' ' && cont[2]==' '){
                    fields += cont.substr(3);
                } else {
                    in.seekg(pos);
                    break;
                }
            }

            vector<double> vals(ntypes, numeric_limits<double>::quiet_NaN());
            for(int k=0;k<ntypes;k++){
                size_t start = (size_t)k*16;
                if(start+16 > fields.size()) break;
                string chunk = fields.substr(start, 16);
                string vstr = chunk.substr(0,14);
                vals[k] = parse_double_field(vstr);
            }
            out.satObs.emplace(sat, move(vals));
        }
        return true;
    }
    return false;
}

// ---------------- RINEX NAV (GPS only) ----------------
struct GpsEph {
    string prn; // "G01"
    int tocWeek{};
    double tocTow{}; // seconds-of-week
    double af0{}, af1{}, af2{};
    double IODE{}, Crs{}, dn{}, M0{};
    double Cuc{}, e{}, Cus{}, sqrtA{};
    double Toe{}, Cic{}, Omega0{}, Cis{};
    double i0{}, Crc{}, omega{}, Omegadot{};
    double IDOT{}, GPSWeek{}, TGD{};
};

static vector<double> parse_nav_4fields(const string& line){
    // RINEX fixed columns: [5-23][24-42][43-61][62-80] (1-based)
    vector<double> out(4, numeric_limits<double>::quiet_NaN());
    string s = line;
    if((int)s.size() < 80) s += string(80 - s.size(), ' ');
    int idxs[4] = {4,23,42,61};
    for(int i=0;i<4;i++){
        out[i] = parse_double_field(s.substr(idxs[i], 19));
    }
    return out;
}

static unordered_map<string, vector<GpsEph>> parse_rinex3_nav_gps(const string& path){
    ifstream in(path);
    if(!in) throw runtime_error("Cannot open nav: "+path);
    string line;
    while(getline(in, line)){
        if(line.find("END OF HEADER") != string::npos) break;
    }

    unordered_map<string, vector<GpsEph>> ephs;

    while(getline(in, line)){
        line = rtrim(line);
        if(line.empty()) continue;
        char sys = line[0];
        int rec_len = (sys=='R') ? 4 : 8;
        vector<string> rec;
        rec.reserve(rec_len);
        rec.push_back(line);
        for(int i=1;i<rec_len;i++){
            string l2;
            if(!getline(in, l2)) break;
            rec.push_back(rtrim(l2));
        }
        if(sys != 'G' || (int)rec.size() < 8) continue;

        // First line MUST be fixed-width parsed (fields may be adjacent with no spaces)
        string l0 = rec[0];
        if((int)l0.size() < 80) l0 += string(80 - l0.size(), ' ');
        string sat = trim(l0.substr(0,3));
        if(sat.empty()) continue;
        int y  = stoi(l0.substr(4,4));
        int m  = stoi(l0.substr(9,2));
        int d  = stoi(l0.substr(12,2));
        int hh = stoi(l0.substr(15,2));
        int mm = stoi(l0.substr(18,2));
        double ss = parse_double_field(l0.substr(21,2));

        GpsEph e;
        e.prn = sat;
        auto gt = gps_week_tow(y,m,d,hh,mm,ss);
        e.tocWeek = gt.week;
        e.tocTow  = gt.tow;
        e.af0 = parse_double_field(l0.substr(23,19));
        e.af1 = parse_double_field(l0.substr(42,19));
        e.af2 = parse_double_field(l0.substr(61,19));

        auto f1 = parse_nav_4fields(rec[1]);
        e.IODE = f1[0]; e.Crs=f1[1]; e.dn=f1[2]; e.M0=f1[3];
        auto f2 = parse_nav_4fields(rec[2]);
        e.Cuc=f2[0]; e.e=f2[1]; e.Cus=f2[2]; e.sqrtA=f2[3];
        auto f3 = parse_nav_4fields(rec[3]);
        e.Toe=f3[0]; e.Cic=f3[1]; e.Omega0=f3[2]; e.Cis=f3[3];
        auto f4 = parse_nav_4fields(rec[4]);
        e.i0=f4[0]; e.Crc=f4[1]; e.omega=f4[2]; e.Omegadot=f4[3];
        auto f5 = parse_nav_4fields(rec[5]);
        e.IDOT=f5[0];
        e.GPSWeek=f5[2]; // third field is GPS week
        auto f6 = parse_nav_4fields(rec[6]);
        e.TGD=f6[2]; // third field

        ephs[e.prn].push_back(e);
    }

    for(auto &kv : ephs){
        auto &v = kv.second;
        sort(v.begin(), v.end(), [](const GpsEph&a,const GpsEph&b){
            if(a.tocWeek != b.tocWeek) return a.tocWeek < b.tocWeek;
            return a.tocTow < b.tocTow;
        });
    }
    return ephs;
}

static const GpsEph* select_eph(const vector<GpsEph>& list, int week, double tow, double &best_abs){
    best_abs = 1e18;
    const GpsEph* best = nullptr;
    for(const auto& e : list){
        double dt = (week - e.tocWeek)*604800.0 + (tow - e.Toe);
        dt = adjust_week_time(dt);
        double ad = fabs(dt);
        if(ad < best_abs){
            best_abs = ad;
            best = &e;
        }
    }
    return best;
}

static pair<Vec3,double> sat_pos_clock(const GpsEph& eph, int week, double tow){
    // position at tow, clock correction (s)
    double A = eph.sqrtA * eph.sqrtA;
    double n0 = sqrt(MU_GPS / (A*A*A));
    double n = n0 + eph.dn;

    double tk = adjust_week_time(tow - eph.Toe);
    double M = eph.M0 + n * tk;

    // Kepler
    double E = M;
    for(int i=0;i<20;i++){
        double E_new = M + eph.e * sin(E);
        if(fabs(E_new - E) < 1e-12) { E = E_new; break; }
        E = E_new;
    }
    double sinE = sin(E), cosE = cos(E);
    double dtr = F_REL * eph.e * eph.sqrtA * sinE;

    double v = atan2(sqrt(1 - eph.e*eph.e) * sinE, cosE - eph.e);
    double phi = v + eph.omega;

    double sin2 = sin(2*phi), cos2 = cos(2*phi);
    double u = phi + eph.Cuc * cos2 + eph.Cus * sin2;
    double r = A * (1 - eph.e * cosE) + eph.Crc * cos2 + eph.Crs * sin2;
    double i = eph.i0 + eph.IDOT * tk + eph.Cic * cos2 + eph.Cis * sin2;

    double x_p = r * cos(u);
    double y_p = r * sin(u);

    double Omega = eph.Omega0 + (eph.Omegadot - OMEGA_E) * tk - OMEGA_E * eph.Toe;

    double cosO = cos(Omega), sinO = sin(Omega);
    double cosi = cos(i), sini = sin(i);

    Vec3 sat;
    sat.x = x_p * cosO - y_p * cosi * sinO;
    sat.y = x_p * sinO + y_p * cosi * cosO;
    sat.z = y_p * sini;

    // satellite clock
    double dt = (week - eph.tocWeek)*604800.0 + (tow - eph.tocTow);
    dt = adjust_week_time(dt);
    double dt_sv = eph.af0 + eph.af1*dt + eph.af2*dt*dt + dtr - eph.TGD;
    return {sat, dt_sv};
}

// ---------------- WLS solver ----------------
struct EpochSolution {
    Vec3 x;
    double dt_r{};
    double rms{};
    int used{};
};

static bool solve_linear_4x4(double N[4][4], double u[4], double x[4]){
    // Gaussian elimination with partial pivot
    double A[4][5];
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++) A[i][j]=N[i][j];
        A[i][4]=u[i];
    }
    for(int col=0; col<4; col++){
        int piv = col;
        double best = fabs(A[col][col]);
        for(int r=col+1;r<4;r++){
            double v = fabs(A[r][col]);
            if(v > best){ best=v; piv=r; }
        }
        if(best < 1e-20) return false;
        if(piv != col){
            for(int j=col;j<5;j++) swap(A[piv][j], A[col][j]);
        }
        double diag = A[col][col];
        for(int j=col;j<5;j++) A[col][j] /= diag;
        for(int r=0;r<4;r++){
            if(r==col) continue;
            double factor = A[r][col];
            if(fabs(factor) < 1e-30) continue;
            for(int j=col;j<5;j++) A[r][j] -= factor * A[col][j];
        }
    }
    for(int i=0;i<4;i++) x[i]=A[i][4];
    return true;
}

static EpochSolution* solve_epoch_spp(
    const EpochObs& ep,
    const ObsHeader& hdr,
    const unordered_map<string, vector<GpsEph>>& ephs,
    const Vec3& x0,
    double dt0,
    const string& obs_code,
    int max_iter = 10
){
    if(hdr.gps_types.empty()) return nullptr;

    // choose obs index
    int idx = -1;
    for(int i=0;i<(int)hdr.gps_types.size();i++){
        if(hdr.gps_types[i] == obs_code) { idx=i; break; }
    }
    if(idx < 0){
        for(int i=0;i<(int)hdr.gps_types.size();i++){
            if(hdr.gps_types[i].size()>=2 && hdr.gps_types[i][0]=='C' && hdr.gps_types[i][1]=='1'){
                idx=i; break;
            }
        }
    }
    if(idx < 0) return nullptr;

    auto gt = gps_week_tow(ep.t.y, ep.t.m, ep.t.d, ep.t.hh, ep.t.mm, ep.t.ss);
    int week = gt.week;
    double tow_rx = gt.tow;

    // collect sats (GPS only + has ephemeris + has P)
    vector<string> sats;
    vector<double> Ps;
    sats.reserve(ep.satObs.size());
    Ps.reserve(ep.satObs.size());
    for(const auto& kv : ep.satObs){
        const string& sat = kv.first;
        if(sat.empty() || sat[0] != 'G') continue;
        auto itE = ephs.find(sat);
        if(itE == ephs.end()) continue;
        const auto& vec = kv.second;
        if(idx >= (int)vec.size()) continue;
        double P = vec[idx];
        if(!isfinite(P) || P < 1e5) continue;
        sats.push_back(sat);
        Ps.push_back(P);
    }
    if((int)sats.size() < 4) return nullptr;

    Vec3 x = x0;
    double dt_r = dt0;

    for(int iter=0; iter<max_iter; iter++){
        double N[4][4]{};
        double u[4]{};
        int used = 0;

        for(size_t si=0; si<sats.size(); si++){
            const string& sat = sats[si];
            double P = Ps[si];

            const auto& evec = ephs.at(sat);
            double best_abs;
            const GpsEph* eph = select_eph(evec, week, tow_rx, best_abs);
            if(!eph || !(best_abs <= 4*3600.0)) continue;

            // transmit time iteration
            double tau = P / C;
            double t_tx = tow_rx - tau;
            t_tx = fmod(t_tx, 604800.0);
            if(t_tx < 0) t_tx += 604800.0;

            Vec3 sat_pos{};
            double dt_sv = 0.0;
            for(int k=0;k<2;k++){
                auto pc = sat_pos_clock(*eph, week, t_tx);
                sat_pos = pc.first;
                dt_sv = pc.second;

                double rho_corr = P - C*(dt_r - dt_sv);
                tau = rho_corr / C;
                t_tx = tow_rx - tau;
                t_tx = fmod(t_tx, 604800.0);
                if(t_tx < 0) t_tx += 604800.0;
            }
            auto pc = sat_pos_clock(*eph, week, t_tx);
            sat_pos = pc.first;
            dt_sv = pc.second;

            // Earth rotation correction
            double ang = OMEGA_E * tau;
            double ca = cos(ang), sa = sin(ang);
            Vec3 sat_rot{
                ca*sat_pos.x + sa*sat_pos.y,
                -sa*sat_pos.x + ca*sat_pos.y,
                sat_pos.z
            };

            double rho = norm(sat_rot - x);
            if(rho < 1.0) continue;

            double el = elevation_angle(x, sat_rot);
            double trop = simple_tropo_delay(el);

            double P_hat = rho + C*(dt_r - dt_sv) + trop;
            double v = P - P_hat;

            Vec3 diff = x - sat_rot;
            double ux = diff.x / rho;
            double uy = diff.y / rho;
            double uz = diff.z / rho;

            double sinel = max(sin(el), 0.1);
            double w = sinel*sinel;

            double arow[4] = {ux, uy, uz, C};
            for(int r=0;r<4;r++){
                for(int c=0;c<4;c++) N[r][c] += w * arow[r] * arow[c];
                u[r] += w * arow[r] * v;
            }
            used++;
        }

        if(used < 4) return nullptr;

        double dx[4];
        if(!solve_linear_4x4(N, u, dx)) return nullptr;

        x.x += dx[0];
        x.y += dx[1];
        x.z += dx[2];
        dt_r += dx[3];

        double dpos = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        if(dpos < 1e-4 && fabs(dx[3]) < 1e-10) break;
    }

    // post-fit RMS
    vector<double> res;
    int used = 0;
    for(size_t si=0; si<sats.size(); si++){
        const string& sat = sats[si];
        double P = Ps[si];
        const auto& evec = ephs.at(sat);
        double best_abs;
        const GpsEph* eph = select_eph(evec, week, tow_rx, best_abs);
        if(!eph || !(best_abs <= 4*3600.0)) continue;

        double tau = P / C;
        double t_tx = tow_rx - tau;
        t_tx = fmod(t_tx, 604800.0); if(t_tx < 0) t_tx += 604800.0;

        Vec3 sat_pos{};
        double dt_sv = 0.0;
        for(int k=0;k<2;k++){
            auto pc = sat_pos_clock(*eph, week, t_tx);
            sat_pos = pc.first;
            dt_sv = pc.second;
            double rho_corr = P - C*(dt_r - dt_sv);
            tau = rho_corr / C;
            t_tx = tow_rx - tau;
            t_tx = fmod(t_tx, 604800.0); if(t_tx < 0) t_tx += 604800.0;
        }
        auto pc = sat_pos_clock(*eph, week, t_tx);
        sat_pos = pc.first; dt_sv = pc.second;

        double ang = OMEGA_E * tau;
        double ca = cos(ang), sa = sin(ang);
        Vec3 sat_rot{ca*sat_pos.x + sa*sat_pos.y, -sa*sat_pos.x + ca*sat_pos.y, sat_pos.z};

        double rho = norm(sat_rot - x);
        double el = elevation_angle(x, sat_rot);
        double trop = simple_tropo_delay(el);
        double P_hat = rho + C*(dt_r - dt_sv) + trop;
        double v = P - P_hat;
        res.push_back(v);
        used++;
    }

    double rms = numeric_limits<double>::quiet_NaN();
    if(!res.empty()){
        double s=0; for(double v:res) s += v*v;
        rms = sqrt(s / res.size());
    }
    return new EpochSolution{x, dt_r, rms, used};
}

static Vec3 apply_marker_offset(const Vec3& x_apc, const ObsHeader& hdr){
    // marker = apc - (ENU->ECEF * [E,N,H])
    double lat, lon, h;
    ecef_to_llh(x_apc, lat, lon, h);
    auto R = ecef_to_enu_matrix(lat, lon); // ECEF->ENU
    Vec3 enu{hdr.antE, hdr.antN, hdr.antH};
    Vec3 shift_ecef = matT_mul(R, enu);
    return x_apc - shift_ecef;
}

// ---------------- main ----------------
static void usage(){
    cerr << "Usage:\n"
         << "  spp_ls <obs.25o> <nav.25p> [--out out.csv] [--obs_code C1C] [--marker]\n\n"
         << "Notes:\n"
         << "  - GPS only (system 'G')\n"
         << "  - WLS with sin(el)^2 weighting, simple tropo, sat clock + relativity + TGD, Earth rotation correction\n";
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if(argc < 3){
        usage();
        return 1;
    }
    string obs_path = argv[1];
    string nav_path = argv[2];
    string out_path = "spp_out.csv";
    string obs_code = "C1C";
    bool marker = false;

    for(int i=3;i<argc;i++){
        string a = argv[i];
        if(a=="--out" && i+1<argc){ out_path = argv[++i]; }
        else if(a=="--obs_code" && i+1<argc){ obs_code = argv[++i]; }
        else if(a=="--marker"){ marker = true; }
        else {
            cerr << "Unknown arg: " << a << "\n";
            usage();
            return 1;
        }
    }

    ifstream obs_in(obs_path);
    if(!obs_in){
        cerr << "Cannot open obs: " << obs_path << "\n";
        return 1;
    }
    ObsHeader hdr = parse_rinex3_obs_header(obs_in);
    if(!hdr.has_approx){
        cerr << "No APPROX POSITION XYZ found in observation header.\n";
        return 1;
    }
    if(hdr.gps_types.empty()){
        cerr << "No GPS observation types found (SYS / # / OBS TYPES for G).\n";
        return 1;
    }

    unordered_map<string, vector<GpsEph>> ephs;
    try {
        ephs = parse_rinex3_nav_gps(nav_path);
    } catch(const exception& e){
        cerr << e.what() << "\n";
        return 1;
    }
    if(ephs.empty()){
        cerr << "No GPS ephemeris parsed from nav.\n";
        return 1;
    }

    ofstream out(out_path);
    if(!out){
        cerr << "Cannot write out: " << out_path << "\n";
        return 1;
    }

    out << "time,x,y,z,lat_deg,lon_deg,h_m,clk_bias_s,rms_m,nsat\n";

    Vec3 x_init = hdr.approx_xyz;
    double dt_init = 0.0;

    EpochObs ep;
    int solved = 0, total = 0;

    long long nmean = 0;
    long double sx=0, sy=0, sz=0;

    while(read_next_epoch(obs_in, hdr, ep)){
        total++;
        auto sol = solve_epoch_spp(ep, hdr, ephs, x_init, dt_init, obs_code);
        if(!sol) continue;
        solved++;

        Vec3 x = sol->x;
        double clk = sol->dt_r;
        x_init = x;
        dt_init = clk;

        Vec3 x_out = marker ? apply_marker_offset(x, hdr) : x;
        double lat, lon, h;
        ecef_to_llh(x_out, lat, lon, h);

        ostringstream ts;
        ts << setfill('0')
           << setw(4) << ep.t.y << "-" << setw(2) << ep.t.m << "-" << setw(2) << ep.t.d
           << "T" << setw(2) << ep.t.hh << ":" << setw(2) << ep.t.mm << ":";
        ts << fixed << setprecision(3) << setw(6) << ep.t.ss;

        out.setf(std::ios::fixed);
        out << ts.str() << ","
            << setprecision(4) << x_out.x << "," << x_out.y << "," << x_out.z << ","
            << setprecision(10) << (lat*180.0/M_PI) << "," << (lon*180.0/M_PI) << ","
            << setprecision(4) << h << ",";
        out << scientific << setprecision(12) << clk << ",";
        out << fixed << setprecision(4) << sol->rms << "," << sol->used << "\n";

        sx += x_out.x; sy += x_out.y; sz += x_out.z; nmean++;
    }

    out.close();

    cerr << "=== SPP done ===\n";
    cerr << "epochs solved: " << solved << " / " << total << "\n";
    cerr << "output: " << out_path << "\n";

    if(nmean > 0){
        Vec3 mean{(double)(sx/nmean), (double)(sy/nmean), (double)(sz/nmean)};
        double lat, lon, h;
        ecef_to_llh(mean, lat, lon, h);
        cerr.setf(std::ios::fixed);
        cerr << setprecision(4);
        cerr << "mean ECEF (m): " << mean.x << ", " << mean.y << ", " << mean.z << "\n";
        cerr << setprecision(10);
        cerr << "mean BLH: " << (lat*180.0/M_PI) << ", " << (lon*180.0/M_PI) << ", " << h << "\n";
    }

    return 0;
}
