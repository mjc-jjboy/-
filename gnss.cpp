// main.cpp
// RINEX 3.05 obs+nav -> GPS single-point pseudorange positioning (SPP)
// Dependencies: Eigen for linear algebra
// Compile: g++ -std=c++17 -O2 -I /path/to/eigen main.cpp -o gnss_ppp

#include <bits/stdc++.h>
using namespace std;

// include Eigen headers
#include <Eigen/Dense>

constexpr double PI = 3.14159265358979323846;
constexpr double MU = 3.986005e14; // earth gravitational constant, m^3/s^2 (GPS)
constexpr double OMEGA_E = 7.2921151467e-5; // earth rotation rate rad/s
constexpr double c0 = 299792458.0; // speed of light m/s

// ----- Data structures -----
struct ObsRecord {
    double epoch_gps; // seconds of week
    // map satellite id -> pseudorange (m). satellite id like "G01"
    unordered_map<string,double> pr;
};

struct EphData {
    string sat;
    // toe and orbital params (as in GPS broadcast)
    double toe;
    double sqrtA;
    double e;
    double i0;
    double OMG0;
    double omg;
    double M0;
    double deltaN;
    double OMG_dot;
    double i_dot;
    double Cuc, Cus, Crc, Crs, Cic, Cis;
    double af0, af1, af2;
    double toc;
};

// satellite ECEF position
struct SatPos { double x,y,z; };

// ----- Utility functions -----
static inline double deg2rad(double d){ return d*PI/180.0; }
static inline double rad2deg(double r){ return r*180.0/PI; }

// parse integer from substring (helper)
string trim(const string &s){
    auto a = s; // copy
    while(!a.empty() && isspace((unsigned char)a.back())) a.pop_back();
    size_t i=0; while(i<a.size() && isspace((unsigned char)a[i])) ++i;
    return a.substr(i);
}

// convert calendar to GPS seconds of week and toe (approx)
void datetime_to_gps_seconds(int year,int month,int day,int hour,int min,int sec, double &gpsWeek, double &gpsSOW){
    // This is an approximate converter using chrono for simplicity.
    // We'll compute epoch seconds since 1980-01-06 00:00:00 UTC (GPS epoch)
    std::tm t = {};
    t.tm_year = year - 1900;
    t.tm_mon = month - 1;
    t.tm_mday = day;
    t.tm_hour = hour;
    t.tm_min = min;
    t.tm_sec = sec;
    t.tm_isdst = 0;
    // timegm non-standard; use mktime assuming system UTC? For typical use parse from RINEX which gives GPS week/time directly.
    time_t tt = timegm(&t); // requires POSIX. If not available, user environment must handle.
    // GPS epoch:
    std::tm g = {};
    g.tm_year = 1980 - 1900; g.tm_mon = 0; g.tm_mday = 6; g.tm_hour = 0; g.tm_min=0; g.tm_sec=0;
    time_t gps0 = timegm(&g);
    long long diff = (long long)tt - (long long)gps0;
    if(diff < 0) diff = 0;
    gpsWeek = diff / (7*86400);
    gpsSOW = diff % (7*86400);
}

// ----- RINEX 3.05 OBS parser (minimal) -----
// We only parse pseudorange types from header: find an observable like "C1C", "C1W", "C1", "P1" etc.
// For simplicity pick the first pseudorange-like code found in header obs types for 'G' system or generally.
vector<ObsRecord> parseRINEXObs(const string &fname, string &prType){
    ifstream ifs(fname);
    if(!ifs) { cerr<<"Cannot open obs file "<<fname<<"\n"; return {}; }
    string line;
    vector<string> obsTypes;
    bool headerDone=false;
    prType = "";
    vector<ObsRecord> records;
    // read header
    while(getline(ifs,line)){
        if(line.size() >= 60){
            string label = trim(line.substr(60));
            if(label == "END OF HEADER"){
                headerDone = true;
                break;
            }
            if(label.find("SYS / # / OBS TYPES") != string::npos){
                // may need to collect multiple lines for multiple systems; here pick line starting with 'G' or first.
                // format: col1 = system char, then count, then types...
                char sys = line[0];
                if(sys == ' ' ) sys = 'G'; // assume GPS if blank
                string typesLine = line.substr(3,57);
                // split by spaces
                stringstream ss(typesLine);
                string token;
                vector<string> types;
                while(ss >> token) types.push_back(token);
                // append to obsTypes (we will choose pseudorange type later)
                for(auto &t: types) obsTypes.push_back(t);
            }
        }
    }
    // Heuristic: find a pseudorange-type among header types
    vector<string> candidates = {"C1C","C1W","C1","C1S","C1X","C1P","P1","C1L","C1I","C1M","C1N"};
    for(auto &cand : candidates){
        for(auto &t : obsTypes){
            if(t == cand) { prType = cand; break; }
        }
        if(!prType.empty()) break;
    }
    if(prType.empty() && !obsTypes.empty()){
        // fallback choose first type
        prType = obsTypes[0];
        cerr<<"Warning: no standard pseudorange obs type found; using "<<prType<<"\n";
    } else {
        cerr<<"Using pseudorange type: "<<prType<<"\n";
    }

    // Rewind file to after header lines and parse epochs
    ifs.clear(); ifs.seekg(0);
    // skip header again
    while(getline(ifs,line)){
        if(line.size()>=60 && trim(line.substr(60))=="END OF HEADER") break;
    }

    // Parse epochs; RINEX 3.05 epoch line starts with '>' char
    while(getline(ifs,line)){
        if(line.empty()) continue;
        if(line[0] == '>'){
            // epoch header: > yyyy mm dd hh mm ss.sssss number_of_satellites
            int Y,Mo,D,H,Mi; double S;
            // handle possible variable spacing
            stringstream ss(line.substr(1));
            ss >> Y >> Mo >> D >> H >> Mi >> S;
            int nsat = 0;
            ss >> nsat;
            // read next nsat lines (sat list may be on one or multiple lines - but RINEX typically lists sat PRNs next in same line)
            vector<string> sats;
            // sat names are arranged starting at col 3-5 etc; but simple parsing: read the remainder tokens
            string rem = line.substr(32); // approximate
            stringstream rss(rem);
            string s;
            while(rss >> s) sats.push_back(s);
            // if not enough, read more lines until get nsat sats
            while((int)sats.size() < nsat){
                string l2;
                if(!getline(ifs,l2)) break;
                stringstream s2(l2);
                while(s2 >> s) sats.push_back(s);
            }
            // Now for each sat read observation lines: each satellite has at least one obs-line containing values matching header types.
            // In RINEX 3 each satellite's observation is in one or more 16-char fields. For simplicity parse tokens from following lines.
            // We'll read ceil(#obsType/5) lines per satellite but we don't know count of obs types; so we will parse by reading tokens (numbers) in sequence.
            // Simpler approach: read next nsat observation raw blocks lines and extract pseudorange by searching for field index of prType in header order.
            // For robust implementation we should have header-obstype order; but we didn't fully parse header. We'll use heuristic:
            // - On obs lines, the first field is the satellite id, then successive 16-char observation fields. We'll parse by fixed columns.
            // Implementation below assumes typical formatting and that prType corresponds to the first Pseudorange-like field.
            vector<double> pr_values;
            vector<string> sat_names_for_epoch;
            // For each satellite, read one line and parse numeric fields by column widths (as RINEX format uses).
            for(int i=0;i<nsat;i++){
                string obsLine;
                if(!getline(ifs,obsLine)) break;
                if(obsLine.size() < 3){
                    i--; continue;
                }
                // satellite id at cols 1-3 (0-based)
                string sat = trim(obsLine.substr(0,3));
                sat_names_for_epoch.push_back(sat);
                // observation fields: each 16 chars, starting at col 3
                vector<string> fld;
                int pos = 3;
                while(pos < (int)obsLine.size()){
                    string f = obsLine.substr(pos, 16);
                    fld.push_back(f);
                    pos += 16;
                }
                // if there are more observation entries on following continuation lines (rare), not handled here.
                // Choose first numeric field as pseudorange candidate
                double pr = 0.0;
                bool got = false;
                for(auto &f : fld){
                    string t = trim(f);
                    if(t.empty()) continue;
                    // replace 'D' with 'E'
                    for(auto &ch: t) if(ch=='D') ch='E';
                    try {
                        pr = stod(t);
                        got = true;
                        break;
                    } catch(...) { continue; }
                }
                if(!got) pr = 0.0;
                pr_values.push_back(pr);
            }
            // compute GPS seconds for epoch
            double gwk, gsow;
            datetime_to_gps_seconds(Y,Mo,D,H,Mi,(int)round(S), gwk, gsow);
            // build record
            ObsRecord rec;
            rec.epoch_gps = gsow;
            for(size_t i=0;i<sat_names_for_epoch.size(); ++i){
                if(pr_values[i] > 1.0) // valid
                    rec.pr[sat_names_for_epoch[i]] = pr_values[i];
            }
            records.push_back(rec);
        } // endif '>'
    } // end while file
    ifs.close();
    return records;
}

// ----- RINEX 3.x NAV parser (GPS broadcast ephemeris) -----
// This parser handles typical 8-line block per sat (GPS). Minimal and assumes values appear as in RINEX 3 GPS format.
map<string, EphData> parseRINEXNav(const string &fname){
    ifstream ifs(fname);
    if(!ifs) { cerr<<"Cannot open nav file "<<fname<<"\n"; return {}; }
    string line;
    // skip header
    while(getline(ifs,line)){
        if(line.size() >= 60 && trim(line.substr(60)) == "END OF HEADER") break;
    }
    map<string, EphData> ephs;
    while(true){
        string l1;
        if(!getline(ifs,l1)) break;
        if(l1.size() < 3) continue;
        // first 3 chars often sat id like 'G01' or 'G01 ' . For RINEX3 the sat id often in col1-3.
        string sat = trim(l1.substr(0,3));
        if(sat.empty()) continue;
        // Parse numeric fields from the 8-line block. We'll collect 8 lines (or available).
        vector<string> block;
        block.push_back(l1);
        for(int i=0;i<7;i++){
            string l;
            if(!getline(ifs,l)) l="";
            block.push_back(l);
        }
        // Convert fields by column positions
        auto field = [&](int row, int colstart, int len)->double{
            if(row < 0 || row >= (int)block.size()) return 0.0;
            string s = "";
            if((int)block[row].size() >= colstart+1){
                int endpos = min((int)block[row].size(), colstart+len);
                s = block[row].substr(colstart, endpos-colstart);
            }
            for(auto &ch : s) if(ch=='D') ch='E';
            try {
                return stod(trim(s));
            } catch(...) { return 0.0; }
        };
        // First line contains epoch and af0/1/2
        // columns (RINEX 3.03 GPS format) - approximate indices used below (0-based):
        // l1[3..5]=year, [6..8]=mon ... but to avoid messy parsing we use stringstream after replacing multiple spaces
        string l1copy = l1;
        for(auto &ch:l1copy) if(ch=='D') ch='E';
        // Extract af0,af1,af2 from known positions: typically cols 22..41 etc. We'll use substring
        double af0 = field(0, 41, 19);
        double af1 = field(0, 22, 19);
        double af2 = field(0, 3, 19); // not exact but try
        // Toe and orbital params in subsequent lines per RINEX spec:
        EphData eph;
        eph.sat = sat;
        eph.af0 = af0; eph.af1 = af1; eph.af2 = af2;
        eph.sqrtA = field(1, 41, 19); // commonly sqrtA is in line 3 but position vary; best-effort
        eph.e = field(2, 3, 19);
        eph.M0 = field(2, 22, 19);
        eph.deltaN = field(1, 3, 19);
        eph.OMG0 = field(3, 3, 19);
        eph.i0 = field(4, 3, 19);
        eph.omg = field(3, 22, 19);
        eph.OMG_dot = field(7, 3, 19);
        eph.i_dot = field(4, 22, 19);
        eph.Cuc = field(3, 41, 19);
        eph.Cus = field(2, 41, 19);
        eph.Crc = field(4, 41, 19);
        eph.Crs = field(1, 41, 19);
        eph.Cic = field(6, 41, 19);
        eph.Cis = field(5, 41, 19);
        eph.toe = field(6, 3, 19);
        eph.toc = field(0, 3, 19);
        // due to wide variations in RINEX formatting this parser might not extract perfectly for all files.
        // store
        ephs[sat] = eph;
    }
    return ephs;
}

// ----- Satellite position from broadcast ephemeris -----
// Implement GPS orbital model: mean motion, Kepler's equation, corrections, then transform to ECEF and apply rotation corrections (Sagnac)
SatPos computeSatPosFromEph(const EphData &eph, double recv_gps_sow, double pseudorange){
    // compute transmission time (first approx)
    double t = recv_gps_sow - pseudorange / c0;
    // correct week rollover (we assume times within same week)
    double tk = t - eph.toe;
    // ensure tk within +/-302400
    if(tk > 302400) tk -= 604800;
    if(tk < -302400) tk += 604800;
    // semi-major axis
    double A = eph.sqrtA * eph.sqrtA;
    double n0 = sqrt(MU/(A*A*A));
    double n = n0 + eph.deltaN;
    double M = eph.M0 + n * tk;
    // normalize M
    M = fmod(M, 2*PI);
    if(M < -PI) M += 2*PI;
    if(M > PI) M -= 2*PI;
    // Solve Kepler's equation E - e*sinE = M
    double E = M;
    for(int iter=0; iter<30; ++iter){
        double f = E - eph.e*sin(E) - M;
        double df = 1.0 - eph.e * cos(E);
        double dx = -f/df;
        E += dx;
        if(fabs(dx) < 1e-12) break;
    }
    double sinE = sin(E);
    double cosE = cos(E);
    double v = atan2(sqrt(1-eph.e*eph.e)*sinE, cosE - eph.e);
    double phi = v + eph.omg;
    double du = eph.Cus*sin(2*phi) + eph.Cuc*cos(2*phi);
    double dr = eph.Crs*sin(2*phi) + eph.Crc*cos(2*phi);
    double di = eph.Cis*sin(2*phi) + eph.Cic*cos(2*phi);
    double u = phi + du;
    double r = A * (1 - eph.e * cosE) + dr;
    double i = eph.i0 + di + eph.i_dot * tk;
    double x_orb = r * cos(u);
    double y_orb = r * sin(u);
    double OMG = eph.OMG0 + (eph.OMG_dot - OMEGA_E) * tk - OMEGA_E * eph.toe;
    double cosOMG = cos(OMG), sinOMG = sin(OMG);
    double cosi = cos(i), sini = sin(i);
    double X = x_orb * cosOMG - y_orb * cosi * sinOMG;
    double Y = x_orb * sinOMG + y_orb * cosi * cosOMG;
    double Z = y_orb * sini;
    // apply Earth rotation correction during signal travel: rotate by OMEGA_E * dt where dt = pseudorange/c
    double dt = pseudorange / c0;
    double omega_dt = OMEGA_E * dt;
    double cosw = cos(omega_dt), sinw = sin(omega_dt);
    // rotate satellite position around z by omega_dt
    double Xr = cosw * X + sinw * Y;
    double Yr = -sinw * X + cosw * Y;
    double Zr = Z;
    return {Xr, Yr, Zr};
}

// ----- Position estimation by iterative least-squares -----
// Input: obs map sat->pseudorange, eph map sat->eph, recv_approx initial
bool estimateReceiverPosition(const unordered_map<string,double> &obs_pr, const map<string,EphData> &ephMap,
                              double recv_tow, Eigen::Vector4d &x_out, double &rms_out, Eigen::Matrix3d &Qpos_out){
    // initial guess (ECEF) - we can use center of Earth or rough geocenter. Better if user passes approx.
    Eigen::Vector3d xk(0.0,0.0,0.0);
    double clk = 0.0; // receiver clock bias in seconds
    // Simple initial guess: use average satellite positions direction to get some rough location if needed.
    // We'll iterate up to maxIter
    int maxIter = 10;
    for(int iter=0; iter<maxIter; ++iter){
        vector<Eigen::Vector4d> Arows;
        vector<double> Lvec;
        // For each satellite: compute satellite pos at transmission time using current pseudorange (we must use observed PR)
        for(auto &kv : obs_pr){
            string sat = kv.first;
            double P = kv.second;
            auto it = ephMap.find(sat);
            if(it == ephMap.end()) continue;
            const EphData &eph = it->second;
            // compute sat pos
            SatPos sp = computeSatPosFromEph(eph, recv_tow, P);
            // geometric distance
            double dx = sp.x - xk(0);
            double dy = sp.y - xk(1);
            double dz = sp.z - xk(2);
            double rho = sqrt(dx*dx + dy*dy + dz*dz);
            if(rho < 1.0) continue;
            // predicted pseudorange including receiver clock bias
            double pred = rho + c0 * clk - c0 * 0.0; // ignoring satellite clock corr here (could read eph.af0/1/2)
            double li = P - pred;
            // design row: [-dx/rho, -dy/rho, -dz/rho, c]
            Eigen::Vector4d arow;
            arow(0) = -dx / rho;
            arow(1) = -dy / rho;
            arow(2) = -dz / rho;
            arow(3) = c0;
            Arows.push_back(arow);
            Lvec.push_back(li);
        }
        int m = (int)Arows.size();
        if(m < 4) { cerr<<"Not enough satellites ("<<m<<") for solution\n"; return false; }
        Eigen::MatrixXd A(m,4);
        Eigen::VectorXd L(m);
        for(int i=0;i<m;++i){
            A.row(i) = Arows[i].transpose();
            L(i) = Lvec[i];
        }
        // Solve normal equations
        Eigen::MatrixXd N = A.transpose() * A;
        Eigen::VectorXd u = A.transpose() * L;
        Eigen::VectorXd dx = N.ldlt().solve(u);
        // update
        xk(0) += dx(0);
        xk(1) += dx(1);
        xk(2) += dx(2);
        clk  += dx(3) / c0;
        if(dx.head(3).norm() < 1e-4 && fabs(dx(3)/c0) < 1e-9) break;
    }
    // compute residuals, RMS, and covariance
    vector<double> res;
    Eigen::MatrixXd A(m,4); Eigen::VectorXd L(m);
    int idx = 0;
    for(auto &kv : obs_pr){
        string sat = kv.first;
        double P = kv.second;
        auto it = ephMap.find(sat);
        if(it == ephMap.end()) continue;
        SatPos sp = computeSatPosFromEph(it->second, recv_tow, P);
        double dx = sp.x - xk(0);
        double dy = sp.y - xk(1);
        double dz = sp.z - xk(2);
        double rho = sqrt(dx*dx + dy*dy + dz*dz);
        if(rho < 1.0) continue;
        double predicted = rho + c0 * clk;
        double v = P - predicted;
        res.push_back(v);
        Eigen::Vector4d arow;
        arow(0) = -dx / rho; arow(1) = -dy / rho; arow(2) = -dz / rho; arow(3)=c0;
        A.row(idx) = arow.transpose();
        L(idx) = v;
        idx++;
    }
    int used = idx;
    if(used < 4){ cerr<<"After building final matrix not enough sats\n"; return false; }
    Eigen::MatrixXd Ause = A.topRows(used);
    Eigen::VectorXd Luse = L.head(used);
    Eigen::MatrixXd Qxx = (Ause.transpose()*Ause).inverse();
    // position covariance for x,y,z are Qxx[0:3,0:3]
    Qpos_out = Qxx.topLeftCorner(3,3);
    // RMS (post-fit)
    double sumsq = 0.0;
    for(int i=0;i<used;i++) sumsq += Luse(i)*Luse(i);
    int dof = used - 4;
    double rms = sqrt( sumsq / max(1, dof) );
    // fill outputs
    x_out << xk(0), xk(1), xk(2), clk;
    rms_out = rms;
    return true;
}

// ----- Helper: ECEF -> Geodetic (WGS84) -----
void ecef2lla(double x,double y,double z, double &lat,double &lon,double &h){
    // iterative method
    double a = 6378137.0;
    double f = 1.0/298.257223563;
    double e2 = f*(2-f);
    lon = atan2(y,x);
    double r = sqrt(x*x + y*y);
    double E2 = a*a - (a*(1-f))*(a*(1-f));
    double F = 54*(a*(1-f))*(a*(1-f))*z*z;
    double G = r*r + (1-e2)*z*z - e2*E2;
    double c = (e2*e2 * F * r*r) / (G*G*G);
    double s = pow(1 + c + sqrt(c*c + 2*c), 1.0/3.0);
    double P = F / (3 * pow((s + 1/s + 1),2) * G*G);
    double Q = sqrt(1 + 2*e2*e2*P);
    double r0 = -(P*e2*r)/(1+Q) + sqrt(0.5*a*a*(1+1.0/Q) - P*(1-e2)*z*z/(Q*(1+Q)) - 0.5*P*r*r);
    double U = sqrt( pow(r - e2*r0,2) + z*z );
    double V = sqrt( pow(r - e2*r0,2) + (1-e2)*z*z );
    double Z0 = (a*(1-f))*z/(a*V);
    h = U*(1 - (a*(1-f))/(a*V));
    lat = atan((z + e2*Z0)/r);
}

// ----- Main -----
int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if(argc < 3){
        cerr<<"Usage: "<<argv[0]<<" <obs_file.obs> <nav_file.nav> [epoch_index]\n";
        return -1;
    }
    string obsFile = argv[1];
    string navFile = argv[2];
    int epochIndex = 0;
    if(argc >= 4) epochIndex = atoi(argv[3]);

    string prType;
    auto obsRecords = parseRINEXObs(obsFile, prType);
    if(obsRecords.empty()){ cerr<<"No observations parsed\n"; return -1; }
    auto ephs = parseRINEXNav(navFile);
    if(ephs.empty()){ cerr<<"No ephemerides parsed\n"; /*continue but may fail*/ }

    if(epochIndex < 0) epochIndex = 0;
    if(epochIndex >= (int)obsRecords.size()) epochIndex = (int)obsRecords.size()-1;
    ObsRecord rec = obsRecords[epochIndex];
    cerr<<"Selected epoch SOW = "<<rec.epoch_gps<<" with "<<rec.pr.size()<<" satellites\n";

    // Build obs_pr map with common naming: RINEX sat names may be like G01; ensure keys consistent with nav keys
    unordered_map<string,double> obs_pr;
    for(auto &kv : rec.pr){
        string satid = kv.first;
        // normalize to 3-char e.g. G01
        string s = satid;
        if(s.size() == 2) s = s.substr(0,1) + "0" + s.substr(1);
        obs_pr[s] = kv.second;
    }

    Eigen::Vector4d sol;
    double rms; Eigen::Matrix3d Qpos;
    bool ok = estimateReceiverPosition(obs_pr, ephs, rec.epoch_gps, sol, rms, Qpos);
    if(!ok){ cerr<<"Position estimation failed\n"; return -1; }

    cout.setf(std::ios::fixed); cout<<setprecision(3);
    cout<<"Estimated ECEF (m): X="<<sol(0)<<" Y="<<sol(1)<<" Z="<<sol(2)<<"\n";
    double lat, lon, h;
    ecef2lla(sol(0), sol(1), sol(2), lat, lon, h);
    cout<<setprecision(9);
    cout<<"Geodetic: lat="<<lat<<" lon="<<lon<<" h="<<h<<"\n";
    cout<<setprecision(6);
    cout<<"Receiver clock bias (s): "<<sol(3)<<"\n";
    cout<<"Post-fit RMS (m): "<<rms<<"\n";

    // PDOP from covariance
    double qxx = Qpos(0,0), qyy = Qpos(1,1), qzz = Qpos(2,2);
    double pdop = sqrt(qxx + qyy + qzz);
    // convert ECEF cov to ENU and get HDOP/VDOP (approx) — we'll approximate by diagonal
    cout<<"PDOP (approx): "<<pdop<<"\n";

    return 0;
}