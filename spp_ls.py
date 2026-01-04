import argparse
import math
import re
import datetime as dt
import numpy as np

# ---------------- constants ----------------
C = 299792458.0
MU_GPS = 3.986005e14
OMEGA_E = 7.2921151467e-5
F_REL = -4.442807633e-10

A_WGS84 = 6378137.0
F_WGS84 = 1 / 298.257223563
E2_WGS84 = F_WGS84 * (2 - F_WGS84)

GPS_EPOCH = dt.datetime(1980, 1, 6, 0, 0, 0)


# ---------------- time utils ----------------
def gps_week_tow(t: dt.datetime):
    sec = (t - GPS_EPOCH).total_seconds()
    week = int(sec // 604800)
    tow = sec - week * 604800
    return week, tow


def adjust_week_time(dt_sec: float):
    # wrap to [-302400, 302400]
    if dt_sec > 302400:
        dt_sec -= 604800
    elif dt_sec < -302400:
        dt_sec += 604800
    return dt_sec


# ---------------- coord utils ----------------
def ecef_to_llh(xyz):
    x, y, z = xyz
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - E2_WGS84))
    for _ in range(10):
        N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat_new = math.atan2(z, p * (1 - E2_WGS84 * N / (N + h)))
        if abs(lat_new - lat) < 1e-12:
            lat = lat_new
            break
        lat = lat_new
    N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - N
    return lat, lon, h


def ecef_to_enu_matrix(lat, lon):
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    # ECEF -> ENU
    return np.array([
        [-slon,        clon,       0.0],
        [-slat * clon, -slat * slon, clat],
        [ clat * clon,  clat * slon, slat],
    ])


def elevation_angle(rx_xyz, sat_xyz):
    lat, lon, _ = ecef_to_llh(rx_xyz)
    R = ecef_to_enu_matrix(lat, lon)
    enu = R @ (sat_xyz - rx_xyz)
    e, n, u = enu
    horiz = math.hypot(e, n)
    return math.atan2(u, horiz)


# ---------------- troposphere (simple) ----------------
def simple_tropo_delay(el_rad):
    # very simple: 2.3m zenith hydrostatic delay with mapping ~ 1/sin(el)
    sinel = max(math.sin(el_rad), 0.1)
    return 2.3 / sinel


# ---------------- RINEX parsing ----------------
def parse_float(field: str):
    s = field.replace("D", "E").replace("d", "E").strip()
    if not s:
        return np.nan
    return float(s)


def parse_rinex3_obs_header(path):
    obs_types = {}
    approx_xyz = None
    ant_hen = (0.0, 0.0, 0.0)  # H/E/N
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "APPROX POSITION XYZ" in line:
                approx_xyz = (
                    float(line[0:14]),
                    float(line[14:28]),
                    float(line[28:42]),
                )
            if "ANTENNA: DELTA H/E/N" in line:
                ant_hen = (
                    float(line[0:14]),
                    float(line[14:28]),
                    float(line[28:42]),
                )
            if "SYS / # / OBS TYPES" in line:
                sys = line[0]
                n = int(line[3:6])
                types = line[7:60].split()
                while len(types) < n:
                    line2 = next(f)
                    types += line2[7:60].split()
                obs_types[sys] = types[:n]
            if "END OF HEADER" in line:
                break
    return obs_types, approx_xyz, ant_hen


class ObsEpoch:
    def __init__(self, t, sats):
        self.t = t
        self.sats = sats  # dict: sat -> dict(obsType -> float)


def read_rinex3_obs(path, keep_system="G"):
    obs_types, approx_xyz, ant_hen = parse_rinex3_obs_header(path)

    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "END OF HEADER" in line:
                break

        epochs = []
        for line in f:
            if not line:
                break
            if line.startswith(">"):
                year = int(line[2:6])
                mon = int(line[7:9])
                day = int(line[10:12])
                hour = int(line[13:15])
                minute = int(line[16:18])
                sec = float(line[19:29])
                nsat = int(line[32:35])

                micro = int(round((sec - int(sec)) * 1e6))
                t = dt.datetime(year, mon, day, hour, minute, int(sec), micro)

                sats = {}
                for _ in range(nsat):
                    sline = next(f)
                    sat = sline[0:3].strip()
                    sys = sat[0]
                    if sys != keep_system:
                        continue

                    types = obs_types[sys]
                    ntypes = len(types)

                    fields = sline[3:].rstrip("\n")
                    while len(fields) < 16 * ntypes:
                        fields += next(f)[3:].rstrip("\n")

                    vals = []
                    for i in range(ntypes):
                        chunk = fields[i * 16:(i + 1) * 16]
                        vstr = chunk[0:14]
                        try:
                            v = float(vstr) if vstr.strip() else np.nan
                        except:
                            v = np.nan
                        vals.append(v)

                    sats[sat] = {types[i]: vals[i] for i in range(ntypes)}

                epochs.append(ObsEpoch(t, sats))

    return epochs, approx_xyz, ant_hen


def parse_rinex3_nav_gps(path):
    # only parse GPS (G) ephemeris blocks (8 lines per record)
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "END OF HEADER" in line:
                break
        body = f.readlines()

    ephs = {}
    i = 0
    while i < len(body):
        line = body[i]
        if not line.strip():
            i += 1
            continue

        sys = line[0]
        if sys == "R":
            rec_len = 4
        elif sys in ("G", "E", "C", "J", "I", "S"):
            rec_len = 8
        else:
            i += 1
            continue

        rec = body[i:i + rec_len]
        i += rec_len

        if sys != "G":
            continue

        prn = rec[0][0:3].strip()
        year = int(rec[0][4:8])
        mon = int(rec[0][9:11])
        day = int(rec[0][12:14])
        hour = int(rec[0][15:17])
        minute = int(rec[0][18:20])
        sec = int(rec[0][21:23])
        toc = dt.datetime(year, mon, day, hour, minute, sec)

        af0 = parse_float(rec[0][23:42])
        af1 = parse_float(rec[0][42:61])
        af2 = parse_float(rec[0][61:80])

        IODE = parse_float(rec[1][4:23])
        Crs = parse_float(rec[1][23:42])
        dn = parse_float(rec[1][42:61])
        M0 = parse_float(rec[1][61:80])

        Cuc = parse_float(rec[2][4:23])
        e = parse_float(rec[2][23:42])
        Cus = parse_float(rec[2][42:61])
        sqrtA = parse_float(rec[2][61:80])

        Toe = parse_float(rec[3][4:23])
        Cic = parse_float(rec[3][23:42])
        Omega0 = parse_float(rec[3][42:61])
        Cis = parse_float(rec[3][61:80])

        i0 = parse_float(rec[4][4:23])
        Crc = parse_float(rec[4][23:42])
        omega = parse_float(rec[4][42:61])
        Omegadot = parse_float(rec[4][61:80])

        IDOT = parse_float(rec[5][4:23])
        GPSWeek = parse_float(rec[5][42:61])

        TGD = parse_float(rec[6][42:61])

        eph = dict(
            prn=prn, toc=toc,
            af0=af0, af1=af1, af2=af2,
            IODE=IODE, Crs=Crs, dn=dn, M0=M0,
            Cuc=Cuc, e=e, Cus=Cus, sqrtA=sqrtA,
            Toe=Toe, Cic=Cic, Omega0=Omega0, Cis=Cis,
            i0=i0, Crc=Crc, omega=omega, Omegadot=Omegadot,
            IDOT=IDOT, GPSWeek=GPSWeek, TGD=TGD,
        )
        ephs.setdefault(prn, []).append(eph)

    for prn in ephs:
        ephs[prn] = sorted(ephs[prn], key=lambda x: x["toc"])
    return ephs


# ---------------- satellite model ----------------
def select_eph(eph_list, week, tow):
    best, best_abs = None, 1e18
    for e in eph_list:
        dt_sec = adjust_week_time(tow - e["Toe"])
        ad = abs(dt_sec)
        if ad < best_abs:
            best_abs = ad
            best = e
    return best, best_abs


def sat_pos_and_clock_gps(eph, week, tow):
    # satellite position at tow (GPST), and clock correction (s) including relativity and TGD
    A = eph["sqrtA"] ** 2
    n0 = math.sqrt(MU_GPS / A ** 3)
    n = n0 + eph["dn"]

    tk = adjust_week_time(tow - eph["Toe"])
    M = eph["M0"] + n * tk
    ecc = eph["e"]

    # Kepler
    E = M
    for _ in range(20):
        E_new = M + ecc * math.sin(E)
        if abs(E_new - E) < 1e-12:
            E = E_new
            break
        E = E_new

    sinE, cosE = math.sin(E), math.cos(E)
    dtr = F_REL * ecc * eph["sqrtA"] * sinE

    v = math.atan2(math.sqrt(1 - ecc * ecc) * sinE, cosE - ecc)
    phi = v + eph["omega"]
    sin2, cos2 = math.sin(2 * phi), math.cos(2 * phi)

    u = phi + eph["Cuc"] * cos2 + eph["Cus"] * sin2
    r = A * (1 - ecc * cosE) + eph["Crc"] * cos2 + eph["Crs"] * sin2
    i = eph["i0"] + eph["IDOT"] * tk + eph["Cic"] * cos2 + eph["Cis"] * sin2

    x_p = r * math.cos(u)
    y_p = r * math.sin(u)

    Omega = eph["Omega0"] + (eph["Omegadot"] - OMEGA_E) * tk - OMEGA_E * eph["Toe"]

    cosO, sinO = math.cos(Omega), math.sin(Omega)
    cosi, sini = math.cos(i), math.sin(i)

    x = x_p * cosO - y_p * cosi * sinO
    y = x_p * sinO + y_p * cosi * cosO
    z = y_p * sini

    # clock
    toc_week, toc_tow = gps_week_tow(eph["toc"])
    dt_sec = (week - toc_week) * 604800 + (tow - toc_tow)
    dt_sec = adjust_week_time(dt_sec)
    dt_sv = eph["af0"] + eph["af1"] * dt_sec + eph["af2"] * dt_sec * dt_sec + dtr - eph["TGD"]

    return np.array([x, y, z], dtype=float), dt_sv


# ---------------- SPP solver (WLS) ----------------
def solve_epoch_spp(epoch: ObsEpoch, ephs, x0, dt0, obs_code="C1C", max_iter=10):
    week, tow_rx = gps_week_tow(epoch.t)

    sats, P_list = [], []
    for sat, od in epoch.sats.items():
        P = od.get(obs_code, np.nan)
        if np.isnan(P) or P <= 1e5:
            continue
        if sat not in ephs:
            continue
        sats.append(sat)
        P_list.append(P)

    if len(sats) < 4:
        return None

    x = np.array(x0, dtype=float)
    dt_r = float(dt0)

    for _ in range(max_iter):
        A_rows, y_rows, w_rows = [], [], []

        for sat, P in zip(sats, P_list):
            eph, absdt = select_eph(ephs[sat], week, tow_rx)
            if eph is None or absdt > 4 * 3600:
                continue

            # transmit time iteration (2 steps)
            tau = P / C
            t_tx = (tow_rx - tau) % 604800.0

            for __ in range(2):
                sat_pos, dt_sv = sat_pos_and_clock_gps(eph, week, t_tx)
                rho_corr = P - C * (dt_r - dt_sv)
                tau = rho_corr / C
                t_tx = (tow_rx - tau) % 604800.0

            sat_pos, dt_sv = sat_pos_and_clock_gps(eph, week, t_tx)

            # Earth rotation correction
            ang = OMEGA_E * tau
            ca, sa = math.cos(ang), math.sin(ang)
            sat_rot = np.array([ca * sat_pos[0] + sa * sat_pos[1],
                                -sa * sat_pos[0] + ca * sat_pos[1],
                                sat_pos[2]], dtype=float)

            rho = np.linalg.norm(sat_rot - x)
            el = elevation_angle(x, sat_rot)
            trop = simple_tropo_delay(el)

            P_hat = rho + C * (dt_r - dt_sv) + trop
            v = P - P_hat

            ux, uy, uz = (x - sat_rot) / rho
            A_rows.append([ux, uy, uz, C])
            y_rows.append(v)

            sinel = max(math.sin(el), 0.1)
            w_rows.append(sinel * sinel)  # weight ~ sin^2(el)

        if len(A_rows) < 4:
            return None

        A = np.array(A_rows, dtype=float)
        y = np.array(y_rows, dtype=float)
        w = np.array(w_rows, dtype=float)
        sw = np.sqrt(w)

        Aw = A * sw[:, None]
        yw = y * sw

        dx, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
        x = x + dx[:3]
        dt_r = dt_r + dx[3]

        if np.linalg.norm(dx[:3]) < 1e-4 and abs(dx[3]) < 1e-10:
            break

    # post-fit RMS
    rms = float(np.sqrt(np.mean(y * y))) if len(y_rows) else np.nan
    return x, dt_r, rms, len(A_rows)


def apply_antenna_offset_marker(x_apc, ant_hen):
    # ant_hen is (H, E, N) from marker to antenna phase center
    H, E, N = ant_hen
    lat, lon, _ = ecef_to_llh(x_apc)
    R = ecef_to_enu_matrix(lat, lon)      # ECEF->ENU
    shift_ecef = R.T @ np.array([E, N, H]) # ENU->ECEF
    return x_apc - shift_ecef


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("obs", help="RINEX3 observation file (.25o)")
    ap.add_argument("nav", help="RINEX3 nav file (.25p)")
    ap.add_argument("--out", default="spp_out.csv", help="output CSV")
    ap.add_argument("--obs_code", default="C1C", help="pseudorange code for GPS (default C1C)")
    ap.add_argument("--marker", action="store_true", help="output marker coord (apply ANTENNA: DELTA H/E/N)")
    args = ap.parse_args()

    epochs, approx_xyz, ant_hen = read_rinex3_obs(args.obs, keep_system="G")
    ephs = parse_rinex3_nav_gps(args.nav)

    if approx_xyz is None:
        raise RuntimeError("No APPROX POSITION XYZ in obs header")

    x_init = np.array(approx_xyz, dtype=float)
    dt_init = 0.0

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("time,x,y,z,lat_deg,lon_deg,h_m,clk_bias_s,rms_m,nsat\n")
        xs = []
        for ep in epochs:
            sol = solve_epoch_spp(ep, ephs, x_init, dt_init, obs_code=args.obs_code)
            if sol is None:
                continue
            x, dtr, rms, ns = sol

            # warm start next epoch
            x_init = x
            dt_init = dtr

            x_out = apply_antenna_offset_marker(x, ant_hen) if args.marker else x
            lat, lon, h = ecef_to_llh(x_out)
            f.write(f"{ep.t.isoformat()},{x_out[0]:.4f},{x_out[1]:.4f},{x_out[2]:.4f},"
                    f"{math.degrees(lat):.10f},{math.degrees(lon):.10f},{h:.4f},"
                    f"{dtr:.12e},{rms:.4f},{ns}\n")
            xs.append(x_out)

    xs = np.array(xs)
    mean_xyz = xs.mean(axis=0)
    lat, lon, h = ecef_to_llh(mean_xyz)

    print("=== SPP done ===")
    print(f"epochs solved: {len(xs)} / {len(epochs)}")
    print("mean ECEF (m):", mean_xyz)
    print("mean BLH:", math.degrees(lat), math.degrees(lon), h)
    print(f"output saved to: {args.out}")


if __name__ == "__main__":
    main()
