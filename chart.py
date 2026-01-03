# Compile and run the three experiments, then compute NEU errors and plot curves and tables.
import subprocess, os, math, pandas as pd, matplotlib.pyplot as plt

base = "/mnt/data"
# compile
subprocess.run(["g++","-O2","-std=c++17",f"{base}/exp1_ls.cpp","-o",f"{base}/exp1_ls"], check=True)
subprocess.run(["g++","-O2","-std=c++17",f"{base}/exp2_static_kf.cpp","-o",f"{base}/exp2_static_kf"], check=True)
subprocess.run(["g++","-O2","-std=c++17",f"{base}/exp3_dyn_kf.cpp","-o",f"{base}/exp3_dyn_kf"], check=True)

# run
subprocess.run([f"{base}/exp1_ls","--obs",f"{base}/RF5Y0470.25o","--nav",f"{base}/BRDC00IGS_R_20250470000_01D_MN.rnx","--sys","C","--out",f"{base}/exp1.txt"], check=True)
subprocess.run([f"{base}/exp2_static_kf","--obs",f"{base}/RF5Y0470.25o","--nav",f"{base}/BRDC00IGS_R_20250470000_01D_MN.rnx","--sys","C","--out",f"{base}/exp2.txt"], check=True)
subprocess.run([f"{base}/exp3_dyn_kf","--obs",f"{base}/9uy70470.25o","--nav",f"{base}/BRDC00IGS_R_20250470000_01D_MN.rnx","--sys","C","--out",f"{base}/exp3.txt"], check=True)

# read approx position from obs header
def read_approx_xyz(path):
    with open(path) as f:
        for line in f:
            if "APPROX POSITION XYZ" in line:
                x,y,z = map(float,line.split()[:3])
                return x,y,z
    raise ValueError("No approx pos")

ref_xyz = read_approx_xyz(f"{base}/RF5Y0470.25o")

# ECEF to NEU
import numpy as np
def ecef2neu(x,y,z,ref):
    X,Y,Z = ref
    r = math.sqrt(X*X+Y*Y)
    lat = math.atan2(Z, r)
    lon = math.atan2(Y, X)
    d = np.array([x-X,y-Y,z-Z])
    R = np.array([[-math.sin(lat)*math.cos(lon), -math.sin(lat)*math.sin(lon), math.cos(lat)],
                  [-math.sin(lon), math.cos(lon), 0],
                  [math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)]])
    n,e,u = R.dot(d)
    return n,e,u

def load_neu(path):
    data=[]
    with open(path) as f:
        for line in f:
            if line.strip()=="" or line.startswith("#"): continue
            parts=line.split()
            x,y,z = map(float, parts[3:6])
            data.append(ecef2neu(x,y,z,ref_xyz))
    return np.array(data)

neu1 = load_neu(f"{base}/exp1.txt")
neu2 = load_neu(f"{base}/exp2.txt")
neu3 = load_neu(f"{base}/exp3.txt")

# plot
def plot_neu(neu, title, fname):
    plt.figure()
    plt.plot(neu[:,0])
    plt.plot(neu[:,1])
    plt.plot(neu[:,2])
    plt.title(title)
    plt.xlabel("Epoch")
    plt.ylabel("Error (m)")
    plt.savefig(f"{base}/{fname}")
    plt.close()

plot_neu(neu1,"LS NEU Error","ls_neu.png")
plot_neu(neu2,"Static KF NEU Error","kf_static_neu.png")
plot_neu(neu3,"Dynamic KF NEU Error","kf_dyn_neu.png")

# convergence table (threshold 1m for 10 consecutive epochs)
def conv_time(neu,thr=1.0,win=10):
    mag = np.linalg.norm(neu,axis=1)
    for i in range(len(mag)-win):
        if np.all(mag[i:i+win]<thr):
            return i
    return None

table = pd.DataFrame({
    "Experiment":["LS","Static KF","Dynamic KF"],
    "ConvEpoch":[conv_time(neu1),conv_time(neu2),conv_time(neu3)],
    "RMS_3D":[np.sqrt(np.mean(np.sum(neu**2,axis=1))) for neu in [neu1,neu2,neu3]]
})
table
