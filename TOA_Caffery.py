import numpy as np
import matplotlib.pyplot as plt

def caffery_toa_localization(anchors, distances):
    """
    anchors: N×2 array, anchor positions
    distances: N array, measured distances (TOA*c)
    """
    a1 = anchors[0]
    d1 = distances[0]
    A, b = [], []
    for i in range(1, len(anchors)):
        xi, yi = anchors[i]
        x1, y1 = a1
        di = distances[i]
        row = [2*(xi-x1), 2*(yi-y1)]
        rhs = (xi**2 - x1**2) + (yi**2 - y1**2) + (d1**2 - di**2)
        A.append(row)
        b.append(rhs)
    A = np.array(A)
    b = np.array(b)
    # 最小二乘解
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    return sol

# 仿真参数
np.random.seed(0)
anchors = np.array([[0,0],[50,0],[0,50],[50,50],[25,80]])
true_pos = np.array([20, 30])
noise_std = 1.0

# 生成测距
d_true = np.linalg.norm(anchors-true_pos,axis=1)
d_meas = d_true + np.random.randn(len(d_true))*noise_std

# Caffery 方法估计
est_pos = caffery_toa_localization(anchors, d_meas)

# 画图
plt.figure(figsize=(6,6))
plt.scatter(anchors[:,0], anchors[:,1], marker='s', label="Anchors")
plt.scatter(true_pos[0], true_pos[1], marker='*', s=200, label="True")
plt.scatter(est_pos[0], est_pos[1], marker='o', label="Caffery Est")
plt.legend()
plt.axis("equal")
plt.grid(True)
plt.title("TOA Localization (Caffery Method)")
plt.show()

print("True position:", true_pos)
print("Estimated (Caffery):", est_pos)
print("Error norm:", np.linalg.norm(est_pos-true_pos))