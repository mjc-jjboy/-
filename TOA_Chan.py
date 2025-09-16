import numpy as np
import matplotlib.pyplot as plt

def chan_toa_localization(anchors, distances, noise_var=1.0):
    """
    Chan TOA localization in 2D
    anchors: N×2 array
    distances: N array (measured distances)
    noise_var: measurement noise variance
    """
    N = len(anchors)
    x1, y1 = anchors[0]
    d1 = distances[0]

    A, b = [], []
    for i in range(1, N):
        xi, yi = anchors[i]
        di = distances[i]
        row = [x1 - xi, y1 - yi]
        rhs = 0.5 * ((x1**2 + y1**2 - d1**2) - (xi**2 + yi**2 - di**2))
        A.append(row)
        b.append(rhs)
    A = np.array(A)
    b = np.array(b)

    # 权重矩阵（这里简单用 I，也可根据噪声方差设置 diag(d_i^2)）
    W = np.eye(len(b)) / noise_var

    # 加权最小二乘解
    est, *_ = np.linalg.lstsq(W @ A, W @ b, rcond=None)
    return est

# ====== 测试 ======
np.random.seed(1)
anchors = np.array([[0,0],[50,0],[0,50],[50,50],[25,80]])
true_pos = np.array([20, 30])
noise_std = 1.0

# 生成 TOA 距离
d_true = np.linalg.norm(anchors - true_pos, axis=1)
d_meas = d_true + np.random.randn(len(d_true)) * noise_std

# Chan 方法估计
est_pos = chan_toa_localization(anchors, d_meas, noise_var=noise_std**2)

# 画图
plt.figure(figsize=(6,6))
plt.scatter(anchors[:,0], anchors[:,1], marker='s', label="Anchors")
plt.scatter(true_pos[0], true_pos[1], marker='*', s=200, label="True")
plt.scatter(est_pos[0], est_pos[1], marker='o', label="Chan Est")
plt.legend()
plt.axis("equal")
plt.grid(True)
plt.title("TOA Localization (Chan Method)")
plt.show()

print("True position:", true_pos)
print("Estimated (Chan):", est_pos)
print("Error norm:", np.linalg.norm(est_pos-true_pos))