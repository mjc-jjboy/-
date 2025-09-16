import numpy as np
import matplotlib.pyplot as plt


def linear_line_of_position(lines):
  """最小二乘法计算 Linear Line of Position (LLOP)"""
  A, C = [], []
  for (a, b, c) in lines:
    A.append([a, b])
    C.append(-c)
  A, C = np.array(A, dtype=float), np.array(C, dtype=float)
  sol, _, _, _ = np.linalg.lstsq(A, C, rcond=None)
  return sol[0], sol[1]


def generate_lines(true_point, n_lines=5, noise_level=0.5):
  """
  随机生成经过真实点的观测直线，并加噪声

  true_point: (x0, y0) 真实点
  n_lines: 直线条数
  noise_level: 噪声强度
  """
  x0, y0 = true_point
  lines = []

  for _ in range(n_lines):
    # 随机生成直线方向
    angle = np.random.uniform(0, np.pi)
    a = np.cos(angle)
    b = np.sin(angle)

    # 理论直线：a(x-x0)+b(y-y0)=0  => a*x + b*y - (a*x0+b*y0)=0
    c = -(a * x0 + b * y0)

    # 加噪声
    c += np.random.normal(0, noise_level)

    lines.append((a, b, c))

  return lines


def plot_simulation(true_point, lines, estimate):
  """绘制模拟结果"""
  x_vals = np.linspace(true_point[0] - 5, true_point[0] + 5, 400)
  plt.figure(figsize=(7, 7))

  # 绘制观测直线
  for i, (a, b, c) in enumerate(lines):
    if abs(b) > 1e-8:
      y_vals = (-a * x_vals - c) / b
      plt.plot(x_vals, y_vals, label=f"Line {i + 1}")
    else:
      x_const = -c / a
      plt.axvline(x_const, label=f"Line {i + 1}")

  # 绘制真实点
  plt.scatter(*true_point, color="green", s=120, marker="o", label="True Position")

  # 绘制估计点
  plt.scatter(*estimate, color="red", s=120, marker="x", label="Estimated Position")

  plt.xlabel("X")
  plt.ylabel("Y")
  plt.title("Dynamic Simulation of LLOP")
  plt.legend()
  plt.grid(True)
  plt.axis("equal")
  plt.show()


if __name__ == "__main__":
  np.random.seed(42)  # 保持可重复性

  true_point = (3.0, 4.0)  # 真实位置
  lines = generate_lines(true_point, n_lines=6, noise_level=0.3)
  estimate = linear_line_of_position(lines)

  print(f"True position: {true_point}")
  print(f"Estimated position: ({estimate[0]:.3f}, {estimate[1]:.3f})")

  plot_simulation(true_point, lines, estimate)