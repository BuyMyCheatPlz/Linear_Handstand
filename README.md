# 一阶线性倒立摆 LQR 控制器

C++ (Eigen) 实现 + MATLAB 数值验证 + 交互式动画仪表盘。

---

## 物理模型

小车-摆杆系统，小角度线性化 ($\sin\theta \approx \theta,\ \cos\theta \approx 1$)：

$$
\frac{d}{dt}
\begin{bmatrix} \theta \\ \dot{\theta} \end{bmatrix}
=
\begin{bmatrix}
0 & 1 \\[4pt]
\frac{(M+m)g}{Ml} & -d
\end{bmatrix}
\begin{bmatrix} \theta \\ \dot{\theta} \end{bmatrix}
+
\begin{bmatrix}
0 \\[4pt]
-\frac{1}{Ml}
\end{bmatrix}
u
$$

| 符号 | 含义 | 默认值 |
|------|------|--------|
| $M$ | 小车质量 | 1.0 kg |
| $m$ | 摆杆质量 | 0.1 kg |
| $l$ | 摆杆长度 | 0.5 m |
| $g$ | 重力加速度 | 9.81 m/s² |
| $d$ | 阻尼系数 | 0.05 |
| $\theta$ | 摆杆偏角 | — |
| $u$ | 控制力输入 | — |

开环极点: $\lambda \approx +4.62,\ -4.67$ (一个不稳定，需要 LQR 镇定)。

---

## LQR 设计

最小化二次代价函数：

$$
J = \int_0^\infty \left( x^T Q x + u^T R u \right) dt
$$

求解连续时间代数 Riccati 方程 (CARE)：

$$
A^T P + P A - P B R^{-1} B^T P + Q = 0
$$

最优反馈增益：$K = R^{-1} B^T P$，控制律：$u = -Kx$

| 参数 | 值 |
|------|-----|
| $Q$ | $\operatorname{diag}(100,\ 1)$ |
| $R$ | $0.01$ |
| $K$ | $[-111.37,\ -14.51]$ |
| 闭环极点 | $-11.34,\ -17.73$ |

---

## 文件结构

```
├── lqr_inverted_pendulum.h      # C++ 头文件: 模型结构体、LQR 求解器接口
├── lqr_inverted_pendulum.cpp    # C++ 实现: CARE (Hamiltonian+特征分解)、RK4 仿真
├── main.cpp                     # C++ 主程序: 求解+仿真+CSV 导出
├── CMakeLists.txt               # CMake 构建文件
│
├── verify_lqr.m                 # MATLAB 验证: lqr/care 对照、C++ 交叉比对、全图表
├── interactive_pendulum.m       # MATLAB 交互式仪表盘: 动画+扰动+相平面+时域图
│
├── pendulum_sim.csv             # C++ 仿真输出 (运行后生成)
├── lqr_gain_cpp.csv             # C++ LQR 增益 (运行后生成)
├── riccati_P_cpp.csv            # C++ Riccati 矩阵 (运行后生成)
│
└── README.md                    # 本文件
```

---

## 编译与运行 (C++)

### 前置依赖

- **Eigen 3** (header-only): https://eigen.tuxfamily.org/
- C++17 编译器 (GCC / Clang / MSVC)

### 直接编译

```bash
# Windows MSYS2 / Linux / macOS
g++ -std=c++17 -I/path/to/eigen -O2 -Wall \
    -o lqr_pendulum main.cpp lqr_inverted_pendulum.cpp
./lqr_pendulum
```

### CMake 编译

```bash
mkdir build && cd build
# 系统已安装 Eigen:
cmake ..
# 或手动指定 Eigen 路径:
cmake -DEIGEN3_INCLUDE_DIR="C:/path/to/eigen" ..
make
./lqr_pendulum
```

### 输出

程序运行后打印 LQR 求解结果、仿真指标，并生成 3 个 CSV 文件供 MATLAB 读取验证。

```
LQR 增益 K = [-111.37, -14.51]
闭环极点: -11.34, -17.73

初始角度: 15° → 调节时间: 0.43s → 稳态: 0.00°
最大控制力: 29.16 N
```

---

## MATLAB 验证

### `verify_lqr.m` — 数值验证与全图表

```matlab
verify_lqr
```

功能：

- 使用 `lqr()` 和 `care()` 独立求解，验证一致性
- 读取 C++ 输出的 CSV，逐点比对 $K$, $P$, 仿真轨迹
- 绘制 6 张图：角度 / 角速度 / 控制力 / 相平面 / C++ 差异 / 极点图
- 权重参数扫描 ($Q$, $R$ 对 $K$ 的影响分析)

### `interactive_pendulum.m` — 交互式仪表盘

```matlab
interactive_pendulum
```

功能：

- 倒立摆实时动画 (含小车位移)
- 相平面轨迹 ($\theta$ vs $\dot{\theta}$)
- 三条时域曲线: $\theta(t)$, $\dot{\theta}(t)$, $F(t)$
- 按钮 / 键盘施加扰动，观察 LQR 恢复过程

| 操作 | 按钮 | 快捷键 |
|------|------|--------|
| 推摆杆 (角速度脉冲) | 👆 推摆杆 | `P` |
| 推小车 (力脉冲, 持续 100ms) | 🚗 推小车 | `C` |
| 重置 | 🔄 重置 | `R` |
| 暂停/继续 | 暂停 | `空格` |
| 调扰动强度 | 滑块 | `↑` `↓` |

视觉反馈：

- 摆杆颜色: 竖直(白) → 倾斜(黄) → 大偏角(红)
- 控制力数字: 小力(绿) → 中力(橙) → 大力(红)
- 橙色 ▽ = 推摆杆时刻，绿色 △ = 推小车时刻
- ⚡ 闪效 = 扰动触发

---

## 算法细节

### CARE 求解 (C++)

采用 **Hamiltonian 矩阵 + 特征分解法**：

1. 构造 $2n \times 2n$ Hamiltonian 矩阵:
   $$
   H = \begin{bmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{bmatrix}
   $$

2. 特征分解 $H V = V \Lambda$，取 $n$ 个稳定特征向量 ($\operatorname{Re}(\lambda) < 0$)

3. 将稳定特征向量分为上下两块 $U_{11}, U_{21} \in \mathbb{C}^{n \times n}$

4. Riccati 解: $P = \operatorname{Re}(U_{21} U_{11}^{-1})$

5. 反馈增益: $K = R^{-1} B^T P$

### 仿真 (C++)

闭环系统 $\dot{x} = (A - BK)x$ 使用 **固定步长 RK4** 积分，步长 1ms。

### 性能优化 (MATLAB 交互式仪表盘)

- `animatedline` 增量绘图 (避免每帧全量重绘 ~6000 点)
- 主循环零 `guidata` 开销 (变量直接访问)
- 动画降频至 ~33fps，图表更新 ~16fps
- 环形缓冲区 + `MaximumNumPoints` 自动滚动

---

## 自定义参数

### 更改物理参数

编辑 `main.cpp` 中的 `PendulumParams` 或 `verify_lqr.m` 中的变量：

```cpp
params.M = 1.0;   // 小车质量
params.m = 0.1;   // 摆杆质量
params.l = 0.5;   // 摆杆长度
params.d = 0.05;  // 阻尼
```

### 调整 LQR 权重

```cpp
Q << 100.0, 0.0,    // q11 越大 → 角度恢复越快
       0.0, 1.0;    // q22 越大 → 抑制摆动速度
R << 0.01;          // R 越大 → 控制越保守 (力越小)
```

- $Q_{11} \uparrow$ → $K_1 \uparrow$ → 更快回到平衡，可能需要更大的控制力
- $R \uparrow$ → $K$ 幅值 $\downarrow$ → 调节更慢，控制力更小

---

## 许可

MIT License
