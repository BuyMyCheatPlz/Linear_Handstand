# 一阶线性倒立摆 LQR 控制器

C++ (Eigen) 实现 + MATLAB 数值验证 + 交互式动画仪表盘。

项目包含两个独立的 LQR 求解模块：
- **`LQR_Origin/`** — 离散时间 LQR 参考实现（多种 C++ 写法 + MATLAB/MEX 接口）
- **`LQR_test/`** — 一阶线性倒立摆连续时间 LQR 控制器（主项目）

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
├── .vscode/                          # VS Code 编辑器配置
│   ├── c_cpp_properties.json         #   C/C++ IntelliSense 配置（编译器路径、include 路径）
│   ├── launch.json                   #   调试启动配置（GDB 调试器设置）
│   └── settings.json                 #   C/C++ Runner 扩展配置（编译/链接选项、警告等级）
│
├── LQR_Origin/                       # ★ 离散时间 LQR 参考实现（独立模块）
│   ├── lqr.cpp                       #   纯 C++ LQR 求解器（自实现 Matrix 类，无外部依赖）
│   ├── lqr.exe                       #   lqr.cpp 的编译产物
│   ├── lqr_eigen.cpp                 #   Eigen 版 LQR 求解器（离散 Riccati 迭代，2 个算例）
│   ├── lqr_eigen.exe                 #   lqr_eigen.cpp 的编译产物
│   ├── lqr_matlab_interface.m        #   MATLAB 封装函数，通过 MEX 桥接调用 C++ 求解器
│   ├── lqr_mex.cpp                   #   MATLAB MEX 包装器（自实现 Matrix 类，纯 C++）
│   ├── lqr_mex.mexw64                #   lqr_mex.cpp 编译后的 MEX 二进制（Windows x64）
│   ├── lqr_mex_eigen.cpp             #   MATLAB MEX 包装器（Eigen 版，推荐用于 MATLAB 调用）
│   ├── run_lqr_sim.m                 #   MATLAB 驱动脚本：编译 MEX → 求解 LQR → 闭环仿真 → 绘图
│   └── test.mat                      #   测试数据文件
│
├── LQR_test/                         # ★ 一阶线性倒立摆 LQR 控制器（主项目）
│   ├── lqr_inverted_pendulum.h       #   C++ 头文件：物理模型结构体、LQR 求解器接口、仿真函数声明
│   ├── lqr_inverted_pendulum.cpp     #   C++ 实现：CARE（Hamiltonian + Schur 分解）、RK4 仿真、CSV 导出
│   ├── main.cpp                      #   C++ 主程序：参数配置 → LQR 求解 → 仿真 → 指标输出 → CSV 导出
│   ├── CMakeLists.txt                #   CMake 构建文件（自动查找 Eigen 或手动指定路径）
│   ├── lqr_pendulum.exe              #   编译产物（Release 版）
│   ├── lqr_pendulum_dbg.exe          #   编译产物（Debug 版，含调试符号）
│   │
│   ├── verify_lqr.m                  #   MATLAB 验证脚本：lqr/care 对照、C++ 交叉比对、6 张全图表
│   ├── interactive_pendulum.m        #   MATLAB 交互式仪表盘：动画 + 扰动 + 相平面 + 三条时域曲线
│   │
│   ├── pendulum_sim.csv              #   C++ 仿真轨迹输出（运行 main.cpp 后生成）
│   ├── lqr_gain_cpp.csv              #   C++ 求解的 LQR 反馈增益 K（运行后生成）
│   ├── riccati_P_cpp.csv             #   C++ 求解的 Riccati 矩阵 P（运行后生成）
│   └── test.csv                      #   测试用 CSV 数据
│
└── README.md                         # 本文件
```

---

## 模块详解

### `LQR_Origin/` — 离散时间 LQR 参考实现

该目录是独立的离散时间 LQR 求解器集合，展示了相同算法的多种实现方式：

| 文件 | 说明 |
|------|------|
| `lqr.cpp` | 完全自包含的 C++ 实现：手写 `Matrix` 类（含逆矩阵高斯-约旦消元）、**离散时间 Riccati 迭代**求解 DARE。零外部依赖，可直接 `g++ lqr.cpp -o lqr.exe` 编译运行。适合学习矩阵运算底层实现。 |
| `lqr_eigen.cpp` | 使用 Eigen 库的独立 C++ 实现，算法与 `lqr.cpp` 相同但代码更简洁。内置 2 个算例（2 阶积分器 + 3 阶双输入系统），验证闭环稳定性。推荐作为 Eigen 用法的参考。 |
| `lqr_mex.cpp` | MATLAB MEX 接口（自实现 Matrix 版）：将 `lqr.cpp` 的求解器包装为 `[K, P] = lqr_mex(A, B, Q, R)`，可在 MATLAB 中直接调用 C++ 加速。 |
| `lqr_mex_eigen.cpp` | MATLAB MEX 接口（Eigen 版）：功能同上，但使用 Eigen 库实现。编译后生成 `lqr_mex_eigen.mexw64`，性能优于纯 C++ 版。 |
| `lqr_matlab_interface.m` | MATLAB 侧的封装函数，对 `lqr_mex` 做参数校验和错误提示。 |
| `run_lqr_sim.m` | MATLAB 一键运行脚本：自动编译 MEX → 求解 LQR → 跑闭环仿真 → 绘制状态/控制曲线。 |

**关键区别**：`LQR_Origin/` 求解的是 **离散时间** 代数 Riccati 方程 (DARE: $P = Q + A^T P A - A^T P B (R + B^T P B)^{-1} B^T P A$)，使用**迭代法**；而 `LQR_test/` 求解的是 **连续时间** 代数 Riccati 方程 (CARE)，使用** Hamiltonian + Schur 分解法**。

### `LQR_test/` — 连续时间倒立摆 LQR 控制器（主项目）

#### C++ 核心

| 文件 | 说明 |
|------|------|
| `lqr_inverted_pendulum.h` | 头文件，定义 `PendulumParams`（物理参数）、`LqrGain`（LQR 结果）结构体，声明 `buildPendulumSystem()`（构造 A, B 矩阵）、`solveLQR()`（CARE 求解）、`simulate()`（RK4 仿真）、`exportCSV()`（数据导出）四个核心函数。 |
| `lqr_inverted_pendulum.cpp` | 实现文件：① 构造小车-摆杆线性化状态空间 A(2×2), B(2×1)；② **Hamiltonian 矩阵 + 特征分解法**求解 CARE → 得 P 和 K；③ 闭环系统 **固定步长 RK4** 数值积分（默认 1ms 步长，5 秒）；④ 导出 CSV 供 MATLAB 验证。 |
| `main.cpp` | 主程序入口：定义默认参数（M=1.0, m=0.1, l=0.5, Q=diag(100,1), R=0.01）→ 打印开环/闭环特征值 → 仿真初始角度 15° 的镇定过程 → 输出调节时间、最大控制力等指标 → 导出 3 个 CSV 文件。 |
| `CMakeLists.txt` | CMake 构建配置：C++17 标准，自动查找 Eigen3（支持 `find_package` 和手动指定 `EIGEN3_INCLUDE_DIR`），平台特定编译选项（MSVC 的 `_USE_MATH_DEFINES`，GCC 的 `-Wall -Wextra`）。 |

#### MATLAB 脚本

| 文件 | 说明 |
|------|------|
| `verify_lqr.m` | **数值验证 + 全图表脚本**：① 使用 `lqr()` 和 `care()` 分别独立求解，验证一致；② 读取 C++ 输出的 3 个 CSV，逐点比对 $K$, $P$, 仿真轨迹；③ 绘制 6 张子图（角度/角速度/控制力/相平面/C++ 差异/极点图）；④ 权重参数扫描分析（$Q_{11}$ 和 $R$ 对增益的影响）。 |
| `interactive_pendulum.m` | **交互式实时仿真仪表盘**：① 倒立摆动画（含小车位移、摆杆颜色随偏角变化）；② 相平面轨迹实时绘制；③ 三条时域曲线 $\theta(t)$, $\dot{\theta}(t)$, $F(t)$；④ 按钮 + 键盘快捷键施加扰动（推摆杆 / 推小车），观察 LQR 恢复过程；⑤ 性能优化：`animatedline` 增量更新 + 动画/图表降频 + 环形缓冲区 + 零 `guidata` 开销。 |

#### 生成的数据文件

| 文件 | 内容 | 生成方式 |
|------|------|----------|
| `pendulum_sim.csv` | 仿真轨迹：时间、角度 (rad)、角速度 (rad/s)、控制力 (N) | `main.cpp` 运行后自动生成 |
| `lqr_gain_cpp.csv` | C++ 求解的反馈增益 K (1×2) | `main.cpp` 运行后自动生成 |
| `riccati_P_cpp.csv` | C++ 求解的 Riccati 矩阵 P (2×2) | `main.cpp` 运行后自动生成 |

---

## 编译与运行 (C++)

### 前置依赖

- **Eigen 3** (header-only): https://eigen.tuxfamily.org/
- C++17 编译器 (GCC / Clang / MSVC)

### 主项目 (`LQR_test/`)

#### 直接编译

```bash
# Windows MSYS2 / Linux / macOS
cd LQR_test
g++ -std=c++17 -I/path/to/eigen -O2 -Wall \
    -o lqr_pendulum main.cpp lqr_inverted_pendulum.cpp
./lqr_pendulum
```

#### CMake 编译

```bash
cd LQR_test
mkdir build && cd build
# 系统已安装 Eigen:
cmake ..
# 或手动指定 Eigen 路径:
cmake -DEIGEN3_INCLUDE_DIR="C:/path/to/eigen" ..
make
./lqr_pendulum
```

### 参考模块 (`LQR_Origin/`)

```bash
# 纯 C++ 版（无外部依赖）
cd LQR_Origin
g++ -std=c++17 -O2 lqr.cpp -o lqr.exe && ./lqr.exe

# Eigen 版
g++ -std=c++17 -O2 -I/path/to/eigen lqr_eigen.cpp -o lqr_eigen.exe && ./lqr_eigen.exe

# MEX 版（在 MATLAB 中）
mex -R2018a lqr_mex.cpp
mex -R2018a lqr_mex_eigen.cpp -I/path/to/eigen/include
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

## MATLAB 验证与交互

### `verify_lqr.m` — 数值验证与全图表

```matlab
cd LQR_test
verify_lqr
```

功能：

- 使用 `lqr()` 和 `care()` 独立求解，验证一致性
- 读取 C++ 输出的 CSV，逐点比对 $K$, $P$, 仿真轨迹
- 绘制 6 张图：角度 / 角速度 / 控制力 / 相平面 / C++ 差异 / 极点图
- 权重参数扫描 ($Q$, $R$ 对 $K$ 的影响分析)

### `interactive_pendulum.m` — 交互式仪表盘

```matlab
cd LQR_test
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

### CARE 求解 (C++, 连续时间)

采用 **Hamiltonian 矩阵 + Schur 分解法**（`lqr_inverted_pendulum.cpp`）：

1. 构造 $2n \times 2n$ Hamiltonian 矩阵:
   $$
   H = \begin{bmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{bmatrix}
   $$

2. 实 Schur 分解 $H = U T U^T$，按实部负/正重排对角块

3. 取前 $n$ 列稳定 Schur 向量，分为上下两块 $U_{11}, U_{21} \in \mathbb{R}^{n \times n}$

4. Riccati 解: $P = U_{21} U_{11}^{-1}$

5. 反馈增益: $K = R^{-1} B^T P$

### DARE 求解 (C++, 离散时间)

采用 **Riccati 迭代法**（`LQR_Origin/` 中各文件）：

1. 初始化 $P_0 = Q$
2. 迭代 $P_{k+1} = Q + A^T P_k A - A^T P_k B (R + B^T P_k B)^{-1} B^T P_k A$
3. 收敛判据：$\|P_{k+1} - P_k\|_\infty < 10^{-10}$
4. 增益: $K = (R + B^T P B)^{-1} B^T P A$

### 仿真 (C++)

闭环系统 $\dot{x} = (A - BK)x$ 使用 **固定步长 RK4** 积分，步长 1ms。

### 性能优化 (MATLAB 交互式仪表盘)

- `animatedline` 增量绘图 (避免每帧全量重绘 ~6000 点)
- 主循环零 `guidata` 开销 (变量直接访问)
- 动画降频至 ~33fps，图表更新 ~16fps
- 环形缓冲区 + `MaximumNumPoints` 自动滚动

---

## LQR_Origin 与 LQR_test 的关系

| 维度 | `LQR_Origin/` | `LQR_test/` |
|------|---------------|-------------|
| 时间域 | **离散时间** | **连续时间** |
| 状态方程 | $x_{k+1} = A x_k + B u_k$ | $\dot{x} = A x + B u$ |
| Riccati 方程 | DARE（迭代求解） | CARE（Schur 分解） |
| 求解方法 | Riccati 迭代法 | Hamiltonian + Schur 分解 |
| 应用场景 | 通用离散 LQR 参考 / MEX 加速 | 倒立摆物理系统镇定控制 |
| MATLAB 接口 | MEX 包装器 | 独立 `.m` 验证 + 可视化 |

`LQR_Origin/` 是通用 LQR 求解器的**参考实现和工具箱**（纯 C++、Eigen、MEX 三种形态）；`LQR_test/` 是针对一阶倒立摆这个**具体物理系统**的控制设计、仿真和可视化。

---

## 自定义参数

### 更改物理参数

编辑 `LQR_test/main.cpp` 中的 `PendulumParams` 或 `verify_lqr.m` 中的变量：

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
