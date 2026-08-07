# room_ac_turb.m 公式集（LaTeX 版）

> 本文件整理 `room_ac_turb.m` 中出现的全部公式，逐条编号、标注代码行号，并对关键公式补充物理推导。所有公式采用 LaTeX 行间格式 `$$...$$`，可在支持 KaTeX/MathJax 的 Markdown 渲染器中显示。

---

## 0. 符号表

### 0.1 房间与空气参数

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $L$ | 房间边长 | `model.L` | m |
| $N$ | 每方向网格数 | `model.N` | — |
| $\Delta x$ | 网格尺寸 | `dx` | m |
| $\Delta t$ | 时间步长 | `dt` | s |
| $\rho$ | 空气密度 | `model.rho` | kg/m³ |
| $c_p$ | 空气比热 | `model.cp` | J/(kg·K) |
| $T_0$ | 初始温度 | `model.T0` | K |
| $T_{set}$ | 目标温度 | `model.T_set` | K |
| $C_{cell}$ | 单格热容 | `cellCapacity` | J/K |

### 0.2 射流参数

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $\alpha_{bg}$ | 背景扩散系数 | `model.alpha_bg` | m²/s |
| $\alpha_{jet}$ | 射流附加扩散 | `model.alpha_jet` | m²/s |
| $v_0$ | 射流名义速度 | `model.jet_velocity_nominal` | m/s |
| $\sigma_0$ | 射流初始半宽 | `model.jet_initial_half_width` | m |
| $k_s$ | 射流扩展率 | `model.jet_spread_rate` | — |
| $\lambda$ | 射流衰减长度比 | `model.jet_decay_length_ratio` | — |
| $y_{in}$, $z_{in}$ | 进风口位置 | `model.inlet_y`, `model.inlet_z` | m |
| $\theta$ | 射流角度 | `caseDef.angle_deg` | deg |

### 0.3 空调参数

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $Q_{nom}$ | 额定制冷量 | `model.cooling_capacity_nominal` | W |
| $\dot m_{nom}$ | 额定质量流量 | `model.m_dot_nominal` | kg/s |
| $T_{sup,min}$ | 最低送风温度 | `model.T_supply_min` | K |
| $\text{COP}_0$ | 满载 COP | `model.COP_full_load` | — |
| $\Delta\text{COP}$ | 部分负载 COP 增益 | `model.COP_part_load_gain` | — |
| $P_{fan,min}$ | 风机最小功率 | `model.fan_power_min` | W |
| $P_{fan,nom}$ | 风机额定功率 | `model.fan_power_nominal` | W |
| $f_{r,min}$ | 风机最小比例 | `model.fan_ratio_min` | — |
| $r_{min}$ | 压缩机最小比例 | `model.compressor_ratio_min` | — |

### 0.4 控制器参数

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $K_p$ | 比例系数 | `model.controller_Kp` | 1/K |
| $K_i$ | 积分系数 | `model.controller_Ki` | 1/(K·s) |
| $T_{band}$ | 积分作用带 | `model.controller_integral_band` | K |

### 0.5 热负荷参数

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $UA$ | 围护总传热 | `model.envelope_UA` | W/K |
| $\bar T_{out}$ | 室外平均温度 | `model.T_outdoor_mean` | K |
| $A_{out}$ | 室外温度日幅值 | `model.T_outdoor_amplitude` | K |
| $\phi_{out}$ | 室外温度相位 | `model.T_outdoor_phase` | rad |
| $N_p$ | 人数 | `model.people_count` | — |
| $q_p$ | 每人显热 | `model.people_sensible_each` | W |
| $q_{eq}$ | 设备热负荷 | `model.equipment_load` | W |
| $q_{lit}$ | 照明热负荷 | `model.lighting_load` | W |
| $q_{sol}$ | 太阳辐射峰值 | `model.solar_peak` | W |

### 0.6 状态量

| 符号 | 含义 | 代码变量 | 单位 |
|------|------|---------|------|
| $T_{i,j,k}^n$ | 时刻 $n\Delta t$ 网格温度 | `T` | K |
| $r$ | 压缩机比例 | `capacityRatio` | — |
| $f_r$ | 风机比例 | `fanRatio` | — |
| $\dot Q_c$ | 制冷功率 | `coolingPower` | W |
| $P_{el}$ | 电功率 | `electricPower` | W |
| $T_{sup}$ | 送风温度 | `supplyTemperature` | K |
| $T_{ret}$ | 回风温度 | `returnTemperature` | K |
| $\dot Q_{load}$ | 总热负荷 | `heatLoad` | W |

---

## 1. 网格与时间离散

### 1.1 网格尺寸

$$
\Delta x = \frac{L}{N} \qquad \text{(:70)}
$$

### 1.2 网格中心坐标

$$
x_i = \left(i - 0.5\right)\Delta x, \quad i=1,\dots,N \qquad \text{(:87)}
$$

y、z 方向同理（`:88-89`）。

### 1.3 单格热容

$$
C_{cell} = \rho\, c_p\, \Delta x^3 \qquad \text{(:90)}
$$

> **物理意义**：让一个网格的空气升温 1K 所需的能量。数值上 $C_{cell} = 1.2 \times 1005 \times 0.5^3 = 150.75\ \text{J/K}$。

### 1.4 总步数

$$
N_t = \text{round}\left(\frac{t_{end}}{\Delta t}\right) \qquad \text{(:79)}
$$

### 1.5 采样间隔（步数）

$$
n_{sample} = \max\left(1,\ \text{round}\left(\frac{t_{sample}}{\Delta t}\right)\right) \qquad \text{(:81)}
$$

### 1.6 浮力搅拌间隔（步数）

$$
n_{mix} = \max\left(1,\ \text{round}\left(\frac{1}{\Delta t}\right)\right) \qquad \text{(:84)}
$$

### 1.7 进度打印间隔

$$
n_{prog} = \max\left(1,\ \text{round}\left(\frac{N_t}{10}\right)\right) \qquad \text{(:85)}
$$

---

## 2. 稳定性判据

### 2.1 扩散稳定性上限

$$
\Delta t \le \frac{\Delta x^2}{6\,\alpha_{max}}, \quad \alpha_{max} = \alpha_{bg} + \alpha_{jet} \qquad \text{(:72-77)}
$$

> **推导**：显式 Euler 离散扩散方程 $\partial T/\partial t = \alpha \nabla^2 T$，三维情况下 von Neumann 稳定性分析给出放大因子 $G = 1 - 6\alpha\Delta t/\Delta x^2(\cos\xi-1)$，要求 $|G|\le 1$，最严条件 $\Delta t \le \Delta x^2/(6\alpha)$。代码里若违反直接报错（`:74-77`）。
>
> **数值**：$\Delta x^2/(6\alpha_{max}) = 0.25/(6\times 0.043) = 0.969$ s。默认 $\Delta t=0.2$ s 满足。

### 2.2 对流 CFL 数

$$
\text{CFL} = \frac{\Delta t}{\Delta x}\, \max_{i,j,k}\left(|u_{i,j,k}| + |w_{i,j,k}|\right) \qquad \text{(:126)}
$$

### 2.3 扩散数

$$
D = \frac{6\,\alpha_{max}\,\Delta t}{\Delta x^2} \qquad \text{(:127)}
$$

### 2.4 组合稳定性数

$$
S = \text{CFL} + D \le 0.95 \qquad \text{(:128-129)}
$$

> **推导**：对流+扩散方程 $\partial T/\partial t + u\partial T/\partial x = \alpha\nabla^2 T$，显式格式的稳定性条件近似为 CFL 与扩散数的加权和不超过 1，代码取保守阈值 0.95。若超出，建议步长为：
>
> $$\Delta t_{sugg} = \frac{0.95}{\max(|u|+|w|)/\Delta x + 6\alpha_{max}/\Delta x^2} \qquad \text{(:130)}$$

---

## 3. 热负荷空间分布

> 函数：`buildHeatLoadFields`（`:327-341`）。把全屋的传热、内部热源、太阳辐射分配到每个网格。

### 3.1 外墙贴面数

$$
n_{face,i,j,k} = \mathbb{1}[i=1] + \mathbb{1}[i=N] + \mathbb{1}[j=1] + \mathbb{1}[j=N] + \mathbb{1}[k=1] + \mathbb{1}[k=N] \qquad \text{(:329-332)}
$$

> 角格 $n_{face}=3$，棱格 $n_{face}=2$，面心格 $n_{face}=1$，内部格 $n_{face}=0$。

### 3.2 围护传热分布

$$
G_{i,j,k} = UA \cdot \frac{n_{face,i,j,k}}{\displaystyle\sum_{i,j,k} n_{face,i,j,k}} \qquad \text{(:333)}
$$

> $G$ 是每个网格的"局部传热系数"[W/K]，保证 $\sum G = UA$。

### 3.3 人员活动区权重

$$
w^{int}_{i,j,k} = \frac{\mathbb{1}[z_k \le 2.0\,\text{m}]}{\displaystyle\sum \mathbb{1}[z_k \le 2.0\,\text{m}]} \qquad \text{(:336-337)}
$$

> 内部热源（人/设备/灯）只加在离地 2m 以内的网格，权重归一。

### 3.4 太阳辐射权重

$$
w^{sol}_{i,j,k} = \frac{\mathbb{1}[x_i = x_N] + 0.5\,\mathbb{1}[z_k = z_N]}{\displaystyle\sum \left(\mathbb{1}[x_i = x_N] + 0.5\,\mathbb{1}[z_k = z_N]\right)} \qquad \text{(:339-340)}
$$

> 太阳加在东墙（$x=L$）权重 1，顶棚（$z=L$）权重 0.5。

---

## 4. 射流几何与形状

> 函数：`buildJetFields`（`:343-374`）。建立空调射流的几何中心线、横截面分布、轴向衰减。

### 4.1 射流中心线（无约束）

$$
z_c^{*}(x) = z_{in} + x \tan\theta \qquad \text{(:350)}
$$

### 4.2 撞墙判据

$$
\text{wall}_{jet}(x) = \mathbb{1}\left[z_c^{*}(x) \le 0\right] \vee \mathbb{1}\left[z_c^{*}(x) \ge L\right] \qquad \text{(:351)}
$$

### 4.3 射流中心线（实际，含贴壁）

$$
z_c(x) = \text{clip}\left(z_c^{*}(x),\ 0,\ L\right) \qquad \text{(:352)}
$$

> 撞顶/撞底后射流"夹住"在边界，转为水平贴壁流动（通过后续 `verticalSlope=0` 实现，`:371`）。

### 4.4 射流半宽

$$
\sigma(x) = \sigma_0 + k_s\, x \qquad \text{(:354)}
$$

> 射流越往前越宽（线性扩展），$k_s=0.10$ 表示每米半宽增加 0.1m。

### 4.5 径向距离平方

$$
r^2(x,y,z) = (y - y_{in})^2 + (z - z_c(x))^2 \qquad \text{(:357)}
$$

### 4.6 射流核心高斯分布

$$
W_{core}(x,y,z) = \exp\left(-\frac{r^2}{2\,\sigma(x)^2}\right) \qquad \text{(:358)}
$$

> **推导**：自由湍射流的横截面速度/浓度分布实验上呈高斯型。从动量守恒 + 涡黏性假设可推得 $\bar u/U_c \sim \exp(-r^2/(2\sigma^2))$。这里把速度分布直接用作"冷量分布权重"。

### 4.7 轴向衰减

$$
W_{decay}(x) = \exp\left(-\frac{x}{L\,\lambda}\right) \qquad \text{(:359)}
$$

> 射流强度沿 x 指数衰减，特征长度 $L\lambda = 10\times 0.8 = 8$ m。

### 4.8 综合射流权重

$$
W(x,y,z) = W_{core}(x,y,z) \cdot W_{decay}(x) \qquad \text{(:360)}
$$

### 4.9 射流附加扩散场

$$
\alpha_{jet}^{field}(x,y,z) = \alpha_{jet} \cdot W(x,y,z) \qquad \text{(:361)}
$$

### 4.10 冷量分配权重

$$
\tilde w_{cool}(x,y,z) = 0.002 + W(x,y,z) \qquad \text{(:364)}
$$

$$
w_{cool}(x,y,z) = \frac{\tilde w_{cool}(x,y,z)}{\displaystyle\sum_{i,j,k} \tilde w_{cool}(x_i,y_j,z_k)} \qquad \text{(:365)}
$$

> 0.002 是给远离射流的角落留的"呼吸口"背景冷量。归一化保证 $\sum w_{cool} = 1$，使总冷量严格等于空调制冷功率。

---

## 5. 射流速度场

### 5.1 x 方向射流速度

$$
u_{jet}(x,y,z) = v_0\, W(x,y,z) \qquad \text{(:368)}
$$

### 5.2 x 方向净速度（去均值，保证体积守恒）

$$
u(x,y,z) = u_{jet}(x,y,z) - \langle u_{jet}\rangle_{y,z} \qquad \text{(:369)}
$$

> 其中 $\langle\cdot\rangle_{y,z}$ 表示对 y、z 方向取平均（每个 x 切面）。这保证每个 x 截面 $\sum_{j,k} u = 0$，即射流核心向前推 + 外围等量回流，封闭房间总体积流量为零。

### 5.3 z 方向斜率

$$
s_z(x) = \tan\theta \cdot \mathbb{1}[\neg\,\text{wall}_{jet}(x)] \qquad \text{(:370-371)}
$$

> 未撞墙时 $s_z = \tan\theta$；撞墙后 $s_z=0$（水平贴壁）。

### 5.4 z 方向射流速度

$$
w_{jet}(x,y,z) = u_{jet}(x,y,z) \cdot s_z(x) \qquad \text{(:372)}
$$

### 5.5 z 方向净速度

$$
w(x,y,z) = w_{jet}(x,y,z) - \langle w_{jet}\rangle_{y,z} \qquad \text{(:373)}
$$

---

## 6. 热负荷计算

> 函数：`heatLoadRate`（`:376-388`）。每步计算室外温度、太阳辐射、各网格热源。

### 6.1 日相位

$$
\phi_{daily}(t) = \frac{2\pi\, t}{86400} \qquad \text{(:378)}
$$

### 6.2 室外温度

$$
T_{out}(t) = \bar T_{out} + A_{out}\, \sin\left(\phi_{daily}(t) + \phi_{out}\right) \qquad \text{(:379-380)}
$$

### 6.3 太阳辐射

$$
\dot q_{sol}(t) = q_{sol}\, \max\left(0,\ 0.70 + 0.30\, \sin\left(\phi_{daily}(t) + \frac{\pi}{3}\right)\right) \qquad \text{(:381)}
$$

> 太阳辐射在 350W~500W 之间波动，相位与室外温度错开 $\pi/3 - \phi_{out}$。

### 6.4 内部热源

$$
\dot q_{int} = N_p\, q_p + q_{eq} + q_{lit} \qquad \text{(:382-383)}
$$

> 默认 $4\times 75 + 400 + 200 = 900$ W。

### 6.5 单格热源率

$$
\dot q_{i,j,k}(t) = G_{i,j,k}\,\left(T_{out}(t) - T_{i,j,k}\right) + \dot q_{int}\, w^{int}_{i,j,k} + \dot q_{sol}(t)\, w^{sol}_{i,j,k} \qquad \text{(:385-386)}
$$

### 6.6 全屋总热负荷

$$
\dot Q_{load}(t) = \sum_{i,j,k} \dot q_{i,j,k}(t) \qquad \text{(:387)}
$$

---

## 7. 压缩机控制器

> 函数：`compressorController`（`:390-420`）。仅变频机使用。

### 7.1 定频机输出

$$
r = 1 \qquad \text{(:393)}
$$

> 定频机永远满载。

### 7.2 前馈项

$$
r_{ff} = \text{clip}\left(\frac{\dot Q_{load}}{Q_{nom}},\ 0,\ 1\right) \qquad \text{(:397)}
$$

### 7.3 误差信号

$$
e(t) = T_{ret}(t) - T_{set} \qquad \text{(:196 传入)}
$$

### 7.4 积分更新（条件积分）

若 $|e| \le T_{band}$：

$$
I^{*} = \text{clip}\left(I + K_i\, e\, \Delta t,\ -0.5,\ 0.5\right) \qquad \text{(:399-400)}
$$

否则：

$$
I^{*} = I \qquad \text{(:402)}
$$

### 7.5 抗饱和条件

$$
\text{accept} = \neg\Big[\big(r_{raw} > 1 \wedge e > 0\big) \vee \big(r_{raw} < 0 \wedge e < 0\big)\Big] \qquad \text{(:407-408)}
$$

其中：

$$
r_{raw} = r_{ff} + K_p\, e + I^{*} \qquad \text{(:404)}
$$

若 $\text{accept}=\text{true}$，则 $I \leftarrow I^{*}$（`:409`）。

### 7.6 输出比例

$$
r = \text{clip}\left(r_{ff} + K_p\, e + I,\ 0,\ 1\right) \qquad \text{(:411)}
$$

### 7.7 最小运行比例保护

若 $0 < r < r_{min}$：

$$
r = \begin{cases}
r_{min}, & \text{if } e > -0.25\ \text{K}\ \vee\ r_{ff} > r_{min}/2 \\
0, & \text{otherwise}
\end{cases} \qquad \text{(:413-419)}
$$

---

## 8. 空调性能

> 函数：`airConditionerPerformance`（`:422-443`）。

### 8.1 风机比例

$$
f_r = f_{r,min} + (1 - f_{r,min})\, \sqrt{r} \qquad \text{(:424)}
$$

> 用 $\sqrt{r}$ 而非线性，因风机风量与转速近似一次关系，而转速与功率近似三次关系，取折中。

### 8.2 质量流量

$$
\dot m = \dot m_{nom}\, f_r \qquad \text{(:425)}
$$

### 8.3 风机功率

$$
P_{fan} = P_{fan,min} + (P_{fan,nom} - P_{fan,min})\, f_r^3 \qquad \text{(:426)}
$$

> 风机相似律：$P \propto n^3$。

### 8.4 压缩机停机情形

若 $r \le 0$：

$$
\dot Q_c = 0,\quad \text{COP} = \text{NaN},\quad T_{sup} = T_{ret},\quad P_{el} = P_{fan} \qquad \text{(:428-434)}
$$

### 8.5 请求制冷量

$$
\dot Q_{req} = r\, Q_{nom} \qquad \text{(:436)}
$$

### 8.6 风侧制冷极限

$$
\dot Q_{air} = \dot m\, c_p\, \max\left(0,\ T_{ret} - T_{sup,min}\right) \qquad \text{(:437)}
$$

> **物理意义**：风量 × 比热 × 最大温差 = 这股空气最多能被冷却多少。即使压缩机请求更大制冷，风侧也搬不走更多冷量。

### 8.7 实际制冷功率

$$
\dot Q_c = \min\left(\dot Q_{req},\ \dot Q_{air}\right) \qquad \text{(:438)}
$$

### 8.8 送风温度

$$
T_{sup} = T_{ret} - \frac{\dot Q_c}{\dot m\, c_p} \qquad \text{(:439)}
$$

### 8.9 部分负载 COP

$$
\text{COP} = \text{COP}_0 + \Delta\text{COP}\,(1 - r) \qquad \text{(:440-441)}
$$

> 低负载时 COP 更高（变频机省电的关键）。满载 $r=1$ 时 COP=3.2，停机附近 $r\to 0$ 时 COP=3.7。

### 8.10 电功率

$$
P_{el} = \frac{\dot Q_c}{\text{COP}} + P_{fan} \qquad \text{(:442)}
$$

---

## 9. 对流增量（迎风格式）

> 函数：`advectionIncrement`（`:445-460`）。

### 9.1 对流方程

$$
\frac{\partial T}{\partial t} = -u\frac{\partial T}{\partial x} - w\frac{\partial T}{\partial z} \qquad \text{(:446 注释)}
$$

### 9.2 速度分解

$$
u^+ = \max(u, 0),\quad u^- = \min(u, 0) \qquad \text{(:136-137)}
$$

$$
w^+ = \max(w, 0),\quad w^- = \min(w, 0) \qquad \text{(:138-139)}
$$

### 9.3 x 方向迎风离散

$$
\left.\Delta T_{adv}\right|_x = -\frac{\Delta t}{\Delta x}\left[ u^+_{i}\,(T_i - T_{i-1}) + u^-_{i}\,(T_{i+1} - T_i) \right] \qquad \text{(:449-451)}
$$

> **推导**：对 $\partial T/\partial x$ 用一阶迎风差分。$u>0$ 时信息从左向右传，用上游 $T_{i-1}$；$u<0$ 时反向，用 $T_{i+1}$。从 Taylor 展开 $T_{i\pm1} = T_i \pm \Delta x\, T' + O(\Delta x^2)$，迎风方向的一阶差分截断误差为 $O(\Delta x)$。

### 9.4 z 方向迎风离散

$$
\left.\Delta T_{adv}\right|_z = -\frac{\Delta t}{\Delta x}\left[ w^+_{k}\,(T_k - T_{k-1}) + w^-_{k}\,(T_{k+1} - T_k) \right] \qquad \text{(:453-455)}
$$

### 9.5 投影去均值

$$
\Delta T_{adv} \leftarrow \Delta T_{adv} - \overline{\Delta T_{adv}} \qquad \text{(:459)}
$$

> **原因**：简化速度场不逐格严格无散度（$\nabla\cdot\vec u \ne 0$），对流项可能产生虚假净加热/冷却。投影掉全域均值，保证对流只搬运不创造能量。

### 9.6 速度缩放

实际计算时速度乘以风机比例：

$$
\text{scale} = \frac{\Delta t}{\Delta x}\, f_r \qquad \text{(:205)}
$$

---

## 10. 扩散增量（守恒离散）

> 函数：`diffusionIncrement`（`:462-477`）。

### 10.1 扩散方程

$$
\frac{\partial T}{\partial t} = \nabla\cdot\left(\alpha\, \nabla T\right) \qquad \text{(:463 注释)}
$$

### 10.2 面心扩散系数（射流部分）

$$
\alpha^{face}_{x,i+\frac12} = \frac{1}{2}\left(\alpha_{jet,i}^{field} + \alpha_{jet,i+1}^{field}\right) \qquad \text{(:144)}
$$

y、z 方向同理（`:145-146`）。

### 10.3 背景扩散系数

$$
\alpha_{bg}^{coeff} = \frac{\Delta t}{\Delta x^2}\, \alpha_{bg} \qquad \text{(:142)}
$$

### 10.4 x 方向通量

$$
F_{x,i+\frac12} = \left(\alpha_{bg}^{coeff} + f_r\, \frac{\Delta t}{\Delta x^2}\, \alpha^{face}_{x,i+\frac12}\right)\left(T_{i+1} - T_i\right) \qquad \text{(:466)}
$$

### 10.5 x 方向增量

$$
\left.\Delta T_{diff}\right|_x\Big|_i = F_{x,i-\frac12} - F_{x,i+\frac12} \qquad \text{(:467-468)}
$$

> **推导**：从通量守恒 $\partial T/\partial t = -\nabla\cdot\vec F$，控制体积分得 $\Delta T_i = (F_{in} - F_{out})\Delta t / C_{cell}$。每个面的通量从一侧加、另一侧减，严格保证 $\sum \Delta T_{diff} = 0$（能量守恒）。

### 10.6 y、z 方向

同理（`:470-476`）。

---

## 11. 浮力搅拌

> 函数：`applyBuoyancyMixing`（`:479-497`）。每 1 秒做一次。

### 11.1 垂直温差

$$
\Delta T_k = T_{k+1} - T_k \qquad \text{(:481)}
$$

### 11.2 不稳定判据（上冷下热）

$$
\text{unstable}_k = \mathbb{1}[\Delta T_k < -0.3\,\text{K}] \qquad \text{(:484)}
$$

### 11.3 不稳定混合系数

$$
\text{Ri}_k = \frac{|\Delta T_k|}{2}, \quad \mu_k = \min\left(0.15,\ 0.05\,(1 + \text{Ri}_k)\right) \qquad \text{(:485-486)}
$$

### 11.4 稳定判据（上热下冷）

$$
\text{stable}_k = \mathbb{1}[\Delta T_k > 0.3\,\text{K}] \qquad \text{(:488)}
$$

### 11.5 稳定混合系数

$$
\text{Ri}_k = \frac{\Delta T_k}{2}, \quad \mu_k = \max\left(0.005,\ \frac{0.03}{1 + \text{Ri}_k}\right) \qquad \text{(:489-490)}
$$

> **推导**：Richardson 数 $Ri = \frac{g\,\beta\,\partial T/\partial z}{(\partial u/\partial z)^2}$ 衡量浮力与剪切之比。$Ri>0$ 稳定（分层），$Ri<0$ 不稳定（翻腾）。本模型简化为只用温差近似 $Ri \approx |\Delta T|/2$，并通过 $\mu$ 控制层间交换强度。

### 11.6 界面交换量

$$
E_k = \frac{1}{2}\, \mu_k\, \Delta T_k \qquad \text{(:492)}
$$

### 11.7 温度更新

$$
T_k \leftarrow T_k + E_k, \quad T_{k+1} \leftarrow T_{k+1} - E_k \qquad \text{(:494-495)}
$$

> 一加一减，能量严格守恒。

---

## 12. 能量更新主步

> 代码：`:200-208`。

### 12.1 累计耗电

$$
E_{el}(t) = \int_0^t P_{el}(\tau)\, d\tau \approx E_{el}(t-\Delta t) + P_{el}\, \Delta t \qquad \text{(:200)}
$$

### 12.2 源项温度增量

$$
\Delta T_{src,i,j,k} = \frac{\Delta t}{C_{cell}}\left(\dot q_{i,j,k} - \dot Q_c\, w_{cool,i,j,k}\right) \qquad \text{(:203)}
$$

> **推导**：从能量守恒 $C_{cell}\, dT/dt = \dot q - \dot Q_c\, w_{cool}$，显式 Euler 离分得 $\Delta T = \Delta t \cdot (\dot q - \dot Q_c w_{cool})/C_{cell}$。注意 $w_{cool}$ 归一化保证 $\sum w_{cool}=1$，所以总冷量严格等于 $\dot Q_c$。

### 12.3 温度更新

$$
T^{n+1}_{i,j,k} = T^n_{i,j,k} + \Delta T_{src} + \Delta T_{adv} + \Delta T_{diff} \qquad \text{(:208)}
$$

### 12.4 浮力修正（每 $n_{mix}$ 步）

$$
T^{n+1} \leftarrow \text{applyBuoyancyMixing}\left(T^{n+1}\right) \qquad \text{(:210-212)}
$$

### 12.5 数值发散检查

$$
\text{if }\exists\, (i,j,k):\ T_{i,j,k} \notin \mathbb{R}\ \text{finite},\ \text{then error} \qquad \text{(:214-217)}
$$

---

## 13. 统计量

### 13.1 平均温度

$$
\bar T = \frac{1}{N^3}\sum_{i,j,k} T_{i,j,k} \qquad \text{(:175)}
$$

### 13.2 温度标准差

$$
\sigma_T = \sqrt{\frac{1}{N^3}\sum_{i,j,k}\left(T_{i,j,k} - \bar T\right)^2} \qquad \text{(:176)}
$$

### 13.3 温度极差

$$
\Delta T_{range} = \max T - \min T \qquad \text{(:177)}
$$

### 13.4 最大速度

$$
|\vec v|_{max} = \max_{i,j,k}\sqrt{u^2_{i,j,k} + w^2_{i,j,k}} \qquad \text{(:266)}
$$

### 13.5 扩散系数范围

$$
\alpha_{range} = \left[\alpha_{bg},\ \alpha_{bg} + \max \alpha_{jet}^{field}\right] \qquad \text{(:269-270)}
$$

---

## 14. 动画帧采样

> 函数：`createAnimation`（`:558-613`）。

### 14.1 起始帧

$$
t_{first} = \min\left(5\,\text{min},\ t_{end}\right) \qquad \text{(:564)}
$$

### 14.2 帧步长

$$
s = \max\left(1,\ \left\lceil\frac{N_{cand}}{N_{max}}\right\rceil\right) \qquad \text{(:567)}
$$

其中 $N_{cand}$ 为候选帧数，$N_{max}$ 为最大帧数（默认 120）。

### 14.3 选用帧

$$
\mathcal{F} = \left\{ f_1, f_{1+s}, f_{1+2s}, \dots, f_{end} \right\} \qquad \text{(:568-571)}
$$

---

## 附录：公式依赖流程图

```
[输入参数]
    ↓
[1. 网格离散] → Δx, C_cell
    ↓
[2. 稳定性检查] → Δt 合法性
    ↓
[3. 热负荷分布] → G, w_int, w_sol   (一次性预计算)
    ↓
[4. 射流场] → W, w_cool, u, w       (每工况一次)
    ↓
[2. CFL 检查] → 速度场稳定性
    ↓
═══════════ 主循环 (每 Δt) ═══════════
    ↓
[6. 热负荷] → q_i, Q_load, T_out
    ↓
[7. 控制器] → r                     (变频机)
    ↓
[8. 空调性能] → Q_c, T_sup, COP, P_el, f_r
    ↓
[12.1 累计耗电] → E_el
    ↓
[12.2 源项] → ΔT_src
    ↓
[9. 对流] → ΔT_adv                 (用 f_r 缩放速度)
    ↓
[10. 扩散] → ΔT_diff                (用 f_r 缩放射流扩散)
    ↓
[12.3 温度更新] → T^{n+1}
    ↓
[11. 浮力搅拌] → T^{n+1} 修正        (每 n_mix 步)
    ↓
[12.5 发散检查]
    ↓
[13. 统计] → 记录采样
    ↓
═══════════ 循环结束 ═══════════
    ↓
[14. 动画帧采样]
    ↓
[输出 PNG/GIF/MAT]
```

---

## 附录：关键推导补充说明

### A.1 主步能量更新的守恒性

**总能量守恒证明**：

对全屋求和 $\sum_{i,j,k} \Delta T_{src}$：

$$
\sum \Delta T_{src} = \frac{\Delta t}{C_{cell}}\left(\sum \dot q_{i,j,k} - \dot Q_c \sum w_{cool,i,j,k}\right) = \frac{\Delta t}{C_{cell}}\left(\dot Q_{load} - \dot Q_c\right)
$$

因 $\sum w_{cool} = 1$（`buildJetFields` 归一化）。这说明：**全屋平均温度的变化仅由"热负荷 - 制冷量"决定**，与冷量空间分布无关。冷量分布只影响温度的**空间均匀度**，不影响全屋平均——这是能量守恒的直接结果。

对流项 $\sum \Delta T_{adv} = 0$（已投影去均值，`:459`）；扩散项 $\sum \Delta T_{diff} = 0$（通量守恒离散，每面一加一减）；浮力搅拌 $\sum \Delta T_{buoy} = 0$（一加一减）。三者都不改变全屋平均温度，只重新分配空间分布。

### A.2 迎风格式的截断误差

对 $\partial T/\partial x$ 在 $u>0$ 时用 $T_i - T_{i-1}$ 近似：

Taylor 展开：
$$
T_{i-1} = T_i - \Delta x\, T'_i + \frac{\Delta x^2}{2} T''_i - O(\Delta x^3)
$$

$$
\frac{T_i - T_{i-1}}{\Delta x} = T'_i - \frac{\Delta x}{2} T''_i + O(\Delta x^2)
$$

截断误差 $-\frac{\Delta x}{2} T''_i$ 等效于"数值扩散" $\alpha_{num} = \frac{u\Delta x}{2}$，使一阶迎风具有**数值耗散**特性——会人为增加混合。本模型 $\alpha_{num} = 1.2 \times 0.5 / 2 = 0.3$ m²/s，比 $\alpha_{jet}=0.04$ 大一个量级，但对趋势对比影响有限。

### A.3 CFL 条件的物理意义

CFL 数 $\text{CFL} = u\Delta t/\Delta x$ 表示一个时间步内信息（被气流携带的温度扰动）传播的网格数。物理上要求 $\text{CFL}\le 1$：信息不能在一个步长内"跨过"一个网格，否则显式格式无法"看到"上游变化，导致数值震荡发散。

### A.4 高斯分布的物理依据

自由湍射流的时均速度分布实验上满足：

$$
\frac{\bar u(r)}{U_c} = \exp\left(-\frac{r^2}{2\sigma^2}\right)
$$

其中 $U_c$ 是中心速度，$\sigma$ 是射流半宽。这可从涡黏性假设 $\nu_t = \text{const}$ + 边界层方程 $\bar u \partial \bar u/\partial x = \nu_t \frac{1}{r}\partial/\partial r(r\,\partial \bar u/\partial r)$ 的相似解推出（Tollmien 1926, Görtler 1942）。本模型直接用此分布作为冷量权重，假设冷量随射流核心一起输运。

### A.5 Richardson 数与稳定性

真实 Richardson 数：

$$
Ri = \frac{g\,\beta\,\partial T/\partial z}{(\partial u/\partial z)^2}
$$

- $Ri > 0$：稳定分层（热气在上），湍流被抑制
- $Ri < 0$：不稳定（冷气在上），自然对流翻腾
- $Ri = 0$：中性

本模型忽略剪切项（无详细速度梯度），简化为 $Ri \approx \Delta T/2$，并通过经验公式 $\mu(Ri)$ 控制层间混合强度。这是工程化的"次网格闭合"近似。

### A.6 风机相似律

离心/轴流风机的相似律：

$$
\frac{\dot V_2}{\dot V_1} = \frac{n_2}{n_1},\quad \frac{P_2}{P_1} = \left(\frac{n_2}{n_1}\right)^3
$$

风量 $\dot V$ 正比于转速 $n$，功率 $P$ 正比于 $n^3$。代码中 $f_r$ 等效于转速比，故风量按 $f_r$ 缩放（`:425`），功率按 $f_r^3$ 缩放（`:426`）。

---

*公式集完。共整理 14 大类约 60 条公式，含 6 处关键推导。所有公式均可在 KaTeX/MathJax 渲染器中正常显示。*
