# FEMcode

用于科研与学习的 MATLAB 二维三角形有限元程序库。

[English](README.md)

FEMcode 提供网格拓扑、有限元空间、基函数求值、稀疏组装、边界处理、插值、误差计算和绘图。目前定位为研究代码库，包含可运行的 Poisson 示例和针对性的回归检查，尚未完成对所有元类型的全面验证。

## 环境与快速开始

- MATLAB。尚未确定最低支持版本，也未验证 GNU Octave 兼容性。
- 可选依赖：`mesh_circle_triangle` 和 `mesh_square_triangle_unstructured` 需要外部 DistMesh。仓库未附带 DistMesh；结构网格 Poisson 示例不需要它。

在仓库根目录启动 MATLAB，执行：

```matlab
addpath(pwd);
addpath(fullfile(pwd,'poisson'));
solver_2D_Poisson_d;
```

示例采用每个坐标方向 100 等分的结构三角网格和 P1 元，求解

$$
-\Delta u=-2e^{x+y}\quad\text{在 }(-1,1)^2\text{ 内},
\qquad u=e^{x+y}\quad\text{在边界上}.
$$

精确解为 `exp(x+y)`。程序输出网格尺寸、耗时，以及 L2、H1 和采样最大误差。问题数据见 `poisson/fun_*.m`。

## 已实现的有限元空间

以下标识符有对应的实现分支；列入此表不代表已经完成全面数值验证。

| 类型 | 标识符 |
| --- | --- |
| 分片多项式与加富标量空间 | `P0`、`P1`、`P2`、`bubbleP1` |
| 间断空间 | `P1dc`、`DGP1`、`DG-P2-quad` |
| Crouzeix–Raviart 元 | `CR` |
| Raviart–Thomas 元 | `RT0` |
| Bernardi–Raugel 与 Mardal–Tai–Winther 元 | `BR`、`MTW` |
| Nédélec 元 | `Ned1-1`、`Ned1-2`、`Ned2-1` |
| 研究变体及旧实现 | `CR-P0`、`CR-RT0`、`BR0`、`BR-RT0`、`RT0-MTW` |

标量空间可通过 `dim` 叠加为多分量空间。RT0、Nédélec 等本身为向量值的空间使用 `dim = 2`，向量基函数的各分量共享自由度。使用研究变体前应核对具体实现。目前未提供完整的四边形或三维计算流程。

## 核心流程

```matlab
[T,P] = mesh_square_triangle(0,1,0,1,16,16);
mesh = mesh_topo(T,P);
space = mesh_FE('P1',1,mesh,3,0,[]);
A = assemble_matrix_2D(@(x) ones(1,size(x,2)),[1,1], ...
    mesh,space,space,[1,1,1,0,1,0;1,1,0,1,0,1],3);
```

以上构造网格、标量 P1 空间和 Laplace 刚度矩阵。载荷组装、边界处理及求解过程见 Poisson 示例。

### 网格与自由度

`P` 为 `2 x N_node` 坐标矩阵，`T` 为 `3 x N_elem` 顶点编号矩阵。`mesh_topo` 统一三角形方向并建立连接关系。

| 字段 | 含义 |
| --- | --- |
| `mesh.E` | 四行依次记录两个相邻单元槽位和有向边的两个端点编号 |
| `mesh.TE` | 带符号的全局边编号；局部边顺序为 `(1,2)`、`(2,3)`、`(3,1)` |
| `mesh.ET` | 边在相邻单元中的局部编号，对应 `E(1:2,:)` |
| `mesh.tau`、`mesh.N` | 边的单位切向及其顺时针旋转得到的单位法向 |
| `mesh.s`、`mesh.hmax` | 边长及三角网格的最大边长 |

对边界边，`E(2,:)` 存放唯一相邻单元，`E(1,:)` 存放非正的边界标签。调用 `mesh_topo(T,P)` 时所有边界标签均为零。可通过 `mesh_topo(T,P,f0,f1,...)` 传入边界函数：当一条边的两个端点均满足某函数的零值判定时，赋予相应标签 `0,-1,...`。未匹配的边界边仍标记为零。标签用于选择边，本身不会施加边界条件。

有限元空间构造函数需要六个参数：

```matlab
space = mesh_FE(basis_type,dim,mesh,Gauss_type_2D,Gauss_type_1D,bndy_number);
```

- `Gauss_type_2D`：三角形积分点数，可选 `1`、`3`、`4`、`9`；同时用于缓存基函数值与一阶导数。
- `Gauss_type_1D`：边积分点数，可选 `1`–`5`；设为 `0` 时不生成边缓存。
- `bndy_number`：边缓存使用的边界标签；`[]` 表示包含所有边。

`space.N_node` 实际表示**自由度总数**。`N_lb` 为局部基函数数，`T` 将局部自由度映射到全局编号，`basis_type` 存放字符形式的类型标识符。`space.P` 根据元类型存放坐标或自由度元数据。向量基函数的不同分量可能共享全局编号。应使用 `get_Dbndynodes(mesh,space,label)` 获取边界自由度，不是每种空间都有 `Dbndynodes` 字段。

### 组装与缓存约定

`assemble_matrix_2D` 组装形如

$$
c\int_\Omega f\,(\partial_x^i\partial_y^j u_k)
(\partial_x^m\partial_y^n v_l)\,dx
$$

的项。`mat_info` 每一行为 `[k,l,i,j,m,n]`，`coe_num` 给出各项的常系数，试探空间与检验空间分别传入。`assemble_vector_2D` 组装载荷向量。`*_FE_*` 组装器还接受有限元函数作为系数，参数约定见源文件头部注释。

`treat_Dirichlet(A,b,mesh,space,label,boundary_fun)` 通过替换矩阵行施加给定自由度，通常不保持矩阵对称性。Poisson 示例使用 MATLAB 反斜杠直接求解。

**积分规则必须与基函数缓存一致。** 组装器和使用缓存的误差函数应采用构造空间时的三角形积分规则。若要用不同规则计算已有解的误差，先重建缓存：

```matlab
[space.basis,space.basis_dx,space.basis_dy] = ...
    generate_basis_data_triangle(mesh,space,9);
% 后续使用三角形缓存的计算必须采用规则 9。
```

## 误差计算与验证

- `compute_error` 支持 `L2`、`H1`、`H10`（H1 半范数）、`L_inf` 和 `Hdiv`。精确解接口为 `u_fun(points,dx,dy)`，每行对应一个分量。
- `compute_norm` 支持 `L2`、`H1`、`H1-semi`、`L_inf` 和 `Hdiv-semi`。
- `Hdiv` 误差为完整范数 `sqrt(||e||_L2^2 + ||div(e)||_L2^2)`；`Hdiv-semi` 为 `||div(u_h)||_L2`。两者要求二维向量场。导数逐单元计算，对非相容场得到的是分片量。
- `L_inf` 仅在积分点上采样，不是连续区域上的精确最大值。

在仓库根目录运行针对性的回归检查：

```matlab
addpath(pwd);
addpath(fullfile(pwd,'test'));
test_Hdiv;
```

测试使用向量 P1 场，覆盖散度抵消、非零散度、精确插值和常量误差的 L2 项。`test/` 中其余脚本为探索性脚本，尚未构成全面的回归测试集。

对上述 Poisson 问题单独进行过一次本地网格加密实验，保持三点积分组装，采用九点积分计算误差，结果如下：

| 每方向划分数 | L2 误差 | L2 阶 | H1 误差 | H1 阶 |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 2.08012e-02 | — | 3.71880e-01 | — |
| 16 | 5.18008e-03 | 2.0056 | 1.85299e-01 | 1.0050 |
| 32 | 1.29373e-03 | 2.0014 | 9.25683e-02 | 1.0013 |
| 64 | 3.23353e-04 | 2.0004 | 4.62740e-02 | 1.0003 |
| 128 | 8.08331e-05 | 2.0001 | 2.31357e-02 | 1.0001 |

收敛阶按 `log2(E_N/E_2N)` 计算。仓库尚未附带该实验的自动化收敛测试脚本。默认示例使用三点积分计算误差，因此误差数值不同。该结果验证了这一光滑 Poisson 算例，不代表所有元类型和边界处理均已验证。仓库未提供 Stokes、磁扩散或 ALE 示例，也未配置 CI 工作流。

## 文件导航

| 文件 | 用途 |
| --- | --- |
| `mesh_*.m`、`orientation_consistence.m` | 网格生成、拓扑和空间构造 |
| `basis_*.m`、`generate_basis_data_*.m` | 基函数求值与缓存 |
| `generate_Gauss_*.m`、`local_*.m`、`assemble_*.m` | 积分与组装 |
| `compute_dof.m`、`interpolate.m`、`FE_function_*.m` | 自由度、插值与函数求值 |
| `treat_Dirichlet*.m`、`get_Dbndynodes.m` | 边界处理 |
| `compute_error.m`、`compute_norm.m`、`plot_*.m` | 诊断与绘图 |
| `poisson/`、`test/` | 示例问题与检查 |

## 来源与许可证

程序雏形根据密苏里科技大学何晓明老师的[有限元编程课](https://www.bilibili.com/video/BV1Zv411t7Lj)编写，之后进行了修改与扩展。

仓库附带 [GNU 通用公共许可证第 3 版](LICENSE.md)。
