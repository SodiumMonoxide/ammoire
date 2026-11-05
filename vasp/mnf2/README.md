# MnF₂ 的 VASP → AMMoire 对照

本算例检验无 SOC 的共线交错磁。用 VASP 计算两个自旋通道的电子能谱，检查节点和镜像关系，再与 AMMoire 的 `roig-mnf2` 四带模型比较。体相基准为后续材料 Wannier 化和界面建模提供校验。

## 输入来源

| 量 | 输入或来源 |
|---|---|
| 晶胞 | 金红石 6 原子，(a=b=4.873) Å，(c=3.311) Å，F 内坐标 (u=0.305)，取 Yuan 等文献中实验结构。 |
| 泛函和关联 | PBE、500 eV、Mn 3d 的 Liechtenstein (U=5) eV、(J=0) eV，关闭 SOC，依据 Yuan 等方法及 Fig. 2(b)。 |
| 网格 | `smoke` 为 3×3×4，`grid` 为 4×4×6，两者为本轮选择的 Γ 中心网格。`paper` 为文献的 10×10×14 Monkhorst–Pack 网格。 |
| 本机实施 | VASP 6.5.0 CPU 版，4 MPI ranks，PAW_PBE Mn_pv 02Aug2007 与 F 08Apr2002。赝势版本为本机选择。使用固定结构、`ISYM=-1`、`EDIFF=1e-7`、反向初始 Mn 磁矩。 |
| 探针 | 从各自 SCF 的 CHGCAR 出发，以 `ICHARG=11`、`LMAXMIX=4` 计算 17 个明列的分数倒格坐标。SCF 与探针使用同一结构、PAW、ENCUT 和 U/J，探针点列表与 SCF 积分网格各有用途。 |

[Yuan 等原文](https://doi.org/10.1103/PhysRevB.102.014422) · [公开全文](https://www.colorado.edu/faculty/zunger-matter-by-design/sites/default/files/attached-files/physrevb.102.014422.pdf) · [VASP LMAXMIX](https://vasp.at/wiki/LMAXMIX) · [Roig 等四带模型](https://doi.org/10.1103/PhysRevB.110.144412)

## 运行

在 WSL 或 Linux 中进入本目录，按本机授权安装设置以下变量。每次使用新的输出目录。脚本用临时目录保存 CHGCAR 及其他中间数据，计算结束移除临时拼接的 POTCAR。输出目录记录输入、文本输出和赝势标题及哈希。超时默认每阶段 480 秒。

```bash
export VASP_LAUNCHER=/path/to/vasp-cpu
export VASP_POTCAR_ROOT=/path/to/licensed/PBE
mkdir -p runs
bash ./run_wsl.sh runs/my-smoke smoke
```

`VASP_LAUNCHER` 接受两个位置参数，依次为 MPI ranks 和变体 `std`。已有 `vasp-cpu` 包装器时可直接使用。只有原生 `vasp_std` 时，创建本地包装脚本，以 `exec mpirun -np "$1" /path/to/vasp_std` 启动对应进程数。赝势目录下需有 `Mn_pv/POTCAR` 和 `F/POTCAR`。可用 `VASP_MN_PAW`、`VASP_RANKS` 和 `VASP_TIMEOUT_SECONDS` 覆盖默认值。

从仓库根目录分析新运行。仓库随附的是派生摘要和数值表。下例取各自旋的最高占据带号 27，适用于这里的 Mn_pv 计算。换赝势后先核实带号与占据。固定带号是谱秩，定量拟合需要用 `PROCAR` 和波函数重叠追踪，处理交叉和近简并。

```powershell
python vasp/mnf2/analyse.py vasp/mnf2/runs/my-smoke
python -m ammoire.vasp_bridge vasp/mnf2/runs/my-smoke/probe/EIGENVAL --band-up 27 --branch lower --param a=4.873 --param c=3.311 --out vasp/mnf2/runs/my-smoke/mnf2-vbm-vs-tb.csv
```

## 本机试算结果

两套 SCF 及各自 17 点的冻结密度探针均达到电子迭代阈值。Mn 的 PAW 投影磁矩均为 ±4.695 μB，晶胞净磁矩为 0（输出精度）。未施加 `NUPDOWN` 约束。

| 指标 | 3×3×4 | 4×4×6 |
|---|---:|---:|
| 6 个节点的占据谱最大自旋残差 | 3 μeV | 9 μeV |
| 镜像交换自旋后，整个占据谱的最大残差 | 3 μeV | 9 μeV |
| 所采样点最大占据谱秩劈裂的绝对值 | 655.509 meV | 655.431 meV |
| 固定带号 27 对四带下支的劈裂 RMSE | 51.388 meV | 51.384 meV |

两套网格之间，17 点 × 27 个占据谱秩的劈裂最大变化为 0.129 meV，第 27 带最大变化为 0.018 meV。最大谱秩劈裂位于分数坐标 (0.2, 0.2, 0) 的第 15 带。这些是采样点的谱秩指标。尚需检验 ENCUT、U 和更密的积分网格，才能评估材料数值误差。

以 4×4×6 为例，价带顶谱秩（第 27 带）的 (E_\uparrow-E_\downarrow) 在 (¼,¼,0) 约为 −38.6 meV，镜像后约为 +38.6 meV。四带模型给 −17.2 和 +17.2 meV，重现了变号。在一般点 (0.3,0.15,0.125)，两者分别约为 +3.95 和 −9.71 meV，符号也不同。第 27 带两 Mn 的 d 投影和约为 0.58–0.68，仍须跟踪轨道与波函数才能建立材料带的一一映射。

[劈裂对照图](signed_splitting.png) · [小网格结果](results-smoke/analysis/summary.json) · [加密网格结果](results-grid/analysis/summary.json) · [逐点带对照](results-grid/analysis/fixed_index_tb_comparison.csv) · [镜像和反演谱残差](results-grid/analysis/symmetry_pairs.csv)

## 对照原理

无 SOC 时，Roig 两亚晶格模型的两个自旋块具有质量项 $t_z(\mathbf{k})\pm J$，其中 $t_z=(t_6+t_7\cos k_zc)\sin(k_xa)\sin(k_ya)$。在 $k_x=0$ 或 $k_y=0$ 上，质量项平方相同。翻转 $k_x$ 或 $k_y$ 会交换两个自旋块的能谱。反演保留每个自旋块的能谱。比较 DFT 的自旋分辨谱即可检验这些关系。对单带劈裂作定量拟合时，需要先识别对应轨道。

将两自旋块的下支作差，可得

$$
\Delta E_-(\mathbf{k})=
\frac{-4Jt_z}{\sqrt{|t_x|^2+(t_z+J)^2}+\sqrt{|t_x|^2+(t_z-J)^2}}.
$$

其中 $t_x$ 是混合两个亚晶格的跃迁。分子要求交错交换 $J$ 与同亚晶格的跃迁差 $t_z$ 同时存在。靠近 Γ 时，$t_z$ 的最低阶是 $k_x k_y$，直接给出采样所见的节点和镜像变号。

WSL 上实际运行的 CPU 二进制编译选项中没有 `VASP2WANNIER90`，当前路径中也未找到 `wannier90.x`。因此本轮使用 `EIGENVAL` 和 `PROCAR`，材料 Wannier 输出需要带接口的 VASP 构建。后续有了同一规范的 `hr.dat`、层间矩阵和所需投影算符后，可接到 AMMoire 的莫尔哈密顿量。Hyperspin 材料还需共面非共线磁序、独立的物理自旋投影及 (J_y) 的构造和对易检验。[VASP 接口](https://vasp.at/wiki/LWANNIER90) · [VASP 非共线计算](https://vasp.at/wiki/LNONCOLLINEAR) · [Hyperspin 模型](https://doi.org/10.1103/dnr2-gv31)

## dsh 审阅

本地试算曾用受限 dsh 文本审阅核对自旋配对、磁矩和 CHGCAR 来源。固定电荷带计算可使用独立的探针点，节点和镜像关系根据原模型公式及实际能谱验收。[接力说明](../../HANDOFF.md) · [VASP ICHARG](https://vasp.at/wiki/ICHARG)
