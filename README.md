# AMMoire

交错磁电子紧束缚模型、转角双层连续模型与 hyperspin 守恒诊断。代码用 Python 和 NumPy，绘图另需 Matplotlib。

AMMoire combines electronic tight-binding models with a plane-wave moiré bilayer solver and hyperspin diagnostics. It also compares collinear VASP eigenvalues with bulk tight-binding models.

- [电子 TB 文献回顾](literature_core.md)覆盖对称性模型、微观多亚晶格模型、材料 Wannier、SOC／拓扑和配对变体，分列已发表、录用与预印本。
- [日报、转角模型和 hyperspin 对照](history_and_moire.md)汇总 4 份文章简报和 3 天模型日报，记录已复用的原代码及边界。
- [默认算例能带](benchmark/bands.png)、[数值结果](benchmark/summary.json)和[逐带原始数据](benchmark/bands.csv)。
- [MnF₂ 的 WSL VASP 基准](vasp/mnf2/README.md)及共线 `EIGENVAL` 到 TB 劈裂的对照入口。
- [本地 dsh 接力](HANDOFF.md) · [仓库协作约定](AGENTS.md) · [GitHub 发布步骤](RELEASING.md)

## 运行

在仓库根目录使用 Python 3.10+。绘图依赖由 `plot` 选项安装。安装后也可直接运行 `ammoire models` 和 `ammoire demo`。

```powershell
python -m pip install -e ".[plot]"
python -m ammoire models
python -m ammoire demo --angle 4 --shells 2 --check-cutoff --plot --out results/quickstart
python -m unittest discover -s tests -v
python -m ammoire.wannier
python -m ammoire.vasp_bridge --self-test
```

已有 NumPy 时，可直接运行列表、测试和不带 `--plot` 的 demo。示例输出 `bands.csv` 和 `summary.json`。输入 `--angle 0 --shells 0` 给出均匀堆叠参考，路径用原晶格布里渊区。非零扭角下 `--shells 0` 保留 (G=0)，用于检查动量转移关闭后的极限。偏压定义为顶层减底层能量。角度单位为度，示例能量和长度使用共同的模型单位。

## 模型与计算接口

| 用途 | API 或注册名 |
|---|---|
| 单带 AM、节点与角向对照 | `square-dxy`、`square-dx2-y2`、`square-g`、`hex-i`、`hex-g-3d`、`cubic-i` |
| 轨道分辨的 AM | `lieb-dx2-y2`、`brekke-lieb`、`gondolf-dxy`、九个 Gondolf 层群代表模型、六个 Roig 模型或材料参数组 |
| HAM | `SquareHAM()` 与 `projected_block(h0,h1,h2,h3,h4,h5)`。前者是方格演示参数，后者保留 Ma 等论文式 (6) 的四带结构。 |
| 旧日报模型 | `ammoire.legacy.daily_model(name, **parameters)`。支持 `roig-2d`、`farajollahpour`、`orbital-hall-6`、`wen-kee`、`guo-90`、`gese-90`。 |
| Wannier90 | `ammoire.wannier.WannierModel.from_hr(path, lattice, dimension=2)` |
| 转角与叠层 | `MoireBilayer(top, bottom, twist_deg=..., shells=..., harmonics=..., bias=...)` |
| 能带与算符 | `observables.layer_operator` 提升层内算符，`commutator_norm` 计算残差，`resolve_bands` 在简并子空间中对角化投影后的观测算符。 |

以下构造两个异向 Néel 轴、相同方格原晶格的层。数值矩阵顺序为 `(层, G, 轨道, 自旋或子系统标签)`。单层和双层 `hamiltonian(k)` 都接收笛卡尔动量，可批量输入 `(..., 2)`。三维单层入口用 `(..., 3)`，`slice_kz(kz)` 返回固定体相动量截面。

```python
import numpy as np
from ammoire import MoireBilayer, make_model
from ammoire.moire import tunneling_harmonics

top = make_model("square-dxy", j=0.18, neel=(0, 0, 1))
bottom = top.with_n((0, 0, -1))
harmonics = tunneling_harmonics(top.reciprocal, constant=0.04, first_star=0.02)
stack = MoireBilayer(top, bottom, twist_deg=3.0, shells=2, harmonics=harmonics)
labels, vertices = stack.high_symmetry_path()
energies = np.linalg.eigvalsh(stack.hamiltonian(vertices))
```

`MODEL_REGISTRY[name]` 记录原论文、采用的方程和默认参数性质。代码同时包含直接实现的矩阵和对称性代表参数化。比如 Gondolf 的 Appendix A 参数组可作原文四带对照，其余九个层群入口保留表 II 的主要谐波，自旋无关色散选用代表参数。RuO₂ 等入口沿用来源论文的模型参数，磁序由输入指定。体相固定 (k_z) 截面保留体相自由度，可作数值对照，材料薄膜需要单独构造。

旧日报适配器是可选入口。若持有旧项目，将 `ALTERMAGNET_CONTINUATION` 指向包含 `tb_extensions` 的续接目录，将 `ALTERMAGNET_TB_ROOT` 指向旧模型库根目录。内置预设、测试、莫尔算例和 VASP/Wannier 解析独立运行。适配器为 Roig-2D／Wen–Kee 应用保持谱一致的周期 Bloch 规范变换。组合层前需统一能标、晶格和轨道规范。原模型单位见 `ammoire/legacy.py`，可用 `energy_scale` 和 `energy_unit` 转换。GeSe-90 与 Guo-90 是整套双层基准。

Wannier90 导入接口接收以 Å 为单位的行晶格矩阵，输出 eV。二维输入要求每条跃迁的第三个平移分量均为零。支持常规 `*_hr.dat` 及简并权重，导出时使用 `use_ws_distance=false`，或使用已把距离修正展开的 `write_ndegen_applied=true`。有旁文件时会检查 `*_wsvec.dat` 的标志。自测使用合成文件覆盖复数跃迁、权重、倒格周期和不合法输入。材料层间跃迁及自旋投影仍需由同一 Wannier 规范提供。[Wannier90 文件说明](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat)

## 莫尔近似和 hyperspin 判据

连续模型使用同一母晶格的两个层，各旋转 $\pm\theta/2$。平面波 $G$ 位于莫尔倒格点，层间矩阵满足 $H_{tG,bG'}=T_{G-G'}$，反向块取伴随。可传任意复数矩阵，包括层间不同轨道数、不同自旋通道和第一星各自的相位。方形、矩形和六角母晶格可用。锚点可以设在 Γ 或指定母晶格谷。物理自旋使用实验室坐标，若原模型定义在随层旋转的自旋坐标，可显式打开 `local_spin_frame` 并以实验室基底提供 $T$。HAM 的内层标签不适用物理自旋旋转。

`potential_top`、`potential_bottom` 可输入层内莫尔标量势或矩阵交换调制，键为倒格基底的整数坐标。给出一个 $U_g$ 时自动补齐 $U_{-g}=U_g^\dagger$。若采用层内调制，检查 $Q$ 守恒时也要检查各层 $[U_g,Q_\ell]$。

设 $Q=\operatorname{diag}(I_G\otimes Q_t,I_G\otimes Q_b)$，其中顶层 $Q_t=2J_y/\hbar$，底层 $Q_b=\sigma_z$。层内各自守恒时，双层守恒的充要条件是每个被截断基底实际耦合的谐波满足 $Q_tT_g=T_gQ_b$。`charge_channels` 将输入的 $T_g$ 分成 $(T_g+Q_tT_gQ_b)/2$ 和 $(T_g-Q_tT_gQ_b)/2$，分别连接相同、相反 $Q$ 扇区。图中颜色为 ⟨Q⟩，点大小由 $G=0$ 权重决定。CSV 按每个 $k$ 的能量排序，`band` 是该点的排序序号。通过简并和交叉追踪连续带时需要波函数重叠。

4° 基准的 2 壳模型含 13 个 $G$ 和 78 条带。取 HAM 能量偏移 1.17，将其上支放到 AM 带边附近。只含同扇区跃迁时，路径上 $\max\|[H,Q]\|_F=0$。加入相对幅度 0.7 的跨扇区跃迁后，该值约为 0.804（模型能量单位），中心 $G$ 谱权加权的 $Q$ 方差约为 0.468。截断检查比较 2 壳与 3 壳在 Γ、X、M 的高斯展宽 $G=0$ 谱，宽度 0.04，最大相对 $L_1$ 变化约为 0.29%。这个误差绑定到该展宽、顶点和参数。Frobenius 范数随维数变化，比较不同截断时优先用同一谱权或归一化对易指标。

目前默认层间跃迁为设定的傅里叶系数，磁序固定。动量截断可在迷你布里渊区边界破坏精确周期性，过大的 (G) 也会超出所选谷的有效范围。材料预测还需要校准层间矩阵、层内莫尔势和磁织构。物理自旋、输运和拓扑分析还需要相应的投影算符、规范匹配和收敛验证。相关方法参照 [Bistritzer–MacDonald](https://doi.org/10.1073/pnas.1108174108)、[Balents](https://scipost.org/10.21468/SciPostPhys.7.4.048) 与 [Ma 等](https://doi.org/10.1103/dnr2-gv31)。
