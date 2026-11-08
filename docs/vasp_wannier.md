# 从 VASP 导出 Wannier90 紧束缚模型

AMMoire 从同一 seed 的 `.win`、`_hr.dat` 和 `_wsvec.dat` 读入晶胞、跃迁及轨道对距离修正。可以先与 Wannier90 自己的参考路径核对插值，再用独立 DFT 点检验材料参数化。[Wannier90 文件格式](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/) · [插值约定](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/notes_interpolations/)

## WSL 接口验收

2026-11-08 已在 Ubuntu 24.04 WSL2 核验 VASP 6.5.0 CPU/GPU 的 `std`、`gam` 和 `ncl` 均含 `-DVASP2WANNIER90`。两种后端此前均通过官方 MoS₂ 的无 SOC 与 SOC 插值对照。独立 Wannier90 3.1.0 读取各自的 VASP 导出矩阵完成运行。GPU 使用安装时配置的专用启动器和 MPI 运行时。

本轮用两个 CPU 算例已导出的矩阵在独立目录重跑 `wannier90.x`。均采用默认的轨道对 Wigner–Seitz 修正。AMMoire 比较同一 seed 的 29 个路径点，判据为最大误差不超过 0.02 meV。

| 输入 | Wannier 带数 | RMSE (meV) | 最大误差 (meV) |
|---|---:|---:|---:|
| MoS₂ 无 SOC | 13 | 0.003065 | 0.01030 |
| MoS₂ SOC 自旋量子态 | 26 | 0.004766 | 0.01425 |

这两组数值衡量 AMMoire 与同一 Wannier90 插值的差异，包括文本文件的截断精度。材料对 DFT 的拟合误差需另用独立 DFT 点计算。[本轮验证记录](../validation/release_checks_0.3.0.json)

## 导出均匀网格和参考路径

1. 为材料建立独立运行目录。使用完整均匀 k 网格导出，设 `ISYM=-1`，核对实际 `.win` 的 `mp_grid`、`kpoints` 和 VASP 输出的逐点顺序。若从同一结构和电子参数的自洽密度做固定电荷计算，使用匹配的 `CHGCAR` 和 `ICHARG=11`。含 d 轨道的 DFT+U 对照需要在自洽及固定电荷计算中使用兼容的 `LMAXMIX`，通常为 4。[ISYM](https://vasp.at/wiki/ISYM) · [ICHARG](https://vasp.at/wiki/ICHARG)
2. 导出给独立程序时，设置 `LWANNIER90=.TRUE.`、`LWRITE_MMN_AMN=.TRUE.`，保持 `LWANNIER90_RUN=.FALSE.`。采用库内运行时，打开 `LWANNIER90_RUN`。核对实际 `.win`、`.amn`、`.mmn`、`.eig` 里的网格、带数、投影和能窗。共线 `ISPIN=2` 时分别处理两个自旋 seed，按输出核对 `.1`、`.2` 的对应关系。SOC 使用自旋量子态基底，物理自旋期望需要投影到同一 Wannier 规范的自旋算符。[LWANNIER90](https://vasp.at/wiki/LWANNIER90) · [LWANNIER90_RUN](https://vasp.at/wiki/LWANNIER90_RUN) · [LWRITE_SPN](https://vasp.at/wiki/LWRITE_SPN)
3. 对每个 seed 在 `.win` 设置 `write_hr = true`、`bands_plot = true`、`bands_plot_mode = s-k` 和 `kpoint_path`。保留默认 `use_ws_distance = true` 及对应 `_wsvec.dat`。用 `wannier90.x SEED` 生成 `_hr.dat`、`_wsvec.dat`、`_band.kpt`、`_band.dat`，从 `.wout` 核对 Wannier 化和插值确实完成。`write_ndegen_applied` 由较新版 Wannier90 提供，3.1.0 使用原始旁文件格式。[Wannier90 参数](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/parameters/) · [Wannier90 4.0.3 变更记录](https://github.com/wannier-developers/wannier90/blob/v4.0.3/CHANGELOG.md)

## 导入和验收命令

```bash
python -m ammoire.wannier_validate path/to/wannier90.1 --tolerance 2e-5 --out seed1_validation.json
```

传入实际 seed 前缀，包括可能存在的 `.1` 或 `.2`，随后对另一通道单独运行。阈值和报告均使用 eV，默认阈值为 `1e-4` eV。命令按参考路径的分数倒格坐标转换成笛卡尔坐标，按同一能量零点比较谱，并检查 (H(k+G)=H(k))。超过阈值返回退出码 1。不同的 k 点或能量打印精度可能需要不同阈值。[参考能带格式](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_banddat)

```python
import numpy as np
from ammoire.wannier import WannierModel

model = WannierModel.from_seed("path/to/wannier90.1")
q = np.array([0.12, 0.18, 0.0])  # fractional reciprocal coordinates
energies_eV = np.linalg.eigvalsh(model.hamiltonian(q @ model.reciprocal))
```

`from_seed` 使用 `.win` 的三维行晶格，接收 Å 或 Bohr，并为默认距离修正要求匹配的旁文件。实现将每个轨道对的 (H_{mn}(R)) 按原始 (R) 的简并和该对轨道的等价平移数分别归一，累加到 (R+T)，再检查 Hermitian 关系。旧 `from_hr(path, lattice, dimension=...)` 继续支持显式晶格。只有 `hr.dat` 时需由调用者确认导出模式。二维输入要求有效跃迁的 (R_z=0)。[Wannier90 3.1.0 插值实现](https://github.com/wannier-developers/wannier90/blob/v3.1.0/src/plot.F90)

材料接力时，先固定 DFT 与 Wannier 的共同晶胞、参考能量和自旋通道，再在训练网格外选 k 点核对目标能窗、轨道权重、劈裂及节点。带交叉或折叠带的比较需要波函数或轨道子空间匹配。组合莫尔层还需要一致规范的层间跃迁矩阵。
