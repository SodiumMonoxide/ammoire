# 从自旋 Laue 群到铁电畴间能谱

给定两个极化畴及其磁构型，哪一个空间操作能把初态映到终态，决定了共线交错磁体的自旋谱如何随切换改变。Sun 等建立的分类适用于自旋轨道耦合可忽略、畴间操作保持磁构型且自旋投影守恒的情形。此处的 I、II、III 标记**畴间切换路径**。常用的 d、g、i 波标记描述**单畴**自旋劈裂的动量节点。例见正文式（2）至（9）、图 2 和图 5，以及补充材料图 S1、S5 和 S7。[Sun et al., *Nature Communications* **17**, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2)

## 三种判据

以正畴的晶体 Laue 群 $G$ 和保持自旋子晶格的指数二子群 $H$ 为输入。每个生成元附带一个符号，$+1$ 保持自旋通道，$-1$ 交换通道。程序加入同自旋能谱的动量反演 $-I$，取带符号生成元的闭包，再检查畴操作 $O$ 的有效动量作用。这里的能谱反演也适用于实空间无反演中心的极性晶体。

| 畴操作 | 类型 | 对两个畴能谱的预测 |
|---|---|---|
| $O\in H$ | I 型 | 同一 $k$、同一自旋的能谱保持。 |
| $O\in G\setminus H$ | II 型 | 同一 $k$ 处两自旋通道的能谱互换，有符号劈裂反号。 |
| $O\notin G$ | III 型 | 自旋谱随 $O$ 在动量空间重新映射。群论允许个别动量、能带或张量分量保持。 |

所有矩阵采用同一笛卡尔坐标系，作用于列向量。对按行存储的采样点，正畴 $k$ 的伙伴是负畴的 `k @ O.T`。三个类型均满足 $E_{+,s}(k)=E_{-,s}(Ok)$。在严格二维层内，$M_z$ 对面内动量等效于恒等操作；在三维体相须使用完整矩阵。$C_{3z}$ 表示绕 $z$ 的 120° 主动旋转。复合操作从右向左作用。

## 文献路径和八带模型

```sh
python -m ammoire fe-response
python -m ammoire fe-response --bilayer --out results/fe-bilayer.json
```

第一条命令给出文献中八种“材料、磁序、路径”组合的分类、初末态及矩阵。双层 MnPS₃ 的路径 1 在 L1 磁序下对应 $M_z$（I 型），在 L2 下需要 $M_z M_y$（II 型）才能同时映射结构和磁构型。路径 2、3 对两种磁序均为 III 型。体相 BiFeO₃ 的两条路径分别由反演和绕 $[\bar{1}10]$ 轴的二重旋转描述，类型为 I 和 II。

第二条命令用 `make_model("sun-hex-bilayer")` 构造显式双层 TB。每层包含蜂窝晶格 A、B 两个磁位点，每个位点含两个物理自旋态。层内跃迁为 $t$，层间跃迁为 $0.5t$，局域交换幅度为 $J_{sd}S=t$。几何模型取三角 Bravais 晶格常数 $a$，层间面内位移为 $a/3$。层间近邻的选取和位移长度是对补充图 S1 的代表性实现。这里计算的是指定结构的能带；材料的极化、切换势垒和磁基态需独立确定。

算例从 `stacking=0, layer_sign=1` 出发。整数 `stacking` 指定位移方向（单位 60°），`layer_sign` 指定公共基底中交错磁矩的排列。切换后按新的堆叠重新构建跃迁。程序在相同实验室自旋轴上对角化两个畴。

| 操作 | 终态 `(stacking, layer_sign)` | 群论 | 相同动量的数值结果 |
|---|---|---|---|
| $M_z$ | (3, −1) | I | 保持 |
| $M_z M_y$ | (3, +1) | II | 自旋通道交换 |
| $M_z C_{3z}$ | (5, −1) | III | 通道保持和整体交换均无法描述所采样能谱 |

## 接入已有 TB 模型

```sh
python -m ammoire fe-response examples/fe_square_dxy.json
```

该输入使用已有的 `square-dxy` 模型，将两个畴的有效劈裂系数分别设为 $+0.2$、$-0.2$，固定 Néel 轴。镜像 $M_y$ 改变 $d_{xy}$ 谐波的符号，得到 II 型对照。默认参数下，四个采样点的交换通道残差为零，同通道残差约为 0.201 个模型能量单位。[完整输入](../examples/fe_square_dxy.json)

JSON 包含 `dimension`、带 `matrix` 和 `spin_sign` 的 `generators`，以及畴间 `operation`。这三项可独立用于群论分类。加入 `positive`、`negative` 和 `kpoints` 即可运行能谱检验。每个畴用注册模型名 `model` 和可选的 `parameters` 构造；`kpoints` 为笛卡尔动量，单位与该模型倒晶格一致。可选 `spin_axis` 默认为 $z$ 轴，`atol` 默认为 $10^{-8}$ 个模型能量单位。两个模型须采用共同的能量零点、长度与能量单位。

```python
import numpy as np
from ammoire import make_model
from ammoire.fe_response import domain_response
from ammoire.fe_symmetry import SpinLaueGroup

mirror_y = np.diag([1, -1])
c4z = [[0, -1], [1, 0]]
group = SpinLaueGroup(2, [(c4z, -1), (mirror_y, -1)])
positive = make_model("square-dxy", j=0.2)
negative = make_model("square-dxy", j=-0.2)
report = domain_response(group, mirror_y, positive, negative, [[0.43, 0.71], [1.02, -0.63]])
print(report["response_type"], report["sample_verdict"])
```

自定义模型应提供 `dimension`、`size`、批量 `hamiltonian(k)` 和同一正交基底中的三个物理自旋 Pauli 矩阵 `spin_operators`。检验沿指定轴投影到两个守恒自旋子空间，在每个通道内分别按能量排序，再比较全谱。带交叉处的轨道身份需波函数证据。HAM 的子系统标签不能代替物理自旋。分析非零 $k_z$ 体相截面时，应传入三维模型和三维动量。

## 报告判读

- `response_type` 取决于输入的完整群、磁构型和畴间操作。遗漏群生成元可能导致 III 型误判。群闭包和有限能谱样本均不能补全未知的晶体对称性。
- `positive_domain_group.max_residual` 检验正畴在全部带符号群操作下的能谱。`two_domain_spectra` 同时列出伙伴动量的同自旋残差、同一动量的保持与交换残差，以及两畴的最大劈裂。残差为全部采样点和通道内能量秩上的最大绝对差。
- `consistent_on_sample` 表示样本满足预测且能区分响应；`undistinguished_on_sample` 表示样本缺少区分度，例如只采到了劈裂节点；`inconsistent` 表示至少一项必要关系超出容差。材料判定还需独立确定两畴的结构及磁序，并检验一般动量处的色散。

已有两自旋通道的电导张量时，`ammoire.fe_reference.transform_conductivity(tensors, O)` 可计算另一畴中的 $O\sigma_s O^{\mathsf T}$。输入形状为 `(2,d,d)`，按自旋上、下排列；输出单位同输入。该变换对应在给定自旋分辨电导上应用空间对称性。绝对电导值需另由输运计算或测量给出。
