# 转角材料的共线谱对称伙伴检查

`ammoire.collinear_symmetry` 将静态 `ISPIN=2` 的 `EIGENVAL` 按用户指定的笛卡尔操作配对。报告分别列出换自旋或同自旋的伙伴残差、指定节点的同 k 自旋差，以及选定能量秩的带符号劈裂。占据谱和全谱分别统计，便于判断高空带收敛是否影响特定带边的结论。

输入探针应来自同一晶胞、固定磁结构和可比的能量基准。晶格用 Å 的三个实空间行矢。坐标从每份 `EIGENVAL` 读入，经对应倒格基矢转换。对每个配对，代码先验证 `k_target = k_source @ operation_cartesian.T + reciprocal_shift @ B`，再比较能谱。倒格平移为可选的三个整数。操作与换自旋规则由结构和磁序给定。

下例以 (k_y) 反射为例，假定 `probes/EIGENVAL` 的前 3 个 k 点依次是 (Q)、镜像 (Q')、Γ。按自己的计算替换晶格、配对和每自旋占据带数。单点文件可省略 `k_index`，多点文件用一基编号。

```json
{
  "lattice_angstrom": [[4.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 20.0]],
  "occupied_bands": 12,
  "probes": {
    "Q": {"eigenval": "probes/EIGENVAL", "k_index": 1},
    "Q_reflected": {"eigenval": "probes/EIGENVAL", "k_index": 2},
    "Gamma": {"eigenval": "probes/EIGENVAL", "k_index": 3}
  },
  "pairs": [
    {"source": "Q", "target": "Q_reflected", "operation_cartesian": [[1,0,0],[0,-1,0],[0,0,1]], "swap_spin": true},
    {"source": "Gamma", "target": "Gamma", "operation_cartesian": [[1,0,0],[0,-1,0],[0,0,1]], "swap_spin": true}
  ],
  "nodes": ["Gamma"],
  "bands": [12, 13]
}
```

```bash
python -m ammoire.collinear_symmetry probe_manifest.json --out probe_report.json
```

`max_occupied` 是每个自旋内指定占据前缀的最大绝对残差，`max_all` 使用全部输出秩。报告也提供最大残差对应的秩、符号和自旋，以及 `EIGENVAL` 中占据前缀的最小占据和其余带的最大占据。`selected_rank_splitting_eV` 的符号为 (E_\uparrow-E_\downarrow)。对部分占据或两个自旋占据数不同的体系，固定前缀只是用户定义的比较窗口。需要追踪物理轨道或原胞动量身份时，将此检查和投影或展开结合使用。

本入口由近期转角多材料研究中的伙伴对照提炼。相同代码已用本机七点归档逐数回放，仓库测试使用可公开复现的合成谱。
