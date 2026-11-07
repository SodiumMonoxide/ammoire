# 先前转角研究的派生归档

这些文件来自 2026-11-06 已完成的本地 NiCl₂ 工作，在 2026-11-07 更新仓库时按文件校验和复制；本轮没有重新运行其 VASP、WAVECAR 展开或自旋螺旋计算。文件及原始研究内相对路径见 [source_manifest.json](source_manifest.json)。完整解释见 [转角说明](../../docs/twist_history.md)。

| 文件 | 含义 |
|---|---|
| `angle_scan.csv`、`tb_matched_q.csv` | 先前八角模型扫描和同 Q 的原胞/折叠支对照；本轮新代码重算的共同列已核对。 |
| `symmetry_checks.json` | 先前互补角、胞矢和谱关系检查。 |
| `same_momentum_edges.json`、`unfolded_sector_weights.csv`、`unfolding_checks.json` | 27.80° 实算波函数的层和原胞扇区派生投影及检验。 |
| `material_27_summary.json` | 27.80° 固定结构的恢复自洽、探针、加严对角化和投影摘要；首次不完整中断结果不作材料结论。 |
| `vasp21_convergence.csv` | 21.79° 截止能与自洽网格复核。 |
| `matched_static_moments.json` | 两角一致球半径下的平滑磁化密度积分。 |
| `spiral_summary.json` | 单层平面单 q 变分控制，双线性/双二次交换模型敏感性；不是双层磁相图。 |

只携带小型派生数值。完整原始计算和授权材料由用户在本地研究目录保管；材料重新分析需要对应原始文件。这里的赝波函数扇区范数、PAW 轨道权重和实空间磁化积分有不同定义，不能混作同一个投影。
