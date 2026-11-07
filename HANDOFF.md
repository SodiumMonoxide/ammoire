# 交给本地 dsh 的接力说明

## 2026-11-07 发布接力

本轮任务是发布更新后的 **AMMoire 0.2.0**。本机源目录、仓库名和 Pages 选择由仓库外的 `DSH_UPLOAD.md` 提供。先使用该文件中的发布任务；下面原有材料研究接力保留作后续研究记录。

上传前阅读 [更新说明](docs/update_2026-11-07.md)、[转角证据说明](docs/twist_history.md) 和 [发布步骤](RELEASING.md)。新增公度求解器吸收后续多角度研究；`benchmark/twist_history` 是此前计算的派生归档，`benchmark/commensurate` 与 `benchmark/cairo` 是本轮重算。VASP、波函数展开和自旋螺旋的材料结果在本轮仅复核归档，不宣称重跑。

先运行以下新增检查，再运行下方原有基线：

```sh
python -m examples.commensurate_demo --points 8 --out results/dsh-commensurate
python -m ammoire.cairo_demo --points 12 --out results/dsh-cairo
```

25 个原生单层预设与公度 TB 均不依赖旧日报目录。`legacy.daily_model` 仍是可选的外部适配器。公度单轨道电子模型固定磁序，不能代替材料轨道拟合、磁序优化或 VASP 波函数展开。

## 接力入口

在本地 dsh 客户端打开上传后克隆的仓库目录。把以下文字作为新会话的第一条任务。

> 先读 AGENTS.md、HANDOFF.md 和 README.md。运行 NumPy、VASP 解析及 Wannier90 解析基线。以 MnF₂ 的已验证节点和镜像谱关系为起点，检查本地 VASP/Wannier90 能力，准备带轨道依据的 DFT→TB 对照。先报告可复核的基线和阻塞项，再依据 HANDOFF.md 的条件推进材料计算。每步记录命令、证据和剩余不确定性，输出到新的本地目录。

## 可复核的当前状态

| 基线 | 已验证结果 |
| --- | --- |
| AMMoire | 24 个内置 TB 预设，平面波莫尔双层、层间及层内傅里叶矩阵、Wannier90 `*_hr.dat` 读入。 |
| HAM／AM 示例 | 4° 扭角、2 壳、13 个平面波、78 条带。同扇区跃迁时 `[H,Q]=0`。加入跨扇区跃迁后出现混合。2 壳与 3 壳的 Γ、X、M 中心谱权最大相对差约 0.29%，对应 0.04 的高斯展宽。 |
| MnF₂ 实算 | 六原子金红石原胞，VASP 6.5.0，PBE+U，关闭 SOC。3×3×4 和 4×4×6 的 SCF 与各 17 点冻结密度探针均达到电子阈值。Mn 磁矩 ±4.695 μB，节点和镜像交换自旋的残差不超过 9 μeV。两网格间占据谱劈裂最大变化 0.129 meV。 |
| MnF₂ 对照边界 | Roig 四带模型复现节点、镜像换自旋与反演同自旋关系。固定第 27 带对下支的劈裂 RMSE 约 51 meV，一般三维动量处出现符号差异。第 27 带是最高占据谱秩，轨道对应仍待确认。 |
| Wannier90 接口 | 上次实际运行的 VASP CPU 构建未启用接口。解析器已经通过合成 `hr.dat` 自测，材料插值仍需实测。 |

## 先跑的命令

在仓库根目录使用 Python 3.10+ 和 NumPy。编辑安装为可选步骤。示例会写入指定输出目录。

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
python -m ammoire.vasp_bridge --self-test
python -m ammoire.wannier
python -m ammoire models
python -m ammoire demo --angle 4 --shells 2 --check-cutoff --out results/dsh-baseline
```

历史基准可查 [莫尔摘要](benchmark/summary.json) 和 [MnF₂ 输入及结果](vasp/mnf2/README.md)。仓库携带派生的公开基准，轨道追踪需要用户已有运行目录中的 `PROCAR` 或本地重新生成的文件。现有 VASP 运行目录可用 `python vasp/mnf2/analyse.py <run-directory> --out <new-analysis-directory>` 重新分析。

## 后续工作与验收

1. **优先对齐轨道。** 用 `PROCAR` 的离子及轨道投影筛选目标子空间，在交叉或近简并点用波函数重叠辅助追踪。保留同一 (k) 上的上下自旋、镜像搭档和晶胞约定。记录带号变化、投影权重、节点与镜像残差，并将拟合点与留出的验证点分开。投影总权重和固定谱秩都只能作线索。
2. **核实本地 VASP 构建。** 从实际运行的 VASP 二进制读取 `--cpp-options`，检查版本、Wannier90 可执行文件和授权 PAW 可用性。只有实际编译选项含 `VASP2WANNIER90` 或 `VASP2WANNIER90v2` 且 Wannier90 版本相容时，才安排材料导出。若缺接口，先给出独立构建计划和小型接口测试，等用户同意再改安装或启动重编译。
3. **先做短算例。** 运行前明确资源预算和运行目录。新增 SCF 基准先核实结构、两 Mn 的 AFM 初态与最终磁矩、PAW 顺序、ENCUT、U/J 和版本。两个网格的比较只覆盖本轮点位，数值收敛另做加密网格与 ENCUT 对照，U 的影响另报告。用该 SCF 的 CHGCAR 进行固定密度探针，保持相同晶胞和物理参数，检查 `LMAXMIX`、退出码、EDIFF 和本征值文件。使用论文网格前先确认计算预算。
4. **验收材料 `hr.dat` 后接莫尔。** 记录轨道投影、能窗、截断、长度与能量单位、Wannier 规范及几何。用独立点比较 Wannier 与 DFT 的能量误差和劈裂，验证节点和镜像。导入器要求常规简并权重及兼容的 Wigner–Seitz 导出方式，详见 [README.md](README.md)。用一致规范提供层间矩阵，再检查莫尔截断。HAM 方向另需共面磁序依据、物理自旋投影及所选 (J_y) 的对易检验。

如已具备计算条件，在 WSL 中按 [MnF₂ 说明](vasp/mnf2/README.md) 设置本机 `VASP_LAUNCHER` 和 `VASP_POTCAR_ROOT`，先用 `smoke` 模式测试。本地接力记录宜包含 Git 提交号、运行命令、退出状态、各阶段结果链接、下一项最小工作。使用 dsh 云端审阅时，仅传公开来源与任务必要的非敏感摘要。
