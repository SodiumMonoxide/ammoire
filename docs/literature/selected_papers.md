# 交错磁与 moiré 文献筛选及 TB 更新依据

核验日期：2026-11-07（Asia/Shanghai）。本报告精选 **14 篇一手研究**；它记录文献和更新依据，功能完成情况以仓库变更与测试记录为准。

适合当前工具的更新是把既有 Cairo 日报模型整理为可独立调用的预设，并增加可核验的公度双层几何。Liu/Zeng 提供转角与堆叠的物理分类依据；整数超胞构造和计数须由代码独立验证。几何公度只保证周期性，不能独自证明交错磁性，也不提供材料的跃迁参数。

## 查重边界

已对照旧 TB 索引、续接日报索引（截至 11/07，共 117 条）及 AMMoire 的 `literature_core.md`、`history_and_moire.md`。同时比较 DOI、arXiv 别名和题名；部分旧题名截断，优先用 DOI 关联。完整条目与阅读深度见 `selected_papers.json`。

本次新补查、且在这些记录中未命中的只有 **CrO（2608.01235）** 与 **RbV₂Se₂O 纹理锁定（2604.18337）** 两篇。这里的“新增”只指相对所查记录新增，不表示用户此前从未读过。其余文献是既有模型、候选、仓库文档或此前转角项目已使用的资料。Sheng 的论文虽未在 117 条日报索引中命中，已见于仓库文档，仍归作既有资料。

## 可直接用于本轮更新的来源

### 1. Altermagnetism and strain induced altermagnetic transition in Cairo pentagonal monolayer

Shuyi Li 等，2025。npj Quantum Materials 10, 83 (2025)；既有日报实现。[原始来源](https://doi.org/10.1038/s41535-025-00793-0) · [arXiv 2412.16857](https://arxiv.org/abs/2412.16857)。

Cairo 六位点模型保留两磁性、四非磁性位点；对角应变改变 g/d 波节点。可将当日日报模型整理为独立预设，并校验未应变节点及应变后的低阶分量。

实现边界：单层电子 TB；可作为扭转前的微观输入，论文自身并非 moiré 超胞。 公式入口：arXiv v1 supplement S5-S7。

### 2. Twisted Magnetic Van der Waals Bilayers: An Ideal Platform for Altermagnetism

Yichen Liu、Junxi Yu、Cheng-Cheng Liu，2024。Physical Review Letters 133, 206702 (2024)；既有记录／仓库资料。[原始来源](https://doi.org/10.1103/PhysRevLett.133.206702) · [arXiv 2404.17146](https://arxiv.org/abs/2404.17146)。

通过翻层和相反转角保持交换反向磁层的面内二重旋转；单层对称性约束 d/g/i 波类型。适合作为公度几何和对称性检查的物理来源。

实现边界：公度格点几何必须另行构造并验证；通用整数超胞工具不等同复现文中 Slater–Koster 跃迁或 VOBr DFT。 公式入口：arXiv v1 Eqs.1-3。

### 3. Bilayer stacking A-type altermagnet: A general approach to generating two-dimensional altermagnetism

Sike Zeng、Yu-Jun Zhao，2024。Physical Review B 110, 174410 (2024)；既有记录／仓库资料。[原始来源](https://doi.org/10.1103/PhysRevB.110.174410) · [arXiv 2407.15097](https://arxiv.org/abs/2407.15097)。

对 A 型双层交错磁给出层群堆叠分类和可交换两层的操作。公度与交错磁是两个独立条件，生成超胞之后仍要核对磁序和层间对称性。

实现边界：理论分类指导几何约束；未提取该文完整材料跃迁参数。

### 4. Valley polarization in twisted altermagnetism

San-Dong Guo 等，2024。Physical Review B 110, L220402 (2024)；既有日报实现。[原始来源](https://doi.org/10.1103/PhysRevB.110.L220402) · [arXiv 2406.13950](https://arxiv.org/abs/2406.13950)。

层偏压使不同层的相反自旋谷错开，原文给出可忽略层间耦合的极限。已有日报实现应保留来源和近似标签，避免重新计作新模型。

实现边界：T=0 解耦四带基准；任意转角的旋转色散不能自动解释为完整 moiré 折叠和散射。 公式入口：arXiv v1 Eqs.1-4。

## 交叉方向的材料与机制参照

### 5. Sublayer-resolved altermagnetism in twisted antiferromagnetic bilayers with spin-layer coupling

Kehan Liu 等，2026。Physical Review B 113, 235429 (2026)；既有候选。[原始来源](https://doi.org/10.1103/ln1f-9x8r)。

MnSe 双层的六角翘曲色散与亚层相关 moiré 势共同产生 i 波劈裂；AA 与 AAᴿ 堆叠的对称机制不同。适合后续亚层投影和势相位约束基准。

实现边界：需先获取全文/补充中的系数与规范；本次不将摘要所述连续模型当成已复现。

### 6. Ubiquitous van der Waals altermagnetism with sliding/moire ferroelectricity

Yuxuan Sheng 等，2025。Science China Physics, Mechanics & Astronomy 68, 297511 (2025)；既有记录／仓库资料。[原始来源](https://doi.org/10.1007/s11433-025-2698-5) · [arXiv 2411.17493](https://arxiv.org/abs/2411.17493)。

滑移或扭转堆叠将铁电极化和自旋劈裂联系起来，可为层势与交换调制选题提供依据。应区分作者所称 pseudo-altermagnetism 与由旋转对称性补偿的交错磁。

实现边界：材料设计参照；尚无本次提取的完整电子跃迁参数。

### 7. Moiré-induced altermagnetism from nonmagnetic constituents

Jingtian Shi、Maxim Khodas、Ivar Martin，2026。预印本（本次未核实到正式期刊）；既有候选。[原始来源](https://arxiv.org/abs/2609.03271) · [arXiv 2609.03271](https://arxiv.org/abs/2609.03271)。

非磁方格层经结构弛豫、交替畴壁旋涡及电子相互作用产生 d/g 波。它要求弹性弛豫与自洽 Hartree–Fock，不能用预先指定的交换场代表自发磁序形成。

实现边界：研究路线；当前固定哈密顿量工具不覆盖完整自洽形成机制。 公式入口：arXiv v2 End Matter S6-S11。

### 8. Altermagnetism in twisted van der Waals homostructures

Junying Chen 等，2026。Nature Physics (2026)；既有记录／仓库资料。[原始来源](https://doi.org/10.1038/s41567-026-03440-y)。

正交 CrPS₄ 同质结构的磁光、Raman 与 DFT 证据构成实验参照。实验厚度和 90° 几何应明确记录，不能把小角连续模型能带当成这些光谱的直接复现。

实现边界：实验背景与材料选取；无本次新增的拟合电子模型。

## 候选与独立研究分支

### 9. Altermagnetic Magnons in Twisted van der Waals Antiferromagnets

Qirui Cui 等，2026。Nano Letters 26, 5078–5085 (2026)；既有记录／仓库资料。[原始来源](https://doi.org/10.1021/acs.nanolett.6c00198)。

转角改变层内交换分布，在 CrPS₄/CrI₃ 产生交错磁磁振子劈裂。它适合未来独立自旋波分支；电子费米子 TB 的本征求解不能直接替代玻色 BdG 问题。

实现边界：玻色磁振子，非当前电子 TB 新预设。 公式入口：spin Hamiltonian Eq.1。

### 10. Twist-induced magnetic topological phase transition in stacked altermagnetic CrO

Zi-Hao Ding 等，2026。预印本（本次未核实到正式期刊）；本次相对所查记录新增。[原始来源](https://arxiv.org/abs/2608.01235) · [arXiv 2608.01235](https://arxiv.org/abs/2608.01235)。

公度扭转 CrO 的磁性与拓扑相变给出通用动量处的 Weyl 点。建议后续导入可核验 Wannier 数据并做全 BZ 节点搜索；本次仅补充候选记录。

实现边界：材料 DFT/拓扑研究；当前通用几何与谐波不是该材料的复现。

### 11. Direct observation of quadruple spin-texture locking in a 2D d-wave altermagnet

Dan Mu 等，2026。预印本（本次未核实到正式期刊）；本次相对所查记录新增。[原始来源](https://arxiv.org/abs/2604.18337) · [arXiv 2604.18337](https://arxiv.org/abs/2604.18337)。

RbV₂Se₂O 自旋 STM/QPI 报告自旋纹理锁定，并将条纹归因于自旋密度波 moiré。它扩展了观察方向，但该 moiré 不是双层几何扭角，结论按预印本处理。

实现边界：实验候选；需完整散射/纹理模型，当前无可复用跃迁表。

## 外围 moiré 方法基准

### 12. Relativistic Mott transition and high-order van Hove singularity in twisted double bilayer WSe2: Mean-field and functional renormalization group study

Bilal Hawashin 等，2026。Physical Review B 113, 075107 (2026)；既有日报实现。[原始来源](https://doi.org/10.1103/tchf-w8h7) · [arXiv 2509.09398](https://arxiv.org/abs/2509.09398)。

ABBA WSe₂ 的 Γ 谷四层连续模型可检验截断收敛、鞍点与拟合方法。保留为外围可复用基准；单粒子部分自旋简并，不能计成交错磁模型。

实现边界：复用既有日报非相互作用模型；不声称实现 Mott 相图或 fRG。 公式入口：arXiv v2 Eqs.1-5; honeycomb TB fit Eq.6。

## 中间转角尝试所用的关键文献

### 13. Nonrelativistic Spin-Momentum Coupling in Antiferromagnetic Twisted Bilayers

Ran He、Dan Wang、Nannan Luo、Jiang Zeng、Ke-Qiu Chen、Li-Ming Tang，PRL **130**, 046401 (2023)。[期刊](https://doi.org/10.1103/PhysRevLett.130.046401) · [出版社全文](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.130.046401/fulltext)。

此前 NiCl₂ 转角项目的材料参照。其无 SOC 双层自旋劈裂提供对称性背景；本轮通用单轨道参数仍是机制演示，不能当作该文的材料能带拟合。本次核验出版社摘要与全文开篇，列为此前研究已用文献。

### 14. Super-moiré spin textures in twisted two-dimensional antiferromagnets

King Cho Wong、Ruoming Peng、Eric Anderson 等，Nature Nanotechnology **21**, 359–365 (2026)，2026-02-02 在线发表。[期刊原文](https://www.nature.com/articles/s41565-025-02103-y)。

此前用户提供文章及后续磁序控制的参照。CrI₃ 双双层的磁纹理尺度可超过结构莫尔周期，其趋势由交换、各向异性与 DMI 等竞争决定。用于区分几何周期和磁响应长度；不把固定共线 NiCl₂ 电子谱或单层螺旋控制称为该磁纹理的复现。

## 本轮取舍

1. **Cairo 六位点原生预设**：复用 11/07 已核验的 S5–S7，重点验证 Hermitian 性、自旋分块、g 波节点以及应变后的 d 波分量。它是微观单层输入，可与几何工具组合，但不把组合自动称为某篇扭转材料的复现。
2. **公度几何工具**：在旋转与堆叠文献指导下，复用此前三角晶格八组已验证整数超胞、共同晶格和原胞位点映射；检查两层公度闭合、唯一位点数及倒格正交关系。具体整数构造属于软件实现，不声称照搬 Liu/Zeng 的全部模型。
3. **保留既有模型与外围基准**：Guo 已实现，维持解耦近似和来源版本标签。WSe₂ 适合后续截断收敛与鞍点测试；其正常态模型不是交错磁。MnSe、Shi、CrO 进入后续候选，先补齐完整系数、规范或自洽求解，再宣称复现。

以上实施优先级是根据现有工具范围作出的判断。本文没有把元数据核验等同于全文模型复现；各条阅读深度保存在 JSON。

## 作者元数据修正

旧仓库文档存在以下需要修正的署名；这些修正不改变已有数值模型。

| 旧文档简称 | 核验后 | 一手证据 |
| --- | --- | --- |
| “Jin 等，PRB 110, 174410” | **Sike Zeng、Yu-Jun Zhao** | [APS 文章作者栏](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.174410) |
| “Järvinen 等，Nano Lett. 26, 5078” | **Qirui Cui 等** | [PMC 原文作者栏](https://pmc.ncbi.nlm.nih.gov/articles/PMC13107460/) |
| “González-Hernández 等，Nano Lett. 6c01391” 用作首作者简称 | **Alberto M. Ruiz 等**；合作者为 Andrei Shumilin、Rafael González-Hernández、José J. Baldoví | [PMC 原文作者栏](https://pmc.ncbi.nlm.nih.gov/articles/PMC13352942/) |
| Guo 谷极化论文统一写作三人 | **期刊版四人**：San-Dong Guo、Yichen Liu、Junxi Yu、Cheng-Cheng Liu；**arXiv v1 三人**，作为公式来源应保留版本差异 | [期刊作者栏](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.L220402) · [v1](https://arxiv.org/abs/2406.13950v1) |

## 检索与获取限制

以 altermagnetism、twisted bilayer、moiré、tight-binding、continuum、spin-layer 等组合搜索，并回查 DOI/arXiv。只采纳原始论文、出版社和作者提交的 arXiv 信息；二手搜索命中仅用于定位。

MnSe 的完整电子模型未取得；Mu 的 HTML 全文不可用，本次仅核验摘要和作者。Sheng 的期刊页面获取失败，发表 DOI 由 arXiv 的 related DOI 核实，作者顺序按 arXiv v1 记录。CrO 与 Mu 均未发现可核实的期刊归属，按预印本处理。对这些条目未作逐图复现承诺。

