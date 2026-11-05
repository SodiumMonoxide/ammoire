# 交错磁电子紧束缚模型文献目录

检索截至 **2026-11-05**。优先收录同行评审原始论文中的显式电子哈密顿量、跃迁表、Wannier 构造，以及提供电子跃迁骨架的相互作用模型。按电子模型层次分组，同一模型在不同应用中的变体保留可查出处。用 “altermagnet*”、“tight-binding”、“Hubbard”、“Wannier”、“layer group” 与材料、晶格名交叉查找，并沿 2019 年的前史和核心论文引文回查。缺少 altermagnet 一词的旧文、正文外的补充材料和不可检索的数据仍可能漏收。

表中“可构造”指原文给出了足以搭建相应模型的表达式或数据入口。“材料基准”指可比对能带、波函数、节点或响应。某个模型是否随本工具提供可执行预设，以软件模型列表为准。

## 1. 对称性模型族和公式入口

| 原始论文 | 模型族及复现入口 |
| --- | --- |
| Hayami、Yanagi、Kusunose，JPSJ **88**, 123702（2019） | 术语确立前的共线模型。四亚晶格烧绿石和四方晶格，分别用团簇多极矩刻画各向异性跃迁及偶宇称自旋劈裂。文中给出具体晶格示例和 32 个晶体点群的对应关系。[论文](https://doi.org/10.7566/JPSJ.88.123702) · [全文](https://arxiv.org/abs/1908.08680) |
| Hayami、Yanagi、Kusunose，PRB **102**, 144441（2020） | 从团簇和键的增强多极矩系统构造实空间电子哈密顿量。可对照共线偶宇称劈裂与其他磁序的不同模型机制。[论文全文](https://doi.org/10.1103/PhysRevB.102.144441) |
| Šmejkal、González-Hernández、Jungwirth、Sinova，Sci. Adv. **6**, eaaz8809（2020） | 金红石棋盘格的极简实空间构造和 RuO₂ 的 Wannier 模型。磁性位点、配体位置及 SOC 决定晶体 Hall 响应，适合校验 Néel 方向和 SOC 的实现。[论文及补充材料](https://doi.org/10.1126/sciadv.aaz8809) |
| Šmejkal、Sinova、Jungwirth，PRX **12**, 031042（2022） | 非相对论自旋群分类，把交错磁劈裂划分为平面型和体型，最低阶可为 d、g、i 波。KRu₄O₈ 的晶场与电子模型是微观校验例。[论文](https://doi.org/10.1103/PhysRevX.12.031042) |
| Roig 等，PRB **110**, 144412（2024） | 40 个中心对称空间群与二重 Wyckoff 位点的电子模型，跨单斜、正交、四方、三方、六方和立方晶系。表 III 收录顺磁极限有 BZ 边界节点面的非对称型群，表 IV 收录相应无节点面模型。两表列出亚晶格跃迁和允许的 SOC 项，表 V–VI 给材料参数。[论文](https://doi.org/10.1103/PhysRevB.110.144412) · [含公式及参数表的全文](https://par.nsf.gov/servlets/purl/10637117) |
| Gondolf 等，PRB **111**, 174436（2025） | 限定于原始晶格且有反演对称二重 Wyckoff 位点的二维层群。附录 A 给一个含同亚晶格各向异性跃迁的 dxy 模型和单带近似，附录 B 表 II 给整个层群模型族。表 I 给复现参数。[论文](https://doi.org/10.1103/PhysRevB.111.174436) · [含二维层群表的全文](https://par.nsf.gov/servlets/purl/10637130) |
| Ezawa，PRB **111**, 125420（2025） | 两自旋低能模型与二维、三维周期格点化，覆盖 d、g、i 波交错磁，也单列奇宇称 p、f 波磁体。方格二维 g 波见式 (101)，体六方 g 波见式 (108)，六方二维 i 波小动量目标见式 (109) 起，立方三维 i 波见式 (124)。式 (116) 的二维动能中含印刷上的 k_z 项，二维实现需使用二维晶格动能并声明此处约定。[论文](https://doi.org/10.1103/PhysRevB.111.125420) · [公式全文](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.111.125420) |
| Ezawa，PRB **113**, L241301（2026） | 非线性自旋热输运的晶格化谐波变体，含体 g 波和二维 i 波。式中节点与坐标方向可作为检验。f 波属于另行讨论的奇宇称磁体。[论文](https://doi.org/10.1103/fnwg-lfzg) |

### Gondolf 二维模型对应表

用 cᵢ = cos kᵢ、sᵢ = sin kᵢ 表示。表中“各向异性”是乘 τ_z 的因子，每行另有允许的跃迁系数。论文使用含胞内坐标的 Bloch 规范表示简化谐波，另在式 (A2) 给出周期规范。以下合并了有相同因子的条目。[原文表 II](https://par.nsf.gov/servlets/purl/10637130#page=8)

| 层群 / 位点 | 各向异性 | 最低阶 |
| --- | --- | --- |
| L15 / 2a–2b，L16 / 2a–2b，L17 / 2a–2b，L40 / 2a–2b，L44 / 2a–2b，L63 / 2b | sₓsᵧ | dxy |
| L51 / 2c | a(cₓ−cᵧ)+bsₓsᵧ | d 的两分量 |
| L61 / 2c | cₓ−cᵧ | dx²−y² |
| L63 / 2a | sₓsᵧ(cₓ−cᵧ) | g |
| L7、L38、L41、L42 的表列位点 | 0 | 在该最小基底中无此自旋劈裂 |

附录 A 的 dxy 原型在周期规范下有 ε₀ = −2t₂(cₓ+cᵧ)−4t₃cₓcᵧ−μ，H_AB = −t₁(1+e⁻ⁱᵏˣ)(1+e⁻ⁱᵏʸ)，t_z = 4t₄sₓsᵧ，以及交错场 Nτ_zσ_z。参数 (t₁,t₂,t₃,t₄,μ,N) = (0.425,0.05,−0.025,−0.075,0.3,0.2) eV。单带和双亚晶格版本采用不同的单位胞和 BZ。[附录 A 与表 I](https://par.nsf.gov/servlets/purl/10637130#page=8)

## 2. 微观晶格、轨道与关联电子模型

| 原始论文 | 晶格、电子自由度及复现入口 |
| --- | --- |
| Brekke、Brataas、Sudbø，PRB **108**, 224421（2023） | 二维 Lieb 型最小微观模型，保留磁性与非磁性位点、电子跃迁及局域交换。电子和磁振子呈 d 波劈裂，超导配对来自磁振子介导相互作用。可先复现正常态电子块。[论文](https://doi.org/10.1103/PhysRevB.108.224421) · [全文](https://arxiv.org/abs/2308.08606) |
| Das、Leeb、Knolle、Knap，PRL **132**, 263402（2024） | 光学方格的 t–t′–δ Hubbard 模型。两磁性亚晶格之间由最近邻连接，同亚晶格的次近邻各向异性交替。Hartree–Fock 有金属和绝缘 d 波区，适合按固定交错场和自洽场分别对照。[论文](https://doi.org/10.1103/PhysRevLett.132.263402) |
| Ferrari、Valentí，PRB **110**, 205140（2024） | Shastry–Sutherland 四位点 Hubbard 模型，正交二聚键配合交错序形成 d 波电子劈裂。变分蒙特卡洛谱函数提供单粒子与关联结果对照。[论文](https://doi.org/10.1103/PhysRevB.110.205140) · [全文](https://arxiv.org/abs/2408.00841) |
| Bose、Vadnais、Paramekanti，PRB **110**, 205120（2024） | 方形–八边形晶格多轨道 t–J 模型。半满附近有强耦合交错序，另有准一维 van Hove 奇点促进的弱耦合 d 波区。投影约束与配对需要单独求解。[论文](https://doi.org/10.1103/PhysRevB.110.205120) |
| Giuli、Mejuto-Zaera、Capone，PRB **111**, L020401（2025） | 相互作用驱动的多位点方格电子模型。可对照由巡游磁性产生的 AM 与预设交错交换的电子块。[论文](https://doi.org/10.1103/PhysRevB.111.L020401) · [全文](https://arxiv.org/abs/2410.00909) |
| Antonenko、Fernandes、Venderbos，PRL **134**, 096703（2025） | 二维 Lieb dx²−y² 与三维金红石 dxy。补充材料式 (S8)–(S10) 给二维无 SOC 两磁位点模型，式 (S33)–(S34) 补回配体，式 (S59)–(S63) 给三维跃迁。保留 SOC 时可检验二维镜面 Chern 能带和三维镜面保护 Weyl 节线。[论文](https://doi.org/10.1103/PhysRevLett.134.096703) · [含公式的补充材料](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.134.096703/supplement.pdf) |
| Zhu、Huo 等，PRL **134**, 166701（2025） | 从自旋团簇搭建 d/g 等不同对称的晶格模型。补充材料给 1/5 缺位方格四位点 Hubbard 平均场和多种团簇构造。[论文](https://doi.org/10.1103/PhysRevLett.134.166701) · [含跃迁及团簇的补充材料](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.134.166701/altermagnet_suppl-v2.pdf) |
| Engineering Altermagnetic States，PRL **135**, 036701（2025） | 筛查 34 个方形铺砌，给出 Lieb、fes（方形–八边形）和 tts（3.3.4.3.4）三种候选网格。适合按结构选取不同磁性位点和界面跳跃。[论文](https://doi.org/10.1103/v38b-5by1) |
| Kaushal、Franz，PRL **135**, 156502（2025） | 修饰 Lieb Hubbard，反 CuO₂ 布局。每胞 2、4 电子的 Mott 交错序和掺杂的 dx²−y² 金属，通过无约束 Hartree–Fock 和精确对角化检验。[论文](https://doi.org/10.1103/s31h-hk2v) · [全文](https://arxiv.org/abs/2412.16421) |
| Inverse Design of 2D Altermagnetic Metal–Organic Framework Monolayers，JACS Au **5**, 381–387（2025） | 用 Hückel 非键分子轨道逆向构造二维金属有机框架的交错序，为明确配体选择的模型提供参照。[论文](https://doi.org/10.1021/jacsau.4c01150) |
| Altermagnetism and strain induced altermagnetic transition，npj QM（2025） | 二维 Cairo 五边形，显式保留磁性和非磁性原子。晶格畸变从 g 波引入 d 波分量，可用来检验应变允许的低阶项。[论文](https://doi.org/10.1038/s41535-025-00793-0) · [全文](https://arxiv.org/abs/2412.16857) |
| Dürrnagel 等，PRL **135**, 036502（2025） | Lieb 金属中由亚晶格干涉驱动的巡游交错磁相变，提供与轨道有序机制不同的微观实现。[论文](https://doi.org/10.1103/2g3v-z76q) · [全文](https://arxiv.org/abs/2412.14251) |
| Rao、Knolle、Classen，PRB（2025） | 两种方格 d 波模型的自旋极化 van Hove 点及竞争不稳定性。可分开检查 dxy 和 dx²−y² 的奇点位置、填充和弱耦合配对假设。[论文](https://doi.org/10.1103/zfk4-rxh1) · [全文](https://arxiv.org/abs/2505.02786) |
| He 等，PRB（2025） | 重访 t–t′–δ 方格 Fermi–Hubbard 半填充相图。同一个电子跃迁骨架的关联校验，避免将每一篇相图计算重复计为新单粒子模型。[论文](https://doi.org/10.1103/4mv8-tb66) |
| Huo 等，PRB **113**, 024420（2026） | Lieb 晶格自旋团簇实现 d/g 波多种磁胞，和 SOC 联合产生高阶拓扑。条带或开边界才可对照边缘、角模结果。[论文](https://doi.org/10.1103/9wcm-pmr2) |
| Zhang 等，PRB **113**, 035411（2026） | Cairo 五边形金属有机框架 TM₂(TCNQ)₂ 的 p–d 模型，显式分析 [110] 应变的 g 到 d 波及等离激元。[论文](https://doi.org/10.1103/1zbl-4q6j) |
| d’Ornellas 等，PRB **113**, 024426（2026） | 非晶格点每位点含两个方向性轨道，自旋旋转不变的相互作用自发产生自旋–轨道对称性破缺。适合实空间检验，需要单独定义动量谱函数。[论文](https://doi.org/10.1103/zqnk-153f) |
| Blanco Peces、Merino，PRB **113**, 165143（2026） | Kagome Dirac 填充的 Hubbard 模型可由库仑相互作用形成交错序，兼讨论 Lieb–kagome。保留自发扩大磁胞的自由度对复现关键。[论文](https://doi.org/10.1103/x7xc-dkk9) |
| Dong 等，PRB **113**, 245117（2026） | 方格单轨道扩展 Hubbard t–U–V。由相互作用自发产生 d 波，不预设方向性晶场或双轨道。[论文](https://doi.org/10.1103/k5vw-c9ks) |
| Yu、Lado，Adv. Sci. e78167（2026） | 石墨烯 antidot 超晶格，单层及 Bernal 双层。原子 p_z 电子紧束缚 + Hubbard 平均场验证平面 i 波，晶格镜面对换自旋亚晶格。这是六角形平面 i 波材料设计的直接模型参照。[论文](https://doi.org/10.1002/advs.78167) |
| Petermann 等，PRR **8**, 033359（2026） | 双曲规则格双分格紧束缚，次近邻跳跃产生交错劈裂。Abelian 双曲准动量至少四维，计算采用其双曲平移结构。[论文](https://doi.org/10.1103/h2ys-c3rz) |

## 3. 材料拟合、配体和实验基准

| 原始论文 | 复现入口与适用范围 |
| --- | --- |
| Naka 等，Nat. Commun. **10**, 4305（2019） | κ 型有机晶体棋盘状二聚体的分子轨道 Hubbard 模型，Methods 式 (1) 给转移积分的组织方式。相邻层的反相堆叠可使自旋流相消，适合作为双层设定的校验。[论文](https://doi.org/10.1038/s41467-019-12229-y) |
| Ma 等，Nat. Commun. **12**, 2846（2021） | V₂Se₂O 型 C 对称配对自旋–谷锁定、应变与自旋流。原始单层材料与谷基准，可用于交叉核对六成员拟合。[论文](https://doi.org/10.1038/s41467-021-23127-7) |
| Guo 等，Mater. Today Phys.（2023） | MAGNDATA 多材料 DFT 调研，另有作者公开配套数据。可用 CoF₂、FeSO₄F、LiFe₂F₆、CrNb₄S₈、RuO₂、CrSb 等筛选基准。导入单个体系前需核实数据究竟是能带还是 Wannier 跃迁。[论文](https://doi.org/10.1016/j.mtphys.2023.100991) · [数据](https://zenodo.org/records/7448505) |
| Bhowal、Spaldin，PRX **14**, 011019（2024） | MnF₂ 的多极矩、晶场和非相对论 d 波模型。配体与局域高阶多极矩的微观解释可校验金红石玩具模型。[论文](https://doi.org/10.1103/PhysRevX.14.011019) |
| Attias、Levchenko、Khodas，PRB **110**, 094425（2024） | 掺杂 FeSb₂ 的材料紧束缚和有限 SOC Hall 响应。按自旋群（无 SOC）及磁群（有限 SOC）分别限制跃迁。[论文](https://doi.org/10.1103/PhysRevB.110.094425) · [全文](https://par.nsf.gov/servlets/purl/10590908) |
| Roig 等，PRB **110**, 144412（2024） | RuO₂ 与 MnF₂ 的金红石原型式 (5)–(8)，FeSb₂ 式 (31)，κ-Cl 式 (F1)–(F2)，表 V 的单轨道参数及表 VI 的 RuO₂ 双轨道参数。CrSb / MnTe 用 SG194 的 B₁g g 波式 (32)–(34)，立方 SG223 的 A₂g i 波见式 (35)–(36)。CrSb / MnTe 的单轨道式用于对称性，低能轨道简并需额外保留。[公式及参数表](https://par.nsf.gov/servlets/purl/10637117) |
| Yang 等，Nat. Commun. **16**, 1442（2025） | CrSb 的 3D ARPES 与基于 Wannier 的紧束缚。Sb 5p 介导 Cr 3d 的第三近邻通道控制大劈裂，k_z 的水平节点与三张垂直节点面可用于验证。数据来自作者仓库。[论文](https://doi.org/10.1038/s41467-025-56647-7) · [源数据](https://doi.org/10.6084/m9.figshare.25865749) |
| Yu、Suh、Roig、Agterberg，Nat. Commun. **16**, 2950（2025） | κ-Cl 的普通态紧束缚与重合 van Hove patch 模型。九种二维非对称型空间群可提供相应的简并结构，适合检查填充与弱耦合不稳定性。[论文](https://doi.org/10.1038/s41467-025-57970-9) |
| Li 等，Commun. Phys. **8**, 311（2025） | 对称化 Cr d 与 Sb p Wannier 哈密顿量，拟合窗口约 −2 到 2 eV。用来计算体 Weyl 点的手性和表面费米弧。[论文 Methods](https://doi.org/10.1038/s42005-025-02232-9) |
| Autieri 等，PRB **112**, 014412（2025） | 八面体转动下 Sr₂RuO₄ 的轨道选择性 t₂g 模型。理论单层为 g 波，体层间跃迁允许 dxy 分量。可对照层间作用怎样降低单层谐波的约束。[论文](https://doi.org/10.1103/ssxp-gz9l) · [全文](https://arxiv.org/abs/2501.14378) |
| Naka 等，JPSJ **94**, 083702（2025） | κ-Cl 第一性原理约束下的应变相关跃迁与有效 Hubbard / Heisenberg 模型。剪切可使不等价位点产生 s 波磁化分量，适合区分应变导致的均匀与交错响应。[论文](https://doi.org/10.7566/JPSJ.94.083702) · [作者公开稿](https://arxiv.org/abs/2505.07327) |
| Microscopic origin，Curr. Appl. Phys. **79**, 29–33（2025） | MnF₂ 对照的简洁双亚晶格模型，原子交换、配体笼畸变和各向异性跃迁共同给出自旋劈裂。[论文](https://doi.org/10.1016/j.cap.2025.07.004) |
| Wang 等，Adv. Funct. Mater.（在线 2025，卷期 2026） | 筛出 MnS₂、CoS₂、MnC₂、CoPSe 的五边形单层。论文用 Wannier90 构建 MnS₂ 纳米盘模型以分析角态，可作为 Cairo 模型的独立材料对照。[论文](https://doi.org/10.1002/adfm.202505145) |
| Switchable Altermagnetism Induced by Polyhedral Rotation Distortion，Nano Lett.（2026） | MnX₂（X=O、S、Se、Te）配位多面体转动与应变使共线 AFM / AM 切换。补充材料式 (S1) 起给每自旋的两亚晶格紧束缚及分壳层跃迁。[论文及补充材料](https://doi.org/10.1021/acs.nanolett.6c01573) |
| Altermagnetism, kagome flat band, and Weyl fermion states，npj QM（2026） | 插层 XY₄Z₈ 族的 DFT、SOC 下 Weyl 态和 Methods 中的 Y 网络紧束缚。网络最小模型解释平带来源，完整自旋劈裂比较采用全材料电子结构。[论文](https://doi.org/10.1038/s41535-026-00885-5) |
| Lange 等，PRX **16**, 031066（2026） | 堆叠 Lieb 电子模型与表面自旋群；体反铁磁的边界可产生二维交错劈裂。NaMnP 和 FeGe₂ 分别提供 d 与 g 波表面参照。[论文](https://doi.org/10.1103/6j6l-kb6x) · [全文](https://arxiv.org/abs/2602.08773) |
| Long 等，Nature（2026） | CrSb 量子振荡与三维阶参量基准。实球谐 yz(3x²−y²) 可校验节点和符号。此项是材料验证数据。[论文](https://doi.org/10.1038/s41586-026-10902-z) |
| Yang 等，Adv. Sci. e77945（2026） | 双层 MnPTe₃ 的光电流与层间混成。每自旋基底含四个有效磁性轨道与六个有效配体轨道，适合从显式磁位点–配体隧穿构造双层校验。[论文](https://doi.org/10.1002/advs.77945) |

材料标签需要区分样品和边界。RuO₂ 与 FeSb₂ 的模型可作为条件性结构原型，具体样品的基态应与实验分别核对。KV₂Se₂O 的表面 ARPES d 波与体反铁磁证据已有表面自旋群解释，移植到双层时要明确所用的是体切片、表面还是孤立单层。[CrSb 论文关于 RuO₂ 的讨论](https://doi.org/10.1038/s41467-025-56647-7) · [FeSb₂ 磁矩实验](https://doi.org/10.1103/PhysRevB.111.075141) · [表面自旋群与 KV₂Se₂O](https://doi.org/10.1103/6j6l-kb6x)

## 4. SOC、拓扑、光学和配对变体

| 原始论文 | 可复现模型要点 |
| --- | --- |
| Rao、Mook、Knolle，PRB **110**, 024425（2024） | Γ 附近二维 d 波 + Rashba + 外场的低能参照，用来比较 Chern 相和光学电导对磁场方向的依赖。[论文](https://doi.org/10.1103/PhysRevB.110.024425) |
| Field-sensitive dislocation bound states，PRB **110**, 165141（2024） | 方格 d 波交错项 + Rashba + 外交换，另构造 s 波超导 BdG。用于比较缺陷的弱拓扑指数及外场方向。[论文](https://doi.org/10.1103/PhysRevB.110.165141) |
| Sorn，PRB **111**, L161109（2025） | 金红石电子格点加交错场。以边界切割得到修饰 SSH 链族，检查反手性表面态与畴壁态的终止面依赖。[论文](https://doi.org/10.1103/PhysRevB.111.L161109) · [全文](https://arxiv.org/abs/2411.11949) |
| Parshukov、Wiedmann、Schnyder，PRB **111**, 224406（2025） | 四位点方格的自旋极化拓扑交叉，原文给电子 TB 与独立磁振子模型。可对照破缺对称性的质量项与电荷或自旋 Chern 数。[论文](https://doi.org/10.1103/PhysRevB.111.224406) · [全文](https://arxiv.org/abs/2403.09520) |
| Orbital-spin locking，PRB **112**, L020401（2025） | 显式轨道的 d 波微观 TB。用于动量相关自旋–轨道成分与光学跃迁，轨道基底和偶极矩应随能带一并传递。[论文](https://doi.org/10.1103/bzzy-ngcs) |
| Quasisymmetry-Constrained Spin Ferromagnetism，PRL **135**, 016703（2025） | 对照最小电子模型与 DFT 的 SOC 驱动弱净自旋磁化。SOC 开启后应分别检查自旋与轨道响应。[论文](https://doi.org/10.1103/839n-rckn) |
| Constraints on superconducting pairing，PRB **112**, 014516（2025） | 对比方格亚晶格晶场各向异性四带原型、低能投影和相互作用诱导 AM。配对选择规则对被投影掉的亚晶格权重敏感。[论文](https://doi.org/10.1103/zylh-rqxl) |
| Berry curvature-induced transport signature，npj QM（2025） | 方格最近邻和次近邻自旋相关跃迁，构造简并点附近的低能模型和 Berry 曲率分布。[论文](https://doi.org/10.1038/s41535-025-00805-z) |
| Topological altermagnetic Josephson junctions，npj QM（2026） | d 波交错磁结区加两侧配对，引入晶轴相对结方向。用正常四带模型生成 Nambu 块，研究拓扑相与方向控制。[论文](https://doi.org/10.1038/s41535-026-00874-8) |
| Liu、Sun，Chin. Phys. B（2026） | 方格两亚晶格、交错 d 波交换及亚晶格势，得到不同 Chern 数的相图。拓扑复现需保留其矩阵结构与边界条件。[论文](https://doi.org/10.1088/1674-1056/ae42bb) |
| Light-induced topological phase transitions，PRB（2026） | 方格四带模型含 p_z 与 d_xz±i d_yz 轨道、反转质量、d 波亚晶格各向异性，加入周期光驱动后分析 Floquet 拓扑与热输运。[论文](https://doi.org/10.1103/8x4v-gf3y) |
| Arucu，Sci. Rep. **16**, 32096（2026） | d 波 t₂g 多轨道紧束缚，使用完整投影的 L_x、L_y、L_z 与 Kubo 算符。可校验轨道 Hall 和自旋／轨道 splitter 响应。[论文](https://doi.org/10.1038/s41598-026-75295-5) |

## 5. 录用及预印本追踪

以下按查到的原始来源状态单独记录。录用页提供 DOI，但检索时未定位到相应的期刊正式出版页。

| 状态 | 模型 |
| --- | --- |
| APS 2026-09-24 录用 | Cheng 等，V₂Se₂O 族六种层状材料的对称性约束拟合，含 SOC、应变、谷及层自由度。适合后续对接材料化莫尔模型，复现时需固定所用预印本或期刊版本。[APS 录用页](https://journals.aps.org/prmaterials/accepted/10.1103/8czs-74vv) · [预印本](https://arxiv.org/abs/2602.09465) |
| 预印本 | Camerano 等，以 Co₁/₄NbSe₂ 的 Wannier 轨道剔除试验分解配体通道。[arXiv](https://arxiv.org/abs/2605.21284) |
| 预印本 | AgCrP₂S₆ 链嵌入二维，外场打破层外对称性，低能 TB 把劈裂关联到第三近邻链间各向异性跳跃。[arXiv](https://arxiv.org/abs/2607.16856) |
| 预印本 | 金属自由三角烯蜂窝框架的 d 波模型，以配位方向控制 π 轨道最近邻重叠。[arXiv](https://arxiv.org/abs/2604.17386) |
| 预印本 | Chen 等，V₂STeO 与 VO 族的自旋、谷、动量锁定和门控交错势，含螺旋／手性输运模型。[arXiv](https://arxiv.org/abs/2603.06487) |
| 预印本 | 二维方格双轨道工程比较 d、g 波的起源。[arXiv](https://arxiv.org/abs/2605.25128) |

## 6. 本工具内置预设与原文的对应

| 预设 | 对应来源和模型边界 |
| --- | --- |
| `square-dxy`、`square-dx2-y2`、`square-g`、`hex-i`、`hex-g-3d`、`cubic-i` | 两分量周期谐波模型。方格 dx²−y² 采用代表性格点谐波，六角平面 i 波采用二维动能。三维 g/i 预设可取固定 k_z 切片。[Ezawa 格点化](https://doi.org/10.1103/PhysRevB.111.125420) |
| `lieb-dx2-y2` | 两磁性位点的四带版本，取主文使用的局域交错交换，默认关闭 SOC，可启用式 (S39) 的 SOC。此版本省略配体显式轨道、J′ 磁多极项和净磁矩。[式 (S8)–(S10)、(S39)](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.134.096703/supplement.pdf) |
| `brekke-lieb` | 两磁性位点及一个非磁性位点的六带正常态电子模型。使用周期规范和 Fig. 2 的能量参数，单位胞边长取文中最近邻距离的两倍。[式 (21)](https://arxiv.org/abs/2308.08606) |
| `gondolf-dxy` | 二维四带 dxy benchmark，采用附录 A 的跃迁和表 I 参数。[附录 A 与表 I](https://par.nsf.gov/servlets/purl/10637130#page=8) |
| `gondolf-l15-2ab`、`gondolf-l16-2ab`、`gondolf-l17-2ab`、`gondolf-l40-2ab`、`gondolf-l44-2ab`、`gondolf-l51-2c`、`gondolf-l61-2c`、`gondolf-l63-2a`、`gondolf-l63-2b` | 表 II 九个可有非零各向异性的二维条目，采用一般跃迁系数。表中 τ_z=0 的条目可用禁用各向异性后的简并极限对照。[附录 B 与表 II](https://par.nsf.gov/servlets/purl/10637130#page=8) |
| `roig-ruo2`、`roig-mnf2`、`roig-fesb2`、`roig-kappa-cl` | 表 V 的单轨道材料原型，保留两磁位点和自旋。RuO₂ 双轨道及其余表 III–IV 构型在上方保留复现入口。[式 (5)–(6)、(30)–(31)、附录 F 和表 V](https://par.nsf.gov/servlets/purl/10637117) |
| `roig-194`、`roig-223` | SG194 六方体 g 与 SG223 立方体 i 的四带模型。参数用于研究对称性及其切片，材料定量拟合还需要相应轨道与跃迁。[式 (32)–(36)](https://par.nsf.gov/servlets/purl/10637117) |

## 7. 统一模型的自由度和验证量

- 单轨道投影常写成 H(k) = ε(k)σ₀+Δf(k) n·σ，适合扫描扭角、谐波方向和反向 Néel 轴。二维方格可用 cos kₓ−cos kᵧ 和 sin kₓsin kᵧ 作为两个 d 分量，乘积给一个 g 波。体六方的典型 g 波在 Γ 附近含 k_z 与三次平面谐波的乘积，因此在 k_z=0 时消失。平面六角晶格的六阶谐波和立方三维六阶谐波是两种不同的 i 波例子。[Ezawa 格点化](https://doi.org/10.1103/PhysRevB.111.125420) · [Roig 六方与立方公式](https://par.nsf.gov/servlets/purl/10637117)
- 双亚晶格 H(k) = ετ₀ + Re(h_AB)τₓ − Im(h_AB)τᵧ + t_zτ_z + Nτ_z n·σ 可区分晶格各向异性和局域交换。无 SOC、同一 k 时，E_{s,±} = ε ± √(|h_AB|²+(t_z+sN)²)。t_z=0 或 N=0 的极限给自旋简并。该能谱还能交叉校验实空间与动量空间矩阵的傅里叶符号。[Roig 式 (1)](https://par.nsf.gov/servlets/purl/10637117) · [Gondolf 附录 A](https://par.nsf.gov/servlets/purl/10637130)
- 相同的 Δf(k) 可能来自不同的磁位点、配体和轨道。界面隧穿、亚晶格选择性及轨道矩阵元需要完整 Bloch 向量。若换用含胞内坐标的 Bloch 规范，层间矩阵也要同步变换。固定 Bloch 规范后，检查 H(k+G) 的规范等价性、晶格旋转联同自旋翻转时的谱关系、节点，以及关掉交换和各向异性时的简并。[Lieb 配体投影](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.134.096703/supplement.pdf) · [轨道–自旋锁定](https://doi.org/10.1103/bzzy-ngcs)
- 有 SOC 或层间磁轴不共线时，逐能带同时记录三个自旋分量及层和轨道投影。在严格简并子空间内优先使用投影算符或迹平均。相互作用模型的磁序、配对和占据依赖于求解器，固定交错场下的电子能带只能检验该假设下的单粒子块。[镜面 Chern 能带和 Weyl 节线](https://doi.org/10.1103/PhysRevLett.134.096703) · [修饰 Lieb Hubbard](https://doi.org/10.1103/s31h-hk2v)
