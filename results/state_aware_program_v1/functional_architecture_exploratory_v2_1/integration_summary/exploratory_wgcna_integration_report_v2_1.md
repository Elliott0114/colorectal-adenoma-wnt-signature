# 腺瘤上皮共表达模块与扰动整合分析报告（v2.1）

## 结论

WGCNA值得纳入稿件，但其最合适的角色是**对1,843-gene programme进行探索性功能分解**，而不是建立新的signature或提供独立外部验证。

增量发现是：模块的扰动响应将programme分成两个不同步的部分。M02、M03和M06构成腺瘤上调响应组；M04、M05和M10构成高度相关的成熟功能下降响应组，M09与后者方向相近但一致性较弱。ASCL2敲除强烈逆转前三个上调模块，却没有同步恢复成熟功能下降模块。这为“APC–WNT/ASCL2干预只能部分逆转上皮身份重塑”提供了模块层证据。

## 主要结果

### 1. 网络结构与内部重复

- 8,221个可检验基因在ABS、GOB和TAC三个状态的供者独立残差表达中形成11个非灰色consensus modules。
- 7个模块在held-out common effect中达到方向一致的富集标准：M02、M03、M04、M05、M06、M09和M10。
- 这7个模块在三个held-out状态中均获得Zsummary ≥ 2；各模块最低Zsummary为3.14–45.06。
- M03、M04、M05、M06、M09和M10与1,843-gene programme或重复通路community显著重叠；M02为更宽泛的腺瘤上调模块。
- 模块边界对consensus quantile有敏感性。仅M02在q=0和q=0.5两种审计中均达到Jaccard ≥ 0.50；其他模块应解释为探索性功能单元，而不是唯一的固定基因集合。

### 2. 五队列外部投影

GSE8671最初误用了仅含282行的历史compact-signature映射。v2.1使用完整GPL570注释后，各模块覆盖率为90.7%–96.5%，五个外部队列均可评分。

- 没有模块达到“REML–Knapp–Hartung 95% CI完全位于预期方向且所有leave-one-cohort-out汇总同向”的严格外部门槛。
- M02：汇总方向性效应0.324，95% CI −0.283至0.931；五次leave-one-cohort-out汇总均为正。
- M06：汇总方向性效应0.172，95% CI −0.377至0.721；五次leave-one-cohort-out汇总均为正。
- M09呈稳定的反向趋势，汇总效应−0.915，95% CI −2.033至0.204，I²=94.3%。
- 51对FFPE组织中，M02和M06分别有82.4%和72.5%的患者出现预期方向变化（Wilcoxon P=1.20×10⁻⁶和3.95×10⁻⁴）。

因此，模块不能作为可迁移的低维检测signature。M02和M06是相对一致的探索性腺瘤上调模块，而成熟功能下降模块的效应更依赖队列和组织背景。

### 3. 正交数据

- M02和M06在Becker snRNA-seq中同向（效应1.171和1.021；P=2.67×10⁻⁷和1.29×10⁻⁴），并在六张空间切片中全部同向。
- M06在CRC Atlas中同向（效应0.289，P=0.0168）。
- M04、M05、M09和M10在CRC Atlas中呈现预期的成熟功能下降，但在Becker或空间数据中部分反向，提示癌与腺瘤、单细胞与组织混合背景之间存在情境差异。
- 16个测量型sentinel获得至少一种蛋白组或空间支持；全部标记为`exploratory_only`。没有regulatory node通过全部预设门槛。

### 4. 模块级扰动响应

九个遗传或药理扰动比较形成两个主要响应组：

- M02–M03、M02–M06和M03–M06的扰动谱相关系数分别为0.958、0.909和0.950。
- M02/M03/M06与M04/M05/M10的相关系数为−0.823至−0.942。
- ASCL2 knockout对M02、M03和M06的方向性逆转分别为0.909、1.244和1.423；对M04、M05、M09和M10分别为−0.752、−1.271、−1.220和−1.183。

这说明ASCL2适合作为**检验腺瘤上调调控臂是否可逆的实验探针**，但不能被表述为整个1,843-gene programme的总开关。后续类器官实验应同时测量：

- 上调调控臂：ASCL2以及M03/M06读出，如EPHB2、ZNRF3和MYO9B；
- 成熟功能臂：M04代表基因CA2、FABP1或HMGCS2。

若ASCL2抑制只降低前者而未恢复后者，该结果本身就是对“部分逆转”模型的实验验证。

## 稿件分流建议

### 正文

加入一段简洁的exploratory consensus co-expression结果，重点写：

1. programme被分解为方向一致且held-out保存的模块；
2. 成熟功能下降模块连接到已重复的代谢/线粒体通路；
3. 扰动谱显示上调调控臂与成熟功能臂并不同步逆转。

不在摘要中增加WGCNA，不把模块称为biomarker、portable signature或therapeutic target。

### Figure 4候选

- a：跨分区通路重复热图；
- b：代表性running-enrichment曲线；
- c：11个模块在discovery与held-out八个排序情境中的方向性富集；
- d：模块与重复功能community的连接；
- e：模块×扰动响应热图。

### 补充材料候选

- consensus gene dendrogram与模块色带；
- 五个外部队列、paired FFPE、Becker、RNA–ATAC、CRC Atlas和空间数据的完整方向矩阵；
- 参数敏感性、所有模块成员、kME、外部效应和16个sentinel的完整表。

## 文件

- 主图候选：`figures/communications_biology_v3.0/functional_architecture_wgcna_candidate_v2_1/figure4_functional_architecture_wgcna_candidate_v2_1.pdf`
- 补充图候选：`figures/communications_biology_v3.0/functional_architecture_wgcna_candidate_v2_1/figureS_wgcna_structure_and_context_candidate_v2_1.pdf`
- 整合结果：`module_integrated_evidence.tsv`
- 扰动谱相关：`module_perturbation_profile_correlations.tsv`
- 完整模块验证：`../module_validation.tsv`
- 实验蛋白候选：`../protein_priorities.tsv`

## Material Passport

- **输入来源：** 冻结的8,221-gene ranking、1,843-gene programme、discovery/held-out donor-aware pseudobulk、五个公共腺瘤队列、51对公共FFPE组织、公共RNA–ATAC、CRC Atlas、空间转录组、蛋白组及公开扰动数据。
- **数据级别：** 公共去标识化表达数据和项目内派生表；本分析不处理直接身份信息。
- **网络建模：** 供者为独立单位；每个状态先去除normal/adenoma主效应，再对供者残差建signed consensus WGCNA；power=15，consensus quantile=0.25，q=0和0.5仅作敏感性审计。
- **校正记录：** v2的GSE8671映射不完整；v2.1在看到纠正后结果前冻结完整GPL570输入修正规则，原v2结果保留为审计。
- **输出可追溯性：** 输入、软件版本与输出SHA-256记录于`exploratory_wgcna_integration_manifest.json`及各阶段manifest。
- **证据状态：** post-result exploratory functional decomposition；不改变已冻结programme及主要统计结论。
