# TSPO 卵巢癌多组学 —— 代码说明

对应当前版本论文初稿（456 行）与 9 张主图 + 补充材料。

## 目录结构（按论文图序）

```
code/
  Fig1_表达与网络/          2 个   TSPO 表达格局、共表达网络、功能模块
  Fig2_bulk共表达与免疫/      5 个   免疫浸润、WGCNA、高/低表达差异、NK 与 CTL
  Fig3_来源分解/            4 个   髓系与上皮来源的方差分解、跨队列复现
  Fig4_上皮代谢/            6 个   上皮混合模型、虚拟敲除、亚群构成、细胞类型表达谱
  Fig5_髓系配体与CAF/       4 个   髓系配体、CAF 活化、CPTAC 蛋白、化疗配对
  Fig6_髓系轨迹/            3 个   伪时序、化疗时点、细胞类型水平
  Fig7_空间转录组/          2 个   切片内相关、三区划分、区域差异表达
  Fig8_空间蛋白组/          2 个   CODEX 区室结构、MIF 区域归属
  Fig9_上游与药物/          7 个   甲基化、转录因子、oncoPredict、生存分析
  补充图/                  3 个   补充图 S1–S4、未分类 spot 再分型
  作图样式/                11 个   统一样式与九张图的绘制
  总控/                    4 个   环境、流程、论文数字汇总
  _archive/                6 个   论文中已删除的分析（分子对接、GDSC 验证、早期合并版）
```

## 运行顺序

前置：`00_setup.R`（定义全部路径）→ `01_build_matrix.R` → `99_session_info.R`

主分析按图序：

```
02 → 17 → 23 → 18 → 19 → 16        图 1、图 2
35 → 31 → 29 → 30                  图 3
03 → 04 → 20 → 21 → 41 → 44        图 4
05 → 24 → 10 → 42                  图 5
06 → 07 → 26                       图 6
08 → 22                            图 7
09 → 28                            图 8
25 → 13 → 14 → 11 → 12 → 20 → 43   图 9
```

作图：`作图样式/fig_theme.R` → `fig01.R … fig09.R` → `补充图/figS.R`

汇总：`27_manuscript_numbers.R` 把全部数字汇成一张表供正文核对

## 本轮新增的四个脚本

| 脚本 | 内容 | 输出 |
|:--|:--|:--|
| `41_celltype_expression.R` | 12 类细胞的配体与受体谱；CD74/CXCR4 专项；TSPO 与炎症小体基因的逐病人相关 | T43a–T43d |
| `42_chemo_paired_programs.R` | 化疗初治 vs 间歇减灭术配对比较（2 区室 × 11 程序）；受体变化；TSPO~SPP1 耦合 | T44, T44b, T44c |
| `43_program_survival.R` | 程序评分与 OS/PFI 的 Cox，含年龄阳性对照；两来源四分组 | T45, T45b, T45c |
| `44_vko_table.R` | 虚拟敲除的模块 Z 值、表达匹配对照、两靶点差值 | T46 |

## 五个容易踩的坑

### 1. Seurat v5 的分层对象

`GSE266577_raw.rds` 与 `X2_myeloid_ready.rds` 均为 v5 分层对象。
取数据前须 `JoinLayers()`，且 `GetAssayData()` 会报错，须改用
`SeuratObject::LayerData(obj, assay = "RNA", layer = "data")`。

### 2. 配对比较的 z 标准化位置

`42_chemo_paired_programs.R` 必须在**整类细胞上** z 标准化一次再取子集均值。
若逐子集 z 标准化，均值恒为 0，配对差将全部为零——这是最容易出错的一步。

### 3. 生存分析的阳性对照

任何生存分析都必须先证明该队列能测出已知的预后因素。
`43_program_survival.R` 以年龄为阳性对照：OS 终点 HR=1.33（p=2.5×10⁻⁴），
终点有效，因此 TSPO 与全部程序评分的阴性可判为真实阴性；
PFI 终点上年龄 p=0.52，该终点判别力不足，不作结论。

### 4. 生存文件的实际终点

`v2_6.8_survival_data.csv` 的 `time`／`event` 实为 **PFI**，
数值与标准 Xena 表的 PFI 完全一致（TSPO 1.02，p=0.802）。
分析 OS 须用 `D:/Rstudy/tcga/TCGA_OV_survival.tsv`。

### 5. 样本条码格式

表达矩阵列为 `TCGA-61-1910-01`（连字符），
oncoPredict 的 IC50 矩阵行名为 `TCGA.20.1683`（点号）。
匹配前须分别 `gsub("-","")` 与 `gsub("\\.","")` 再取前 12 位。

## 数据路径

```
D:/Rstudy/tcga/results_v2/EXPR_clean.rds          TCGA-OV 表达，18,747 × 304
D:/Rstudy/tcga/TCGA_OV_survival.tsv               标准 Xena 生存表（OS/PFI/DFI）
D:/Rstudy/tcga/OV_clinicalMatrix.tsv              年龄与病理
D:/Rstudy/tcga2/GSE266577_processed/...            137,083 细胞，22 类，32 例
D:/Rstudy/tcga/TSPO_spatial_processed.rds          空间转录组，16,626 × 4,386
D:/Rstudy/tcga/TSPO_trajectory/X2_myeloid_ready.rds  髓系 9,098 细胞，5 供体
D:/Rstudy/tcga/results_v2/sc_vko_*.csv             虚拟敲除输出
```

## 运行环境

```
R 4.5.3
Seurat 5.5.0      v5 分层对象
limma, lme4, WGCNA, slingshot, decoupleR, dorothea, GSVA
clusterProfiler, oncoPredict 1.2, pROC, CellChat
ggplot2 4.0.3, patchwork 1.3.2, circlize 0.4.18
```

## 合并版

`TSPO_全部代码_合并版.R` 是 53 个在用脚本按图序拼接成的单文件（773 KB，16,404 行），
便于检索与通读；实际运行请用分类文件夹中的独立脚本。
## 行数构成

| 部分 | 脚本 | 行数 |
|:---|:---:|:---:|
| 分析脚本（总控 + Fig1–9 + 补充图） | 39 | 9,636 |
| 作图脚本（fig_theme + fig01–09 + 样式） | 11 | 6,768 |
| **在用合计** | **50** | **16,404** |
| 归档（已删除的分析 + 历史拼接版） | 6 | 22,808 |

作图脚本单文件偏长（fig04_epi.R 1,654 行、fig01.R 746 行），因为每张主图都包含
十来个面板的独立绘制、坐标轴与图例微调；分析脚本平均约 250 行。

`_archive/历史合并文件/` 下是四个早期拼接版（共 22,158 行），
内容与编号脚本重复，不参与合并，仅作留档。