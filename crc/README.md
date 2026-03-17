# CRC Single-Cell Analysis Workspace

本目录是一个结直肠癌单细胞转录组分析工作区，核心任务包括：

- 两个公开数据集的预处理与整合
- 大类细胞注释
- CD8 T 细胞精细分群
- macrophage/myeloid 精细分群
- macrophage 与 CD8 间的 CellChat 通讯分析

这份 README 的目标不是解释算法细节，而是让你快速知道：

- 每个文件夹放什么
- 主流程从哪一步开始
- 最重要的结果在哪里看
- 新增文件以后应该往哪里放

## 1. 推荐阅读顺序

如果你只是想快速看结果，建议按这个顺序：

1. `results/integration/figures/`
   看整合后的全局 UMAP
2. `results/annotation/figures/`
   看 major cell type 注释
3. `results/cd8/`
   看 CD8 精细分群、marker 和 TYW 基因表达
4. `results/macrophage/`
   看 myeloid/macrophage 子群结果
5. `results/communication/mac_t/`
   看 macrophage 与 CD8 的互作汇总、overview 和 network 图

## 2. 顶层目录总览

```text
CRC/
├── README.md                         项目总指南
├── config/                           参数与 marker panel 配置
├── data/                             标准化数据入口
│   ├── raw/                          原始数据软链接入口
│   ├── interim/                      QC 后的中间 h5ad
│   └── processed/                    后续分析直接使用的 h5ad
├── GSE146771/                        原始数据目录
├── GSE178341/                        原始数据目录
├── logs/                             各步骤运行日志
├── reference/                        参考图片、错误记录、外部示例
│   ├── examples/                     样式参考图
│   └── errors/                       临时错误记录
├── results/                          全部分析输出
├── schemas/                          项目结构/housekeeping 相关说明
├── src/                              主分析脚本与公共工具
├── workflows/                        一键运行脚本
├── archive/                          不参与主流程的归档文件
│   ├── vendor_zips/                  手动下载的源码压缩包
│   └── scratch/                      临时文件/自动生成的杂项输出
├── requirements.txt                  Python 依赖
├── .venv/                            Python 虚拟环境
└── .r_libs/                          本地 R 包库
```

## 3. 关键目录说明

### `config/`

- `params.yaml`
  - 主流程参数，如 UMAP、聚类、CD8/myeloid 细分参数
- `marker_panels.yaml`
  - major cell type、CD8 状态等 marker 基因集合

### `data/`

标准化的数据层，建议以后都优先从这里理解项目状态。

- `data/raw/`
  - 指向原始数据目录的软链接
  - 当前已有：
    - `data/raw/GSE146771 -> ../../GSE146771`
    - `data/raw/GSE178341 -> ../../GSE178341`
- `data/interim/`
  - 预处理和 QC 后的中间对象
  - 当前关键文件：
    - `gse146771_qc.h5ad`
    - `gse178341_qc.h5ad`
- `data/processed/`
  - 后续所有分析直接读取的对象
  - 当前关键文件：
    - `crc_integrated_clusters.h5ad`
    - `crc_integrated_annotated.h5ad`
    - `cd8_subclusters.h5ad`
    - `cd8_merged.h5ad`
    - `macrophage_subclusters.h5ad`

说明：

- 当前 `src/01_prepare_qc.py` 和 `src/05_refresh_qc_figures.py` 仍直接读取顶层 `GSE146771/` 和 `GSE178341/`。
- `data/raw/` 是为了让目录语义更清晰，便于后续逐步统一到标准入口。

### `results/`

这是最重要的目录，按分析主题拆分。

#### `results/qc/`

- `tables/`
  - QC 汇总表
- `figures/`
  - QC 图

#### `results/integration/`

- 整合、聚类和 marker 相关输出
- `figures/integrated_umap.png`
  - 全局整合 UMAP，总入口图之一

#### `results/annotation/`

- major cell type 注释结果
- `figures/major_annotation_umap.png`
  - 大类群注释总图
- `subgroup_refinement_summary.md`
  - CD8 / macrophage 细分与合并标准说明

#### `results/cd8/`

- CD8 精细分群主目录
- 关键文件：
  - `cd8_cluster_annotation.csv`
  - `cd8_cluster_markers.csv`
  - `cd8_subgroup_counts.csv`
  - `figures/cd8_umap.png`
  - `figures/cd8_state_guided_umap.png`
  - `figures/cd8_marker_bubbleplot.png`
  - `figures/cd8_merged_umap.png`
- `tyw_genes/`
  - TYW 家族相关分析结果

#### `results/macrophage/`

- macrophage / myeloid 精细分群主目录
- 关键文件：
  - `myeloid_cluster_annotation.csv`
  - `figures/macrophage_umap.png`

#### `results/communication/`

通讯分析主目录，分成两个层次：

- `cellchat/`
  - CellChat 全量对象和全局输出
  - 关键文件：
    - `cellchat.rds`
    - `cellchat_checkpoint.rds`
    - `cellchat.png`
    - `cellchat_outgoing.png`
    - `cellchat_all_interactions.csv`
    - `cellchat_myeloid_cd8_interactions.csv`
- `mac_t/`
  - 面向“macrophage vs CD8”这个具体问题整理后的结果展示目录
  - 关键文件：
    - `macrophage_to_cd8_heatmap.png`
    - `cd8_to_macrophage_heatmap.png`
    - `macrophage_to_cd8_top_pathways.png`
    - `cd8_to_macrophage_top_pathways.png`
    - `cellchat_overview.png`
    - `cellchat_outgoing_overview.png`
    - `*_circle.png`
    - `*_hierarchy.png`
    - `README.md`

### `src/`

按步骤编号组织的主脚本目录。

- `01_prepare_qc.py`
  - 读取两个原始数据集，做基础 QC，输出 `data/interim/*.h5ad`
- `02_integrate_cluster.py`
  - 整合、降维、聚类
- `03_annotate_major.py`
  - major cell type 注释
- `04_cd8_subtype.py`
  - CD8 子群识别与 marker 提取
- `05_refresh_qc_figures.py`
  - 只刷新 QC 图
- `06_refresh_umaps.py`
  - 刷新整合图和注释图
- `07_state_umap_and_bubbleplots.py`
  - CD8 状态 UMAP 和 bubble plot
- `08_cd8_tyw_expression.py`
  - TYW 基因表达相关图
- `09_refine_subgroups_and_interactions.py`
  - CD8/macrophage 合并与互作准备
- `10_export_cellchat_inputs.py`
  - 导出 CellChat 输入矩阵和 metadata
- `11_run_cellchat.R`
  - 运行 CellChat、生成 incoming/outgoing overview 和 network 图
- `12_summarize_mac_t_cellchat.py`
  - 将 CellChat 结果整理成 `mac_t/` 下的热图、表格和 pathway 总结图
- `crc_sc_integration/`
  - 公共工具函数

### `logs/`

按步骤拆分日志，适合排错时直接定位。

- `logs/prepare/`
- `logs/integrate/`
- `logs/annotation/`
- `logs/cd8/`

### `reference/`

用于存放：

- 外部参考图片
- 希望模仿的可视化样式
- 临时报错文本

当前例如：

- `cellchat.png`
- `chat2.png`
- `chat3.png`
- `error.txt`

建议以后继续把“参考图”放这里，不要混入 `results/`。

第二轮整理后，`reference/` 进一步拆成：

- `reference/examples/`
  - 放样式参考图、外部截图、希望模仿的版式
- `reference/errors/`
  - 放临时报错文本或问题记录

### `archive/`

这是本次整理后新增的归档区，用来收纳不参与当前主流程、但也不想删除的文件。

- `archive/vendor_zips/`
  - 手动下载的压缩包
  - 当前包括：
    - `CellChat-main.zip`
    - `presto-master.zip`
- `archive/scratch/`
  - 临时或自动生成但不属于正式结果的文件
  - 当前包括：
    - `Rplots.pdf`

## 4. 当前主流程怎么跑

### Python 环境

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 分步运行

```bash
source .venv/bin/activate
python src/01_prepare_qc.py
python src/02_integrate_cluster.py
python src/03_annotate_major.py
python src/04_cd8_subtype.py
python src/07_state_umap_and_bubbleplots.py
python src/09_refine_subgroups_and_interactions.py
python src/10_export_cellchat_inputs.py
Rscript src/11_run_cellchat.R
python src/12_summarize_mac_t_cellchat.py
```

### 一键运行

```bash
bash workflows/run_crc_pipeline.sh
```

## 5. 最值得优先查看的结果

如果你只想看最终结论，优先打开这些文件：

- `results/integration/figures/integrated_umap.png`
- `results/annotation/figures/major_annotation_umap.png`
- `results/cd8/figures/cd8_merged_umap.png`
- `results/macrophage/figures/macrophage_umap.png`
- `results/communication/mac_t/macrophage_to_cd8_heatmap.png`
- `results/communication/mac_t/cd8_to_macrophage_heatmap.png`
- `results/communication/mac_t/cellchat_overview.png`
- `results/communication/mac_t/cellchat_outgoing_overview.png`
- `results/communication/mac_t/macrophage_to_cd8_mhc_i_circle.png`
- `results/communication/mac_t/cd8_to_macrophage_ccl_circle.png`

## 6. 新文件建议放哪里

以后新增文件时，建议遵循下面的约定：

- 新原始数据：放顶层原始数据目录，或先放外部盘后通过 `data/raw/` 建软链接
- 新中间对象：放 `data/interim/`
- 新正式分析对象：放 `data/processed/`
- 新结果图/表：放对应 `results/<topic>/`
- 新参考图片：放 `reference/`
- 新下载压缩包：放 `archive/vendor_zips/`
- 临时导出的杂项文件：放 `archive/scratch/`

不建议再把以下类型文件直接放在项目根目录：

- 压缩包
- 临时 pdf
- 与正式流程无关的截图
- 某一步的单独测试输出

## 7. 当前目录清晰度结论

这次整理后，目录可以按三层理解：

1. 输入层：`GSE*/`、`data/raw/`
2. 处理层：`src/`、`config/`、`workflows/`、`data/interim/`、`data/processed/`
3. 输出层：`results/`、`logs/`、`reference/`、`archive/`

以后如果继续扩展项目，优先保持这个结构，不要再把新的结果或下载文件散落到根目录。

## 8. 第二轮整理补充

本轮额外做了三件事：

- `reference/` 拆分为 `examples/` 和 `errors/`
- 新增 [results/communication/README.md](/Users/aurorasxh/codex_test/scrnaseq/CRC/results/communication/README.md) 作为通讯分析总导航
- `src/01_prepare_qc.py` 与 `src/05_refresh_qc_figures.py` 统一改为通过 `data/raw/` 读取原始数据

这意味着之后如果要替换原始数据路径，优先维护 `data/raw/` 的软链接即可，不需要再改顶层数据目录引用。
