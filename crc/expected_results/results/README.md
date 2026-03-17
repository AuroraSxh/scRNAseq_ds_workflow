# Results Guide

`results/` 是本项目所有正式分析产物的总出口。

## 目录索引

```text
results/
├── README.md
├── qc/                           预处理和质控结果
├── integration/                  跨数据集整合结果
├── annotation/                   major cell type 注释结果
├── cd8/                          CD8 T 细胞精细分群结果
├── macrophage/                   myeloid/macrophage 精细分群结果
├── communication/                CellChat 与 macrophage-CD8 通讯结果
└── bio-foundation-housekeeping/  项目初始化与结构化报告
```

## 推荐查看顺序

1. `qc/`
   先确认两个数据集的过滤前后质量
2. `integration/`
   看整合后的全局 UMAP
3. `annotation/`
   看 major cell type 注释是否合理
4. `cd8/`
   看 CD8 子群划分、marker 和 TYW 相关图
5. `macrophage/`
   看 macrophage / myeloid 子群
6. `communication/`
   看 macrophage 与 CD8 的通讯总结

## 各目录一句话说明

- `qc/`
  - 质控统计表和 QC 图
- `integration/`
  - 整合对象的聚类、marker 和全局嵌入图
- `annotation/`
  - 大类细胞注释及 marker 支撑图
- `cd8/`
  - CD8 细分、合并、marker 和 TYW 相关分析
- `macrophage/`
  - macrophage/myeloid 子群注释和可视化
- `communication/`
  - CellChat 全量输出以及 macrophage-CD8 定向整理结果
- `bio-foundation-housekeeping/`
  - 项目结构和初始化报告

## 最值得优先打开的文件

- `integration/figures/integrated_umap.png`
- `annotation/figures/major_annotation_umap.png`
- `cd8/figures/cd8_merged_umap.png`
- `macrophage/figures/macrophage_umap.png`
- `communication/mac_t/macrophage_to_cd8_heatmap.png`
- `communication/mac_t/cellchat_overview.png`
