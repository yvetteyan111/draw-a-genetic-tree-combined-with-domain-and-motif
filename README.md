
# 系统发育树-结构域-基序一体化 SVG 绘图脚本说明

本脚本用于批量绘制每个代谢通路步骤（step）的综合可视化图，包含以下内容：

- 系统发育树（Phylogenetic Tree）
- 蛋白质结构域（Domain Architecture）
- MEME 基序（Motif Architecture）
- 可点击的序列标签
- 可点击的 Pfam 结构域块
- 可点击的 Motif 块
- 自动输出解析后的 Domain 和 Motif 表格

---

# 1. 脚本文件

主绘图脚本：
/home/ziyan/enzyem_pipeline/scripts/draw_trees/draw_one_step.R

批量运行脚本：
/home/ziyan/enzyem_pipeline/scripts/draw_trees/batch_draw_tree_domian_motif.R

---

# 2. 输入文件要求

每个 step 文件夹中必须包含以下文件：
stepX/
├── merged_aligned_trimmed.fasta.treefile
├── domain_hits.txt
├── merged.fasta
└── meme_results/
    └── meme.xml

# 3. 输出文件

脚本会在每个 step 文件夹下自动创建：
draw_tree/

并输出以下文件：
stepX_tree_domain_motif.svg
parsed_domain_table.csv
parsed_motif_table.csv

# 4. 主函数

核心函数：
draw_one_step(step_dir)

使用示例：
source("/home/ziyan/enzyem_pipeline/draw_one_step.R")
draw_one_step("/home/ziyan/enzyem_pipeline/results/visualization/step1")

# 5. 批量运行

执行批量脚本：
source("/home/ziyan/enzyem_pipeline/scripts/batch_draw_tree_domian_motif.R")

该脚本会自动扫描目录：
/home/ziyan/enzyem_pipeline/results/visualization/

下的所有 step 文件夹，并逐个运行：
draw_one_step(step_dir)

最终生成汇总文件：
/home/ziyan/enzyem_pipeline/results/visualization/draw_tree_batch_summary.csv

# 6. get_step_files() 函数

该函数负责统一管理输入和输出文件路径：
get_step_files(step_dir)

# 7. Domain 外链功能

每个结构域块都可以点击跳转到 InterPro 的 Pfam 页面。

## 7.1 链接生成函数
make_domain_href()


## 7.2 处理逻辑

例如：PF08031.18,先去掉版本号，变成：PF08031,然后生成链接：
[InterPro Pfam Entry](https://www.ebi.ac.uk/interpro/entry/pfam/PF08031/?utm_source=chatgpt.com)

## 7.3 SVG 中的实现

每个结构域矩形通过：
wrap_svg_link(rect_body, r$href)
包裹为 SVG 超链接，用户点击后即可打开对应的 InterPro 页面。


# 8. 序列标签跳转

## Known 序列
/protein/known/by-label/?label=...

## Predicted 序列
点击跳转到：
/protein/predicted/by-label/?label=...

由以下函数控制：

make_tip_href(label, is_known)


# 9. Motif 跳转

Motif 块默认跳转到当前 step 的 MEME HTML 页面：
./meme_results/meme.html


由以下函数控制：
make_motif_href()

# 10. 自动布局功能

脚本会根据以下因素自动调整布局：

* 叶节点数量
* 标签长度
* 序列长度
* Domain 和 Motif 宽度
* 图例数量

自动计算：

* 树的缩放比例
* Domain 面板宽度
* Motif 面板宽度
* 字体大小
* 图整体宽度

由函数：
auto_layout_params()

实现。

---

# 11. 输出内容

SVG 图中包含：

* 系统发育树
* Bootstrap 支持值（可选）
* Known 与 Predicted 不同颜色显示
* Domain 结构图
* Motif 结构图
* 图例（Sequence type、Motif、Domain、Species）


# 12. 脚本主要功能概述

| 功能       | 说明                         |
| -------- | -------------------------- |
| 读取树文件    | 解析 IQ-TREE 输出的 `.treefile` |
| 解析结构域    | 读取 HMMER `domain_hits.txt` |
| 解析 Motif | 读取 MEME `meme.xml`         |
| 自动布局     | 根据标签长度和序列长度调整图形            |
| SVG 输出   | 生成高质量矢量图                   |
| 可点击链接    | 支持跳转到本地页面和 InterPro        |
| 批量处理     | 自动遍历所有 step                |

---

# 16. 作者 yvetteyan

