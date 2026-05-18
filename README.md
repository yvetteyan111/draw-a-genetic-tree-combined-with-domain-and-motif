下面是可以直接复制使用的 README：

````markdown
# Phylogenetic Tree, Domain and Motif Visualization Script

## 1. 功能简介

本脚本用于为每个 pathway step 自动生成整合型进化树 SVG 图。图中同时展示：

- 蛋白序列系统发育树
- Known / Predicted 序列标签
- 物种来源颜色
- Pfam domain 结构域分布
- MEME motif 基序分布
- 可点击跳转链接，包括蛋白详情页、Pfam 页面和 MEME 结果页面

脚本适用于每个 step 文件夹已经完成 MAFFT、trimAl、IQ-TREE、HMMER 和 MEME 分析后的结果可视化。

## 2. 输入文件结构

每个 step 文件夹下需要包含以下文件：

```text
stepX/
├── merged_aligned_trimmed.fasta.treefile
├── domain_hits.txt
├── merged.fasta
└── meme_results/
    └── meme.xml
````

其中：

| 文件                                      | 说明                         |
| --------------------------------------- | -------------------------- |
| `merged_aligned_trimmed.fasta.treefile` | IQ-TREE 生成的进化树文件           |
| `domain_hits.txt`                       | HMMER/Pfam 结构域比对结果         |
| `merged.fasta`                          | Known 和 Predicted 蛋白序列合并文件 |
| `meme_results/meme.xml`                 | MEME motif 识别结果            |

## 3. 输出结果

脚本会在每个 step 文件夹下自动创建 `draw_tree/` 子文件夹，并输出：

```text
stepX/
└── draw_tree/
    ├── stepX.svg
    ├── parsed_domain_table1.csv
    └── parsed_motif_table1.csv
```

输出文件说明：

| 文件                         | 说明            |
| -------------------------- | ------------- |
| `stepX.svg`                | 最终整合型进化树图     |
| `parsed_domain_table1.csv` | 解析后的 domain 表 |
| `parsed_motif_table1.csv`  | 解析后的 motif 表  |

## 4. 特殊 step 命名规则

如果 step 文件夹名为：

```text
step8.19
```

脚本会自动识别其中的两个 step 编号：

```text
step8
step19
```

并生成两个 SVG 文件：

```text
step8.19/
└── draw_tree/
    ├── step8.svg
    └── step19.svg
```

这两个 SVG 使用同一套树、domain 和 motif 数据，但跳转链接中的 `step=` 参数不同。

例如：

```text
/protein/known/by-label/?step=step8&label=...
/protein/known/by-label/?step=step19&label=...
```

该功能由 `expand_step_codes()` 和 `for (link_step_code in step_codes_for_links)` 循环实现。

## 5. 主要函数说明

### `get_step_files(step_dir)`

根据 step 文件夹路径自动生成输入文件和输出文件路径。

### `check_required_files(files)`

检查必须文件是否存在，包括：

* treefile
* domain_hits.txt
* meme.xml

如果缺失，脚本会停止运行。

### `expand_step_codes(step_name)`

从 step 文件夹名中提取所有数字。

例如：

```r
expand_step_codes("step8.19")
```

返回：

```r
c("step8", "step19")
```

### `read_domain_table(domain_file)`

读取 HMMER domtblout 格式的 domain 结果，并筛选：

```r
domain_score >= 20
i_Evalue <= 1e-5
```

### `resolve_domain_overlap(df)`

处理同一条序列上重叠的 domain，只保留 i-Evalue 更好的 domain。

### `read_meme_xml(meme_file)`

解析 MEME XML 文件，提取每条序列上的 motif 起止位置。

### `draw_one_step(step_dir)`

核心绘图函数。输入一个 step 文件夹路径，生成对应 SVG 图和解析后的 CSV 文件。

## 6. 使用方法

在 R 中运行：

```r
draw_one_step("/home/ziyan/enzyem_pipeline/results/visualization/step8.19")
```

如果需要批量处理多个 step 文件夹：

```r
base_dir <- "/home/ziyan/enzyem_pipeline/results/visualization"

step_dirs <- list.dirs(base_dir, recursive = FALSE)

results <- lapply(step_dirs, function(step_dir) {
  tryCatch(
    draw_one_step(step_dir),
    error = function(e) {
      message("处理失败：", step_dir)
      message("错误信息：", e$message)
      NULL
    }
  )
})

results_df <- dplyr::bind_rows(results)
write.csv(results_df, file.path(base_dir, "tree_visualization_summary.csv"), row.names = FALSE)
```

## 7. 可调整参数

脚本开头可以调整以下参数：

| 参数                         | 作用                         |
| -------------------------- | -------------------------- |
| `SHOW_TIP_LABEL`           | 是否显示叶节点标签                  |
| `SHOW_NODE_SUPPORT`        | 是否显示节点支持率                  |
| `TREE_LINEWIDTH`           | 进化树线条粗细                    |
| `USE_TOPOLOGY_ONLY`        | 是否忽略 branch length，只显示拓扑结构 |
| `KNOWN_COLOR`              | Known 序列颜色                 |
| `MAX_SPECIES_IN_LEGEND`    | 图例中最多显示的物种数量               |
| `MEME_POSITION_ZERO_BASED` | MEME motif 位置是否从 0 开始      |

## 8. SVG 交互功能

生成的 SVG 支持点击跳转：

| 点击对象           | 跳转目标                        |
| -------------- | --------------------------- |
| Known 序列标签     | Known protein detail 页面     |
| Predicted 序列标签 | Predicted protein detail 页面 |
| Domain 方块      | InterPro Pfam 页面            |
| Motif 方块       | 本地 MEME HTML 页面             |

## 9. 注意事项

1. treefile 中的 tip label 必须和 `domain_hits.txt`、`meme.xml`、`merged.fasta` 中的序列名一致。
2. Known 序列标签需要以 `Known|` 或 `*Known|` 开头。
3. Predicted 序列标签一般应包含 `IMP_ID|Species|PREDICTED|Confidence`。
4. 如果某个 step 名中包含多个数字，例如 `step8.19`，脚本会自动生成多个 SVG 文件。
5. 输出 SVG 中的链接依赖 Django 后端对应路由已经正确配置。

## 10. 适用场景

该脚本适合用于展示每个反应 step 的同源序列结果。它可以将序列进化关系、domain 保守性和 motif 保守性整合在同一张图中，便于用户判断预测酶序列是否与已知酶具有一致的进化和结构特征。

```
```
