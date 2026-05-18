# =========================================================
# 0. 安装并加载包
# =========================================================
cran_pkgs <- c("ape", "ggplot2", "dplyr", "readr", "stringr", "xml2", "ggnewscale", "scales")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bio_pkgs <- c("ggtree")
for (p in bio_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(ape)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(xml2)
library(ggnewscale)
library(ggtree)
library(scales)

# =========================================================
# 1. 输入输出路径
# =========================================================
tree_file   <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/merged_aligned_trimmed.fasta.treefile"
domain_file <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/domain_hits.txt"
meme_file   <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/meme_results/meme.xml"
protein_fasta <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/merged.fasta"

folder_path <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/draw_tree"

if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print("文件夹创建成功！")
} else {
  print("文件夹已存在，无需创建。")
}


out_svg        <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/draw_tree/step6_ec.svg"
out_domain_csv <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/draw_tree/parsed_domain_table.csv"
out_motif_csv  <- "/home/ziyan/enzyem_pipeline/results/visualization_new/step6/draw_tree/parsed_motif_table.csv"

# =========================================================
# 2. 参数
# =========================================================
SHOW_TIP_LABEL    <- TRUE
SHOW_NODE_SUPPORT <- FALSE
SUPPORT_MODE      <- "full"   # "full" / "ufboot" / "alrt"
SUPPORT_MIN       <- NA

TREE_LINEWIDTH    <- 1.0
TREE_X_SCALE      <- 2.5
USE_TOPOLOGY_ONLY <- FALSE

TIP_SIZE_PRED  <- 3.0
TIP_SIZE_KNOWN <- 3.4
NODE_SIZE      <- 2.0

BACKBONE_LINEWIDTH <- 0.6
DOMAIN_LINEWIDTH   <- 3.0
MOTIF_LINEWIDTH    <- 3.2

MAIN_TITLE_SIZE <- 18
SUBTITLE_SIZE   <- 4

PDF_WIDTH  <- 18
PDF_HEIGHT <- 8

MEME_POSITION_ZERO_BASED <- TRUE
KEEP_ONLY_SPECIFIC_DOMAIN <- TRUE

# =========================================================
# 3. 工具函数
# =========================================================
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

clean_species <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x
}
#去掉accession号小数点后的内容
clean_pfam_accession <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\..*$", "", x)   # 去掉 .18 这种版本号
  x
}

# -------------------------
# 3.1 读取 FASTA 长度
# -------------------------
read_fasta_lengths <- function(fasta_file) {
  if (is.na(fasta_file) || !file.exists(fasta_file)) return(NULL)
  
  lines <- readLines(fasta_file, warn = FALSE)
  if (length(lines) == 0) return(NULL)
  
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) return(NULL)
  
  res <- vector("list", length(header_idx))
  
  for (i in seq_along(header_idx)) {
    s <- header_idx[i]
    e <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    
    header <- sub("^>", "", lines[s])
    header <- clean_label(header)
    
    seq_lines <- lines[(s + 1):e]
    seq_lines <- seq_lines[seq_lines != ""]
    seq <- paste(seq_lines, collapse = "")
    seq <- gsub("\\s+", "", seq)
    
    res[[i]] <- data.frame(
      label = header,
      seq_len = nchar(seq),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(res) %>%
    distinct(label, .keep_all = TRUE)
}

# -------------------------
# 3.2 读取 HMMER 原始 domtblout 文件
# -------------------------
read_domain_table <- function(domain_file, max_evalue = 1e-5, min_score = 20) {
  col_names <- c(
    "target_name", "target_accession", "tlen",
    "query_name", "query_accession", "qlen",
    "full_seq_E_value", "full_seq_score", "full_seq_bias",
    "domain_number", "domain_of", "c_Evalue", "i_Evalue",
    "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
    "description"
  )
  
  df <- readr::read_table(
    file = domain_file,
    comment = "#",
    col_names = col_names,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  out <- df %>%
    dplyr::mutate(
      label = clean_label(target_name),
      start = as.numeric(ali_from),
      end   = as.numeric(ali_to),
      domain = trimws(query_name),
      accession = clean_pfam_accession(query_accession),
      domain_score = as.numeric(domain_score),
      i_Evalue = as.numeric(i_Evalue)
    )
  
  out %>%
    dplyr::filter(
      !is.na(label), label != "",
      !is.na(start), !is.na(end),
      !is.na(domain), domain != "",
      !is.na(i_Evalue),
      !is.na(domain_score),
      domain_score >= min_score,
      i_Evalue <= max_evalue,
    ) %>%
    dplyr::select(label, start, end, domain, accession, domain_score, i_Evalue) %>%
    dplyr::distinct()
}

resolve_domain_overlap <- function(df) {
  
  df <- df %>% arrange(label, start, i_Evalue)
  
  result <- list()
  
  for (seq_id in unique(df$label)) {
    sub_df <- df[df$label == seq_id, ]
    
    keep <- rep(TRUE, nrow(sub_df))
    
    for (i in seq_len(nrow(sub_df))) {
      if (!keep[i]) next
      
      for (j in seq_len(nrow(sub_df))) {
        if (i == j || !keep[j]) next
        
        # 判断是否重叠
        overlap <- !(sub_df$end[i] < sub_df$start[j] ||
                       sub_df$end[j] < sub_df$start[i])
        
        if (overlap) {
          # 保留 i_Evalue 更小（更显著）的
          if (sub_df$i_Evalue[i] <= sub_df$i_Evalue[j]) {
            keep[j] <- FALSE
          } else {
            keep[i] <- FALSE
          }
        }
      }
    }
    
    result[[seq_id]] <- sub_df[keep, ]
  }
  
  bind_rows(result)
}

# -------------------------
# 3.3 读取 meme.xml
# -------------------------
read_meme_xml <- function(meme_file, zero_based = TRUE) {
  doc <- read_xml(meme_file)
  
  seq_nodes <- xml_find_all(doc, ".//sequence")
  seq_map <- data.frame(
    sequence_id = xml_attr(seq_nodes, "id"),
    label = xml_attr(seq_nodes, "name"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(label = clean_label(label)) %>%
    dplyr::filter(!is.na(sequence_id), sequence_id != "", !is.na(label), label != "") %>%
    dplyr::distinct(sequence_id, .keep_all = TRUE)
  
  motif_nodes <- xml_find_all(doc, ".//motif")
  motif_info <- data.frame(
    motif_id   = xml_attr(motif_nodes, "id"),
    motif_name = xml_attr(motif_nodes, "name"),
    motif_alt  = xml_attr(motif_nodes, "alt"),
    width      = suppressWarnings(as.numeric(xml_attr(motif_nodes, "width"))),
    stringsAsFactors = FALSE
  )
  
  if (nrow(motif_info) > 0) {
    motif_info$motif_display <- ifelse(
      !is.na(motif_info$motif_alt) & motif_info$motif_alt != "",
      motif_info$motif_alt,
      ifelse(
        !is.na(motif_info$motif_name) & motif_info$motif_name != "",
        motif_info$motif_name,
        motif_info$motif_id
      )
    )
  }
  
  contributing_sites <- xml_find_all(doc, ".//contributing_site")
  
  if (length(contributing_sites) > 0) {
    res_list <- vector("list", length(contributing_sites))
    
    for (i in seq_along(contributing_sites)) {
      node <- contributing_sites[[i]]
      
      seq_id  <- xml_attr(node, "sequence_id")
      pos_raw <- suppressWarnings(as.numeric(xml_attr(node, "position")))
      
      motif_node <- xml_find_first(node, "ancestor::motif[1]")
      motif_id   <- xml_attr(motif_node, "id")
      motif_name <- xml_attr(motif_node, "name")
      motif_alt  <- xml_attr(motif_node, "alt")
      motif_w    <- suppressWarnings(as.numeric(xml_attr(motif_node, "width")))
      
      motif_display <- ifelse(
        !is.na(motif_alt) && motif_alt != "",
        motif_alt,
        ifelse(!is.na(motif_name) && motif_name != "", motif_name, motif_id)
      )
      
      if (is.na(pos_raw) || is.na(motif_w)) next
      
      start <- if (zero_based) pos_raw + 1 else pos_raw
      end   <- start + motif_w - 1
      
      res_list[[i]] <- data.frame(
        sequence_id = seq_id,
        motif = motif_display,
        start = as.numeric(start),
        end   = as.numeric(end),
        stringsAsFactors = FALSE
      )
    }
    
    motif_df <- bind_rows(res_list) %>%
      dplyr::left_join(seq_map, by = "sequence_id") %>%
      dplyr::select(label, motif, start, end) %>%
      dplyr::filter(
        !is.na(label), label != "",
        !is.na(motif), motif != "",
        !is.na(start), !is.na(end)
      ) %>%
      dplyr::distinct()
    
    return(motif_df)
  }
  
  scanned_site_nodes <- xml_find_all(doc, ".//scanned_site")
  
  if (length(scanned_site_nodes) > 0) {
    res_list <- vector("list", length(scanned_site_nodes))
    
    for (i in seq_along(scanned_site_nodes)) {
      node <- scanned_site_nodes[[i]]
      
      motif_id <- xml_attr(node, "motif_id") %||% xml_attr(node, "motif")
      pos_raw  <- suppressWarnings(as.numeric(xml_attr(node, "position") %||% xml_attr(node, "pos")))
      
      parent_scanned_sites <- xml_find_first(node, "ancestor::scanned_sites[1]")
      seq_id <- xml_attr(parent_scanned_sites, "sequence_id")
      
      motif_row <- motif_info %>% dplyr::filter(motif_id == !!motif_id)
      
      motif_display <- motif_row$motif_display[1] %||% motif_id
      motif_w <- suppressWarnings(as.numeric(motif_row$width[1]))
      
      if (is.na(pos_raw) || is.na(motif_w)) next
      
      start <- if (zero_based) pos_raw + 1 else pos_raw
      end   <- start + motif_w - 1
      
      res_list[[i]] <- data.frame(
        sequence_id = seq_id,
        motif = motif_display,
        start = as.numeric(start),
        end   = as.numeric(end),
        stringsAsFactors = FALSE
      )
    }
    
    motif_df <- bind_rows(res_list) %>%
      dplyr::left_join(seq_map, by = "sequence_id") %>%
      dplyr::select(label, motif, start, end) %>%
      dplyr::filter(
        !is.na(label), label != "",
        !is.na(motif), motif != "",
        !is.na(start), !is.na(end)
      ) %>%
      dplyr::distinct()
    
    return(motif_df)
  }
  
  warning("在 meme.xml 中没有找到可识别的 motif 位点结构。")
  data.frame(label = character(), motif = character(), start = numeric(), end = numeric())
}

# -------------------------
# 3.4 bootstrap 标签
# -------------------------
get_support_labels <- function(node_labels_raw, mode = "full", min_value = NA) {
  mode <- match.arg(mode, c("full", "ufboot", "alrt"))
  if (is.null(node_labels_raw)) return(NULL)
  
  node_labels_raw[is.na(node_labels_raw)] <- ""
  
  if (mode == "full") {
    out <- node_labels_raw
  } else if (mode == "ufboot") {
    out <- ifelse(grepl("/", node_labels_raw), sub(".*/", "", node_labels_raw), node_labels_raw)
  } else {
    out <- ifelse(grepl("/", node_labels_raw), sub("/.*", "", node_labels_raw), node_labels_raw)
  }
  
  vals_num <- suppressWarnings(as.numeric(out))
  if (!is.na(min_value)) {
    out <- ifelse(!is.na(vals_num) & vals_num >= min_value, out, "")
  }
  
  out[is.na(out)] <- ""
  out
}

# -------------------------
# 3.5 SVG 工具
# -------------------------
svg_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

wrap_svg_link <- function(body, href = NA_character_) {
  if (is.na(href) || href == "") return(body)
  paste0('<a xlink:href="', svg_escape(href), '">', body, '</a>')
}

hex_col <- function(x) {
  rgb_mat <- grDevices::col2rgb(x)
  grDevices::rgb(rgb_mat[1, ] / 255, rgb_mat[2, ] / 255, rgb_mat[3, ] / 255)
}

build_edge_df <- function(tree_data) {
  child_df <- tree_data %>%
    dplyr::filter(parent != node) %>%
    dplyr::select(node, parent, x_child = x, y_child = y)
  
  parent_df <- tree_data %>%
    dplyr::select(parent_node = node, x_parent = x, y_parent = y)
  
  child_df %>%
    dplyr::left_join(parent_df, by = c("parent" = "parent_node"))
}

make_tip_href <- function(label, is_known, link_step_code = NULL) {
  step_part <- if (!is.null(link_step_code) && !is.na(link_step_code) && link_step_code != "") {
    paste0("step=", utils::URLencode(link_step_code, reserved = TRUE), "&")
  } else {
    ""
  }
  
  if (isTRUE(is_known)) {
    paste0(
      "/protein/known/by-label/?",
      step_part,
      "label=",
      utils::URLencode(label, reserved = TRUE)
    )
  } else {
    paste0(
      "/protein/predicted/by-label/?",
      step_part,
      "label=",
      utils::URLencode(label, reserved = TRUE)
    )
  }
}


make_domain_href <- function(accession) {
  accession <- clean_pfam_accession(accession)
  
  if (is.na(accession) || accession == "") return(NA_character_)
  
  paste0(
    "https://www.ebi.ac.uk/interpro/entry/pfam/",
    utils::URLencode(accession, reserved = TRUE),
    "/"
  )
}

make_motif_href <- function(label, motif) {
  "/static/pathways/berberine/all_step_phylogeny_svg/step6_meme.html"
}

# -------------------------
# 3.x 提取最后一个 EC 号
# -------------------------
extract_last_ec_block <- function(x) {
  x <- clean_label(x)
  # 使用 | 分割字符串
  parts <- unlist(strsplit(x, "\\|"))
  
  # 基础检查：确保至少有 3 个字段
  if (length(parts) < 3) return(NA_character_)
  
  # 去掉第一个字段前后的空格进行判断
  first_field <- trimws(parts[1])
  
  # 判断逻辑：是否以 "Known" 开头
  if (grepl("^Known", first_field)) {
    # 如果是以 Known 开头，取第 4 个字段（需确保存在第 4 个）
    if (length(parts) >= 4) {
      return(trimws(parts[4]))
    } else {
      return(NA_character_)
    }
  } else {
    # 否则取第 3 个字段
    return(trimws(parts[3]))
  }
}
# 用于显示成一列文本
extract_ec_display_text <- function(x) {
  ecs <- extract_last_ec_block(x)
  if (length(ecs) == 0) return(NA_character_)
  paste(ecs, collapse = "; ")
}
# -------------------------
# 3.x EC 调色
# -------------------------
make_ec_colors <- function(ec_levels) {
  ec_levels <- unique(ec_levels[!is.na(ec_levels) & ec_levels != ""])
  n <- length(ec_levels)
  if (n == 0) return(character())
  
  cols <- grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[1:n],
    c = 85,
    l = 62
  )
  setNames(cols, ec_levels)
}


# =========================================================
# 4. 读取树
# =========================================================
tree <- read.tree(tree_file)

if (is.null(tree$tip.label) || length(tree$tip.label) == 0) {
  stop("treefile 读取失败，没有检测到 tip labels。")
}

tree$tip.label <- clean_label(tree$tip.label)
tip_labels <- tree$tip.label

# =========================================================
# 5. 识别 Known / Predicted 并整理 label
# =========================================================
tip_info <- data.frame(
  label = tip_labels,
  is_known = grepl("^\\*?Known\\|", tip_labels),
  stringsAsFactors = FALSE
)

tip_info$species <- ifelse(
  tip_info$is_known,
  clean_species(sub("^\\*?Known\\|[^|]+\\|([^|]+).*$", "\\1", tip_info$label)),
  clean_species(sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", tip_info$label))
)

tip_info$seq_type <- ifelse(tip_info$is_known, "Known", "Predicted")

tip_info$known_id <- ifelse(
  tip_info$is_known,
  sub("^\\*?Known\\|([^|]+)\\|.*$", "\\1", tip_info$label),
  NA
)

# 新增：提取最后一个 EC
tip_info$ec_display <- vapply(tip_info$label, extract_ec_display_text, character(1))

tip_info$display_label <- tip_info$label

species_levels <- tip_info %>%
  dplyr::filter(!is_known) %>%
  dplyr::pull(species) %>%
  unique() %>%
  sort()

make_many_species_colors <- function(species_levels) {
  n <- length(species_levels)
  if (n == 0) return(c())
  
  # 物种较少时，直接用高区分度调色板
  if (n <= 12) {
    cols <- hcl.colors(n, palette = "Dark 3")
    return(setNames(cols, species_levels))
  }
  
  # 物种较多时，使用更大的 HCL 色轮，尽量拉开色相差异
  # 这里不用 hue_pal()，因为大数量时容易相邻颜色太接近
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  
  # 交替使用两组亮度和色度，让相邻颜色不要只靠色相区分
  l_vals <- rep(c(55, 70), length.out = n)
  c_vals <- rep(c(100, 65), length.out = n)
  
  cols <- grDevices::hcl(
    h = hues,
    c = c_vals,
    l = l_vals
  )
  
  setNames(cols, species_levels)
}

species_colors <- make_many_species_colors(species_levels)

KNOWN_COLOR <- "#C0392B"

tree$node.label <- get_support_labels(
  node_labels_raw = tree$node.label,
  mode = SUPPORT_MODE,
  min_value = SUPPORT_MIN
)

# =========================================================
# 6. 读取结构域和 motif
# =========================================================
domain_df <- read_domain_table(domain_file)
motif_df  <- read_meme_xml(meme_file, zero_based = MEME_POSITION_ZERO_BASED)

domain_df$label <- clean_label(domain_df$label)
motif_df$label  <- clean_label(motif_df$label)

cat("read_domain_table 后记录数:", nrow(domain_df), "\n")

domain_df <- domain_df %>% dplyr::filter(label %in% tree$tip.label)
motif_df  <- motif_df  %>% dplyr::filter(label %in% tree$tip.label)

cat("过滤 tree tip 后记录数:", nrow(domain_df), "\n")

domain_df <- resolve_domain_overlap(domain_df)

cat("resolve_domain_overlap 后记录数:", nrow(domain_df), "\n")

write.csv(domain_df, out_domain_csv, row.names = FALSE)
write.csv(motif_df, out_motif_csv, row.names = FALSE)

cat("结构域条目数：", nrow(domain_df), "\n")
cat("motif 条目数：", nrow(motif_df), "\n")
cat("domain 匹配到树 tip 的序列数：", length(unique(domain_df$label)), "\n")
cat("motif 匹配到树 tip 的序列数：", length(unique(motif_df$label)), "\n")
cat("Known 序列数：", sum(tip_info$is_known), "\n")
# =========================================================
# 7. 计算蛋白长度
# =========================================================
seq_len_df <- read_fasta_lengths(protein_fasta)

if (is.null(seq_len_df)) {
  seq_len_df <- bind_rows(
    domain_df %>% dplyr::transmute(label, pos = end),
    motif_df  %>% dplyr::transmute(label, pos = end)
  ) %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(seq_len = max(pos, na.rm = TRUE), .groups = "drop")
}

if (is.null(seq_len_df) || nrow(seq_len_df) == 0) {
  seq_len_df <- data.frame(
    label = tree$tip.label,
    seq_len = 100,
    stringsAsFactors = FALSE
  )
} else {
  max_len_all <- max(seq_len_df$seq_len, na.rm = TRUE)
  seq_len_df <- dplyr::full_join(
    data.frame(label = tree$tip.label, stringsAsFactors = FALSE),
    seq_len_df,
    by = "label"
  )
  seq_len_df$seq_len[is.na(seq_len_df$seq_len)] <- max_len_all
}

# =========================================================
# 8. 颜色
# =========================================================
domain_levels <- sort(unique(domain_df$domain))
motif_levels <- unique(motif_df$motif)

motif_num <- suppressWarnings(as.numeric(stringr::str_extract(motif_levels, "\\d+")))

motif_levels <- motif_levels[order(motif_num, motif_levels, na.last = TRUE)]

domain_colors <- if (length(domain_levels) > 0) {
  setNames(hcl.colors(length(domain_levels), "Set 2"), domain_levels)
} else {
  NULL
}

motif_colors <- if (length(motif_levels) > 0) {
  setNames(hcl.colors(length(motif_levels), "Dark 2"), motif_levels)
} else {
  NULL
}

species_colors_hex <- if (length(species_colors) > 0) {
  setNames(vapply(species_colors, hex_col, character(1)), names(species_colors))
} else {
  character()
}

if (length(species_colors_hex) > 0) {
  species_colors_hex <- species_colors_hex[order(names(species_colors_hex))]
}

domain_colors_hex <- if (!is.null(domain_colors) && length(domain_colors) > 0) {
  setNames(vapply(domain_colors, hex_col, character(1)), names(domain_colors))
} else {
  character()
}

motif_colors_hex <- if (!is.null(motif_colors) && length(motif_colors) > 0) {
  setNames(vapply(motif_colors, hex_col, character(1)), names(motif_colors))
} else {
  character()
}

known_color_hex <- hex_col(KNOWN_COLOR)
known_legend_hex <- c("Known sequence" = known_color_hex)

if (length(species_colors_hex) > 0) {
  species_colors_hex <- species_colors_hex[order(names(species_colors_hex))]
}

predicted_ecs <- unlist(lapply(
  tip_info$label[!tip_info$is_known],
  extract_last_ec_block
))

ec_levels <- sort(unique(predicted_ecs))
ec_colors <- make_ec_colors(ec_levels)

ec_colors_hex <- if (length(ec_colors) > 0) {
  setNames(vapply(ec_colors, hex_col, character(1)), names(ec_colors))
} else {
  character()
}
# =========================================================
# 9. 建树图
# =========================================================
if (USE_TOPOLOGY_ONLY) {
  p0 <- ggtree(tree, linewidth = TREE_LINEWIDTH, branch.length = "none")
} else {
  p0 <- ggtree(tree, linewidth = TREE_LINEWIDTH)
}

p0$data$x <- p0$data$x * TREE_X_SCALE
tree_data <- p0$data

tip_pos <- tree_data %>%
  dplyr::filter(isTip) %>%
  dplyr::select(label, y, x_tip = x) %>%
  dplyr::left_join(tip_info, by = "label")

domain_plot_df <- domain_df %>%
  dplyr::left_join(tip_pos %>% dplyr::select(label, y), by = "label") %>%
  dplyr::filter(!is.na(y))

motif_plot_df <- motif_df %>%
  dplyr::left_join(tip_pos %>% dplyr::select(label, y), by = "label") %>%
  dplyr::filter(!is.na(y))

seq_len_plot_df <- seq_len_df %>%
  dplyr::left_join(tip_pos %>% dplyr::select(label, y), by = "label") %>%
  dplyr::filter(!is.na(y))

# =========================================================
# 10. 布局参数
# =========================================================
tree_max_x <- max(tree_data$x, na.rm = TRUE)
seq_max_len <- max(seq_len_plot_df$seq_len, na.rm = TRUE)
if (!is.finite(seq_max_len) || seq_max_len <= 0) seq_max_len <- 100

max_label_chars <- max(nchar(tip_info$display_label), na.rm = TRUE)
BASE_UNIT <- max(tree_max_x, 1)

LABEL_OFFSET <- 0.05
TEXT_GAP <- max(4.0, max_label_chars * 0.035)

EC_PANEL_GAP <- max(0.8, BASE_UNIT * 0.06)
PANEL_GAP    <- max(0.8, BASE_UNIT * 0.06)

# EC 列单独放宽，给多个 EC 留空间
EC_PANEL_WIDTH <- max(5.5, BASE_UNIT * 0.35)
EC_PANEL_START <- tree_max_x + TEXT_GAP

# Domain 列从 EC 列右边重新开始
DOMAIN_PANEL_START <- EC_PANEL_START + EC_PANEL_WIDTH + EC_PANEL_GAP
DOMAIN_PANEL_WIDTH <- max(3.5, BASE_UNIT * 0.42)

# Motif 列再接在 Domain 后面
MOTIF_PANEL_START <- DOMAIN_PANEL_START + DOMAIN_PANEL_WIDTH + PANEL_GAP
MOTIF_PANEL_WIDTH <- max(3.8, BASE_UNIT * 0.42)

PLOT_X_MAX <- MOTIF_PANEL_START + MOTIF_PANEL_WIDTH


aa_to_domain_x <- function(pos) {
  DOMAIN_PANEL_START + (pos / seq_max_len) * DOMAIN_PANEL_WIDTH
}

aa_to_motif_x <- function(pos) {
  MOTIF_PANEL_START + (pos / seq_max_len) * MOTIF_PANEL_WIDTH
}

domain_backbone_df <- seq_len_plot_df %>%
  dplyr::mutate(
    x = aa_to_domain_x(1),
    xend = aa_to_domain_x(seq_len)
  )

motif_backbone_df <- seq_len_plot_df %>%
  dplyr::mutate(
    x = aa_to_motif_x(1),
    xend = aa_to_motif_x(seq_len)
  )

domain_plot_df <- domain_plot_df %>%
  dplyr::mutate(
    x = aa_to_domain_x(start),
    xend = aa_to_domain_x(end)
  )

motif_plot_df <- motif_plot_df %>%
  dplyr::mutate(
    x = aa_to_motif_x(start),
    xend = aa_to_motif_x(end)
  )

tip_label_df <- tip_pos %>%
  dplyr::mutate(
    x_label = x_tip + LABEL_OFFSET,
    x_ec = EC_PANEL_START + EC_PANEL_WIDTH / 2,
    x_guide_end = EC_PANEL_START - 0.15
  )

# =========================================================
# 11. 生成可点击 SVG 所需数据
# =========================================================
LINK_STEP_CODE <- "step6"

tip_label_df <- tip_label_df %>%
  dplyr::mutate(
    href = mapply(
      make_tip_href,
      label,
      is_known,
      MoreArgs = list(link_step_code = LINK_STEP_CODE)
    ),
    species_clean = clean_species(species)
  )

domain_plot_df <- domain_plot_df %>%
  dplyr::mutate(href = vapply(accession, make_domain_href, character(1)))

motif_plot_df <- motif_plot_df %>%
  dplyr::mutate(href = mapply(make_motif_href, label, motif))

edge_df <- build_edge_df(tree_data)

node_support_df <- tree_data %>%
  dplyr::filter(!isTip, !is.na(label), label != "")

tip_label_df$ec <- vapply(
  tip_label_df$label,
  extract_ec_display_text,
  character(1)
)

tip_label_df$ec_color <- ifelse(
  tip_label_df$is_known,
  known_color_hex,
  ec_colors_hex[tip_label_df$ec]
)

tip_label_df$ec_color[is.na(tip_label_df$ec_color) | tip_label_df$ec_color == ""] <- "#666666"
# =========================================================
# 12. 检查未匹配物种并加入 fallback 图例
# =========================================================
tip_label_df$tip_color <- ifelse(
  tip_label_df$is_known,
  known_color_hex,
  species_colors_hex[tip_label_df$species_clean]
)

tip_label_df$tip_color[is.na(tip_label_df$tip_color) | tip_label_df$tip_color == ""] <- "#666666"

unmatched_species <- sort(unique(
  tip_label_df$species_clean[
    !tip_label_df$is_known &
      !(tip_label_df$species_clean %in% names(species_colors_hex))
  ]
))

if (length(unmatched_species) > 0) {
  cat("以下 predicted species 未匹配到图例颜色，将使用灰色 Unmatched：\n")
  print(unmatched_species)
}

species_legend_hex <- species_colors_hex

if (any(tip_label_df$tip_color == "#666666" & !tip_label_df$is_known)) {
  species_legend_hex <- c(species_legend_hex, "Unmatched" = "#666666")
}

# 图例最多展示 10 个物种，其余用省略号表示
MAX_SPECIES_IN_LEGEND <- 10

if (length(species_legend_hex) > MAX_SPECIES_IN_LEGEND) {
  shown_names <- names(species_legend_hex)[1:MAX_SPECIES_IN_LEGEND]
  species_legend_hex <- species_legend_hex[shown_names]
  
  # 省略号本身给一个中性灰色，只表示“后面还有”
  species_legend_hex <- c(species_legend_hex, "..." = "#999999")
}

# =========================================================
# 13. SVG 版面参数
# =========================================================
PLOT_LEFT <- 40
TOP_PAD <- 90
BOTTOM_PAD <- 40
ROW_HEIGHT <- 28

BASE_PLOT_WIDTH <- 2000
PLOT_RIGHT <- PLOT_LEFT + BASE_PLOT_WIDTH
SVG_WIDTH  <- PLOT_RIGHT + 40

TITLE_Y <- 30
SUBTITLE_Y <- 60

LEGEND_ROW_GAP <- 28
LEGEND_ITEM_GAP <- 26
LEGEND_LEFT <- PLOT_LEFT
LEGEND_MAX_WIDTH <- SVG_WIDTH - 40

y_min <- min(tip_pos$y, na.rm = TRUE)
y_max <- max(tip_pos$y, na.rm = TRUE)

tree_height <- TOP_PAD + (y_max - y_min + 1) * ROW_HEIGHT + BOTTOM_PAD

estimate_legend_row_height <- function(title, items, max_width, x_start) {
  title_space <- max(120, nchar(title) * 10)
  xx <- x_start + title_space
  lines_n <- 1
  
  if (length(items) == 0) return(LEGEND_ROW_GAP)
  
  for (nm in items) {
    text_w <- max(40, nchar(nm) * 8.5)
    item_w <- 24 + text_w + LEGEND_ITEM_GAP
    
    if (xx + item_w > max_width) {
      lines_n <- lines_n + 1
      xx <- x_start + title_space
    }
    xx <- xx + item_w
  }
  
  lines_n * LEGEND_ROW_GAP
}

legend_h1 <- estimate_legend_row_height("Sequence type", names(known_legend_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
legend_h2 <- estimate_legend_row_height("EC", names(ec_colors_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
legend_h3 <- estimate_legend_row_height("Motif", names(motif_colors_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
legend_h4 <- estimate_legend_row_height("Domain", names(domain_colors_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
legend_h5 <- estimate_legend_row_height("Species", names(species_legend_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)

legend_total_height <- legend_h1 + legend_h2 + legend_h3 + legend_h4 + legend_h5 + 60
svg_height <- tree_height + legend_total_height

x_to_px <- function(x) {
  PLOT_LEFT + (x / PLOT_X_MAX) * (PLOT_RIGHT - PLOT_LEFT)
}

y_to_px <- function(y) {
  TOP_PAD + (y_max - y) * ROW_HEIGHT
}

# =========================================================
# 14. 转成 SVG 像素坐标
# =========================================================
edge_df <- edge_df %>%
  dplyr::mutate(
    px_parent = x_to_px(x_parent),
    py_parent = y_to_px(y_parent),
    px_child  = x_to_px(x_child),
    py_child  = y_to_px(y_child)
  )

tip_label_df <- tip_label_df %>%
  dplyr::mutate(
    px_tip = x_to_px(x_tip),
    py = y_to_px(y),
    px_label = x_to_px(x_label),
    px_ec = x_to_px(x_ec),
    px_guide_end = x_to_px(x_guide_end)
  )

node_support_df <- node_support_df %>%
  dplyr::mutate(
    px = x_to_px(x),
    py = y_to_px(y)
  )

domain_backbone_df <- domain_backbone_df %>%
  dplyr::mutate(
    px = x_to_px(x),
    pxend = x_to_px(xend),
    py = y_to_px(y)
  )

motif_backbone_df <- motif_backbone_df %>%
  dplyr::mutate(
    px = x_to_px(x),
    pxend = x_to_px(xend),
    py = y_to_px(y)
  )

domain_plot_df <- domain_plot_df %>%
  dplyr::mutate(
    px = x_to_px(x),
    pxend = x_to_px(xend),
    py = y_to_px(y)
  )

motif_plot_df <- motif_plot_df %>%
  dplyr::mutate(
    px = x_to_px(x),
    pxend = x_to_px(xend),
    py = y_to_px(y)
  )

# =========================================================
# 15. 开始写 SVG
# =========================================================
svg_lines <- character()

push_svg <- function(...) {
  svg_lines <<- c(svg_lines, paste0(...))
}

push_svg('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
push_svg(
  '<svg xmlns="http://www.w3.org/2000/svg" ',
  'xmlns:xlink="http://www.w3.org/1999/xlink" ',
  'width="', SVG_WIDTH, '" height="', svg_height, '" ',
  'viewBox="0 0 ', SVG_WIDTH, ' ', svg_height, '">'
)

push_svg('
<style>
.tree-edge{stroke:#000000;stroke-width:2;fill:none;}
.guide-line{stroke:#c8c8c8;stroke-width:1;stroke-dasharray:4 4;fill:none;}
.backbone{stroke:#c8c8c8;stroke-width:1;fill:none;}
.node-support{
  font-family:Arial, Helvetica, sans-serif;
  font-size:13px;
  fill:#444444;
  text-anchor:end;
}
.tip-label{
  font-family:Arial, Helvetica, sans-serif;
  font-size:15px;
  dominant-baseline:middle;
}
.tip-label-known{
  font-family:Arial, Helvetica, sans-serif;
  font-size:16px;
  font-weight:bold;
  dominant-baseline:middle;
}
.panel-title{
  font-family:Arial, Helvetica, sans-serif;
  font-size:18px;
  font-weight:bold;
  text-anchor:middle;
}
.main-title{
  font-family:Arial, Helvetica, sans-serif;
  font-size:24px;
  font-weight:bold;
  text-anchor:middle;
}
.legend-title{
  font-family:Arial, Helvetica, sans-serif;
  font-size:14px;
  font-weight:bold;
}
.legend-text{
  font-family:Arial, Helvetica, sans-serif;
  font-size:13px;
  dominant-baseline:middle;
}
a:hover text { text-decoration: underline; }
a:hover rect { opacity: 0.78; }
a:hover line.clickable-line { opacity: 0.78; }
</style>
')

push_svg(
  '<text class="main-title" x="', SVG_WIDTH / 2, '" y="', TITLE_Y, '">',
  svg_escape("Phylogenetic tree with domain architecture and motif architecture"),
  '</text>'
)
push_svg(
  '<text class="panel-title" x="', x_to_px(EC_PANEL_START + EC_PANEL_WIDTH / 2),
  '" y="', SUBTITLE_Y, '">', svg_escape("EC"), '</text>'
)

push_svg(
  '<text class="panel-title" x="', x_to_px(DOMAIN_PANEL_START + DOMAIN_PANEL_WIDTH / 2),
  '" y="', SUBTITLE_Y, '">', svg_escape("Domain architecture"), '</text>'
)
push_svg(
  '<text class="panel-title" x="', x_to_px(MOTIF_PANEL_START + MOTIF_PANEL_WIDTH / 2),
  '" y="', SUBTITLE_Y, '">', svg_escape("Motif architecture"), '</text>'
)

# ========= 树分支 =========
for (i in seq_len(nrow(edge_df))) {
  r <- edge_df[i, ]
  
  push_svg(
    '<line class="tree-edge" x1="', sprintf("%.2f", r$px_parent),
    '" y1="', sprintf("%.2f", r$py_parent),
    '" x2="', sprintf("%.2f", r$px_parent),
    '" y2="', sprintf("%.2f", r$py_child), '"/>'
  )
  
  push_svg(
    '<line class="tree-edge" x1="', sprintf("%.2f", r$px_parent),
    '" y1="', sprintf("%.2f", r$py_child),
    '" x2="', sprintf("%.2f", r$px_child),
    '" y2="', sprintf("%.2f", r$py_child), '"/>'
  )
}

# ========= bootstrap =========
if (SHOW_NODE_SUPPORT && nrow(node_support_df) > 0) {
  for (i in seq_len(nrow(node_support_df))) {
    r <- node_support_df[i, ]
    push_svg(
      '<text class="node-support" x="', sprintf("%.2f", r$px - 8),
      '" y="', sprintf("%.2f", r$py - 4), '">',
      svg_escape(r$label),
      '</text>'
    )
  }
}

# ========= tip 到右侧的虚线 =========
for (i in seq_len(nrow(tip_label_df))) {
  r <- tip_label_df[i, ]
  push_svg(
    '<line class="guide-line" x1="', sprintf("%.2f", r$px_tip),
    '" y1="', sprintf("%.2f", r$py),
    '" x2="', sprintf("%.2f", r$px_guide_end),
    '" y2="', sprintf("%.2f", r$py), '"/>'
  )
}

# ========= tip 标签 =========
if (SHOW_TIP_LABEL && nrow(tip_label_df) > 0) {
  for (i in seq_len(nrow(tip_label_df))) {
    r <- tip_label_df[i, ]
    
    if (isTRUE(r$is_known)) {
      label_body <- paste0(
        '<g>',
        '<circle cx="', sprintf("%.2f", r$px_tip), '" cy="', sprintf("%.2f", r$py),
        '" r="3.2" fill="', known_color_hex, '"/>',
        '<text class="tip-label-known" x="', sprintf("%.2f", r$px_label),
        '" y="', sprintf("%.2f", r$py), '" fill="', known_color_hex, '">',
        svg_escape(r$display_label),
        '</text></g>'
      )
    } else {
      label_body <- paste0(
        '<g>',
        '<text class="tip-label" x="', sprintf("%.2f", r$px_label),
        '" y="', sprintf("%.2f", r$py), '" fill="', r$tip_color, '">',
        svg_escape(r$display_label),
        '</text></g>'
      )
    }
    
    push_svg(wrap_svg_link(label_body, r$href))
  }
}

# ========= Domain 骨架 =========
for (i in seq_len(nrow(domain_backbone_df))) {
  r <- domain_backbone_df[i, ]
  push_svg(
    '<line class="backbone" x1="', sprintf("%.2f", r$px),
    '" y1="', sprintf("%.2f", r$py),
    '" x2="', sprintf("%.2f", r$pxend),
    '" y2="', sprintf("%.2f", r$py), '"/>'
  )
}

# ========= EC 列 =========
# Known: 整个 EC 文本统一红色
# Predicted: 同一个 EC 号同一个颜色；多个 EC 号拆开分别着色后横向排布

EC_FONT_SIZE <- 14
EC_TEXT_GAP  <- 6

if (nrow(tip_label_df) > 0) {
  for (i in seq_len(nrow(tip_label_df))) {
    r <- tip_label_df[i, ]
    ecs <- extract_last_ec_block(r$label)
    
    if (length(ecs) == 0) next
    
    if (isTRUE(r$is_known)) {
      ec_text <- paste(ecs, collapse = "; ")
      push_svg(
        '<text x="', sprintf("%.2f", r$px_ec),
        '" y="', sprintf("%.2f", r$py),
        '" font-family="Arial, Helvetica, sans-serif" font-size="', EC_FONT_SIZE,
        '" fill="', known_color_hex,
        '" text-anchor="middle" dominant-baseline="middle" font-weight="bold">',
        svg_escape(ec_text),
        '</text>'
      )
    } else {
      # 多个 predicted EC 分段着色
      # 先简单按字符宽度估计整体宽度，使其以 px_ec 为中心
      piece_widths <- sapply(ecs, function(ec) max(28, nchar(ec) * 8.2))
      sep_width <- 10
      total_width <- sum(piece_widths) + sep_width * (length(ecs) - 1)
      x_cursor <- r$px_ec - total_width / 2
      
      for (j in seq_along(ecs)) {
        ec_now <- ecs[j]
        ec_col <- ec_colors_hex[[ec_now]]
        if (is.null(ec_col) || is.na(ec_col) || ec_col == "") ec_col <- "#666666"
        
        text_x <- x_cursor + piece_widths[j] / 2
        
        push_svg(
          '<text x="', sprintf("%.2f", text_x),
          '" y="', sprintf("%.2f", r$py),
          '" font-family="Arial, Helvetica, sans-serif" font-size="', EC_FONT_SIZE,
          '" fill="', ec_col,
          '" text-anchor="middle" dominant-baseline="middle">',
          svg_escape(ec_now),
          '</text>'
        )
        
        x_cursor <- x_cursor + piece_widths[j]
        
        if (j < length(ecs)) {
          push_svg(
            '<text x="', sprintf("%.2f", x_cursor + sep_width / 2),
            '" y="', sprintf("%.2f", r$py),
            '" font-family="Arial, Helvetica, sans-serif" font-size="', EC_FONT_SIZE,
            '" fill="#666666" text-anchor="middle" dominant-baseline="middle">',
            svg_escape(";"),
            '</text>'
          )
          x_cursor <- x_cursor + sep_width
        }
      }
    }
  }
}
# ========= Domain block =========
DOMAIN_BOX_H <- 10

if (nrow(domain_plot_df) > 0) {
  for (i in seq_len(nrow(domain_plot_df))) {
    r <- domain_plot_df[i, ]
    fill_col <- domain_colors_hex[[r$domain]]
    if (is.null(fill_col) || is.na(fill_col) || fill_col == "") fill_col <- "#999999"
    
    rect_body <- paste0(
      '<rect x="', sprintf("%.2f", r$px),
      '" y="', sprintf("%.2f", r$py - DOMAIN_BOX_H / 2),
      '" width="', sprintf("%.2f", max(1, r$pxend - r$px)),
      '" height="', DOMAIN_BOX_H,
      '" rx="2" ry="2" fill="', fill_col,
      '" stroke="', fill_col, '"/>'
    )
    push_svg(wrap_svg_link(rect_body, r$href))
  }
}

# ========= Motif 骨架 =========
for (i in seq_len(nrow(motif_backbone_df))) {
  r <- motif_backbone_df[i, ]
  push_svg(
    '<line class="backbone" x1="', sprintf("%.2f", r$px),
    '" y1="', sprintf("%.2f", r$py),
    '" x2="', sprintf("%.2f", r$pxend),
    '" y2="', sprintf("%.2f", r$py), '"/>'
  )
}

# ========= Motif block =========
MOTIF_BOX_H <- 10

if (nrow(motif_plot_df) > 0) {
  for (i in seq_len(nrow(motif_plot_df))) {
    r <- motif_plot_df[i, ]
    fill_col <- motif_colors_hex[[r$motif]]
    if (is.null(fill_col) || is.na(fill_col) || fill_col == "") fill_col <- "#999999"
    
    rect_body <- paste0(
      '<rect x="', sprintf("%.2f", r$px),
      '" y="', sprintf("%.2f", r$py - MOTIF_BOX_H / 2),
      '" width="', sprintf("%.2f", max(1, r$pxend - r$px)),
      '" height="', MOTIF_BOX_H,
      '" rx="2" ry="2" fill="', fill_col,
      '" stroke="', fill_col, '"/>'
    )
    push_svg(wrap_svg_link(rect_body, r$href))
  }
}

# =========================================================
# 16. 自动换行图例
# =========================================================
draw_legend_row <- function(title, items, colors, x_start, y, shape = "rect",
                            max_width = SVG_WIDTH - 40, line_gap = 22) {
  xx <- x_start
  yy <- y
  
  push_svg(
    '<text class="legend-title" x="', xx, '" y="', yy, '">',
    svg_escape(title), '</text>'
  )
  
  title_space <- max(120, nchar(title) * 10)
  xx <- xx + title_space
  
  if (length(items) == 0) return(invisible(yy))
  
  for (nm in items) {
    col_now <- colors[[nm]]
    text_w <- max(40, nchar(nm) * 8.5)
    item_w <- 24 + text_w + LEGEND_ITEM_GAP
    
    if (xx + item_w > max_width) {
      xx <- x_start + title_space
      yy <- yy + line_gap
    }
    
    if (shape == "circle") {
      push_svg(
        '<circle cx="', xx + 8, '" cy="', yy - 5,
        '" r="5.5" fill="', col_now,
        '" stroke="', col_now, '"/>'
      )
    } else {
      push_svg(
        '<rect x="', xx, '" y="', yy - 10,
        '" width="16" height="10" fill="', col_now,
        '" stroke="', col_now, '"/>'
      )
    }
    
    push_svg(
      '<text class="legend-text" x="', xx + 24, '" y="', yy - 2, '">',
      svg_escape(nm), '</text>'
    )
    
    xx <- xx + item_w
  }
  
  invisible(yy)
}

legend_y <- tree_height + 22

legend_y <- draw_legend_row(
  "Sequence type",
  names(known_legend_hex),
  known_legend_hex,
  LEGEND_LEFT,
  legend_y,
  shape = "circle",
  max_width = LEGEND_MAX_WIDTH,
  line_gap = LEGEND_ROW_GAP
) + LEGEND_ROW_GAP

legend_y <- draw_legend_row(
  "EC",
  names(ec_colors_hex),
  ec_colors_hex,
  LEGEND_LEFT,
  legend_y,
  shape = "rect",
  max_width = LEGEND_MAX_WIDTH,
  line_gap = LEGEND_ROW_GAP
) + LEGEND_ROW_GAP

legend_y <- draw_legend_row(
  "Domain",
  names(domain_colors_hex),
  domain_colors_hex,
  LEGEND_LEFT,
  legend_y,
  shape = "rect",
  max_width = LEGEND_MAX_WIDTH,
  line_gap = LEGEND_ROW_GAP
) + LEGEND_ROW_GAP

legend_y <- draw_legend_row(
  "Motif",
  names(motif_colors_hex),
  motif_colors_hex,
  LEGEND_LEFT,
  legend_y,
  shape = "rect",
  max_width = LEGEND_MAX_WIDTH,
  line_gap = LEGEND_ROW_GAP
) + LEGEND_ROW_GAP


legend_y <- draw_legend_row(
  "Species",
  names(species_legend_hex),
  species_legend_hex,
  LEGEND_LEFT,
  legend_y,
  shape = "rect",
  max_width = LEGEND_MAX_WIDTH,
  line_gap = LEGEND_ROW_GAP
)

push_svg('</svg>')

writeLines(svg_lines, out_svg, useBytes = TRUE)

cat("SVG 已输出：", out_svg, "\n")
cat("Domain CSV 已输出：", out_domain_csv, "\n")
cat("Motif CSV 已输出：", out_motif_csv, "\n")
