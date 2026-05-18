library(ape)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(xml2)
library(ggtree)
library(scales)

SHOW_TIP_LABEL    <- TRUE
SHOW_NODE_SUPPORT <- FALSE
SUPPORT_MODE      <- "full"
SUPPORT_MIN       <- NA
TREE_LINEWIDTH    <- 1.0
USE_TOPOLOGY_ONLY <- FALSE
MEME_POSITION_ZERO_BASED <- TRUE
KNOWN_COLOR <- "#C0392B"
MAX_SPECIES_IN_LEGEND <- 10

get_step_files <- function(step_dir) {
  step_name <- basename(step_dir)
  
  folder_path <- file.path(step_dir, "draw_tree")
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  list(
    step_dir = step_dir,
    step_name = step_name,
    tree_file = file.path(step_dir, "merged_aligned_trimmed.fasta.treefile"),
    domain_file = file.path(step_dir, "domain_hits.txt"),
    meme_file = file.path(step_dir, "meme_results", "meme.xml"),
    protein_fasta = file.path(step_dir, "merged.fasta"),
    folder_path = folder_path,
    out_domain_csv = file.path(folder_path, "parsed_domain_table1.csv"),
    out_motif_csv = file.path(folder_path, "parsed_motif_table1.csv")
  )
}

check_required_files <- function(files) {
  required <- c("tree_file", "domain_file", "meme_file")
  for (nm in required) {
    f <- files[[nm]]
    if (is.null(f) || is.na(f) || !file.exists(f)) {
      stop("Required file not found: ", nm, " = ", f)
    }
  }
}

expand_step_codes <- function(step_name) {
  nums <- stringr::str_extract_all(step_name, "\\d+")[[1]]
  if (length(nums) == 0) return(step_name)
  paste0("step", nums)
}

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

clean_pfam_accession <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\..*$", "", x)
  x
}

hex_col <- function(x) {
  rgb_mat <- grDevices::col2rgb(x)
  grDevices::rgb(
    rgb_mat[1, ] / 255,
    rgb_mat[2, ] / 255,
    rgb_mat[3, ] / 255
  )
}

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
  paste0(
    '<a xlink:href="',
    svg_escape(href),
    '" target="_blank">',
    body,
    '</a>'
  )
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

make_motif_href <- function(label, motif, link_step_code = NULL) {
  if (is.null(link_step_code) || is.na(link_step_code) || link_step_code == "") {
    return(NA_character_)
  }
  
  paste0(
    "/static/pathways/berberine/all_step_phylogeny_svg/",
    link_step_code, 
    "_meme.html"
  )
}

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
    mutate(
      label = clean_label(target_name),
      start = as.numeric(ali_from),
      end = as.numeric(ali_to),
      domain = trimws(query_name),
      accession = clean_pfam_accession(query_accession),
      domain_score = as.numeric(domain_score),
      i_Evalue = as.numeric(i_Evalue)
    )
  
  out %>%
    filter(
      !is.na(label), label != "",
      !is.na(start), !is.na(end),
      !is.na(domain), domain != "",
      !is.na(i_Evalue),
      !is.na(domain_score),
      domain_score >= min_score,
      i_Evalue <= max_evalue
    ) %>%
    select(label, start, end, domain, accession, domain_score, i_Evalue) %>%
    distinct()
}

resolve_domain_overlap <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  df <- df %>% arrange(label, start, i_Evalue)
  result <- list()
  
  for (seq_id in unique(df$label)) {
    sub_df <- df[df$label == seq_id, ]
    keep <- rep(TRUE, nrow(sub_df))
    
    for (i in seq_len(nrow(sub_df))) {
      if (!keep[i]) next
      
      for (j in seq_len(nrow(sub_df))) {
        if (i == j || !keep[j]) next
        
        overlap <- !(sub_df$end[i] < sub_df$start[j] ||
                       sub_df$end[j] < sub_df$start[i])
        
        if (overlap) {
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

read_meme_xml <- function(meme_file, zero_based = TRUE) {
  if (is.na(meme_file) || !file.exists(meme_file)) {
    return(data.frame(label = character(), motif = character(), start = numeric(), end = numeric()))
  }
  
  doc <- read_xml(meme_file)
  
  seq_nodes <- xml_find_all(doc, ".//sequence")
  seq_map <- data.frame(
    sequence_id = xml_attr(seq_nodes, "id"),
    label = xml_attr(seq_nodes, "name"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(label = clean_label(label)) %>%
    filter(!is.na(sequence_id), sequence_id != "", !is.na(label), label != "") %>%
    distinct(sequence_id, .keep_all = TRUE)
  
  motif_nodes <- xml_find_all(doc, ".//motif")
  motif_info <- data.frame(
    motif_id = xml_attr(motif_nodes, "id"),
    motif_name = xml_attr(motif_nodes, "name"),
    motif_alt = xml_attr(motif_nodes, "alt"),
    width = suppressWarnings(as.numeric(xml_attr(motif_nodes, "width"))),
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
      
      seq_id <- xml_attr(node, "sequence_id")
      pos_raw <- suppressWarnings(as.numeric(xml_attr(node, "position")))
      
      motif_node <- xml_find_first(node, "ancestor::motif[1]")
      motif_id <- xml_attr(motif_node, "id")
      motif_name <- xml_attr(motif_node, "name")
      motif_alt <- xml_attr(motif_node, "alt")
      motif_w <- suppressWarnings(as.numeric(xml_attr(motif_node, "width")))
      
      motif_display <- ifelse(
        !is.na(motif_alt) && motif_alt != "",
        motif_alt,
        ifelse(!is.na(motif_name) && motif_name != "", motif_name, motif_id)
      )
      
      if (is.na(pos_raw) || is.na(motif_w)) next
      
      start <- if (zero_based) pos_raw + 1 else pos_raw
      end <- start + motif_w - 1
      
      res_list[[i]] <- data.frame(
        sequence_id = seq_id,
        motif = motif_display,
        start = as.numeric(start),
        end = as.numeric(end),
        stringsAsFactors = FALSE
      )
    }
    
    return(
      bind_rows(res_list) %>%
        left_join(seq_map, by = "sequence_id") %>%
        select(label, motif, start, end) %>%
        filter(
          !is.na(label), label != "",
          !is.na(motif), motif != "",
          !is.na(start), !is.na(end)
        ) %>%
        distinct()
    )
  }
  
  scanned_site_nodes <- xml_find_all(doc, ".//scanned_site")
  
  if (length(scanned_site_nodes) > 0) {
    res_list <- vector("list", length(scanned_site_nodes))
    
    for (i in seq_along(scanned_site_nodes)) {
      node <- scanned_site_nodes[[i]]
      
      motif_id <- xml_attr(node, "motif_id") %||% xml_attr(node, "motif")
      pos_raw <- suppressWarnings(as.numeric(xml_attr(node, "position") %||% xml_attr(node, "pos")))
      
      parent_scanned_sites <- xml_find_first(node, "ancestor::scanned_sites[1]")
      seq_id <- xml_attr(parent_scanned_sites, "sequence_id")
      
      motif_row <- motif_info %>% filter(motif_id == !!motif_id)
      motif_display <- motif_row$motif_display[1] %||% motif_id
      motif_w <- suppressWarnings(as.numeric(motif_row$width[1]))
      
      if (is.na(pos_raw) || is.na(motif_w)) next
      
      start <- if (zero_based) pos_raw + 1 else pos_raw
      end <- start + motif_w - 1
      
      res_list[[i]] <- data.frame(
        sequence_id = seq_id,
        motif = motif_display,
        start = as.numeric(start),
        end = as.numeric(end),
        stringsAsFactors = FALSE
      )
    }
    
    return(
      bind_rows(res_list) %>%
        left_join(seq_map, by = "sequence_id") %>%
        select(label, motif, start, end) %>%
        filter(
          !is.na(label), label != "",
          !is.na(motif), motif != "",
          !is.na(start), !is.na(end)
        ) %>%
        distinct()
    )
  }
  
  warning("No recognizable motif site structure found in meme.xml")
  data.frame(label = character(), motif = character(), start = numeric(), end = numeric())
}

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

build_edge_df <- function(tree_data) {
  child_df <- tree_data %>%
    filter(parent != node) %>%
    select(node, parent, x_child = x, y_child = y)
  
  parent_df <- tree_data %>%
    select(parent_node = node, x_parent = x, y_parent = y)
  
  child_df %>%
    left_join(parent_df, by = c("parent" = "parent_node"))
}

make_many_species_colors <- function(species_levels) {
  n <- length(species_levels)
  if (n == 0) return(c())
  
  if (n <= 12) {
    cols <- hcl.colors(n, palette = "Dark 3")
    return(setNames(cols, species_levels))
  }
  
  hues <- seq(15, 375, length.out = n + 1)[1:n]
  l_vals <- rep(c(55, 70), length.out = n)
  c_vals <- rep(c(100, 65), length.out = n)
  
  cols <- grDevices::hcl(
    h = hues,
    c = c_vals,
    l = l_vals
  )
  
  setNames(cols, species_levels)
}

auto_layout_params <- function(tree_data_raw,
                               tip_info,
                               seq_len_plot_df,
                               domain_df,
                               motif_df) {
  n_tip <- sum(tree_data_raw$isTip, na.rm = TRUE)
  max_label_chars <- max(nchar(tip_info$display_label), na.rm = TRUE)
  tree_max_x <- max(tree_data_raw$x, na.rm = TRUE)
  
  if (!is.finite(n_tip) || n_tip <= 0) n_tip <- 20
  if (!is.finite(max_label_chars)) max_label_chars <- 20
  if (!is.finite(tree_max_x) || tree_max_x <= 0) tree_max_x <- 1
  
  row_height <- if (n_tip <= 30) {
    28
  } else if (n_tip <= 80) {
    22
  } else {
    18
  }
  
  tip_font_size <- if (n_tip <= 30) {
    15
  } else if (n_tip <= 80) {
    12
  } else {
    10
  }
  
  tree_x_scale <- if (tree_max_x < 1) {
    3.2
  } else {
    2.5
  }
  
  domain_panel_width <- max(2.2, tree_max_x * tree_x_scale * 0.42)
  motif_panel_width <- max(2.8, tree_max_x * tree_x_scale * 0.42)
  
  base_plot_width <- if (n_tip <= 30) {
    1700
  } else if (n_tip <= 80) {
    1900
  } else {
    2200
  }
  
  list(
    TREE_X_SCALE = tree_x_scale,
    DOMAIN_PANEL_WIDTH = domain_panel_width,
    MOTIF_PANEL_WIDTH = motif_panel_width,
    ROW_HEIGHT = row_height,
    TIP_FONT_SIZE = tip_font_size,
    BASE_PLOT_WIDTH = base_plot_width
  )
}

draw_one_step <- function(step_dir) {
  
  files <- get_step_files(step_dir)
  check_required_files(files)
  
  step_name <- files$step_name
  step_codes_for_links <- expand_step_codes(step_name)
  
  message("Processing: ", step_name)
  message("Link step codes: ", paste(step_codes_for_links, collapse = ", "))
  
  tree <- read.tree(files$tree_file)
  
  if (is.null(tree$tip.label) || length(tree$tip.label) == 0) {
    stop("treefile has no tip labels")
  }
  
  tree$tip.label <- clean_label(tree$tip.label)
  tip_labels <- tree$tip.label
  
  tip_info <- data.frame(
    label = tip_labels,
    is_known = grepl("^\\*?Known\\|", tip_labels),
    stringsAsFactors = FALSE
  )
  
  tip_info$species <- sapply(tip_info$label, function(x) {
    parts <- strsplit(x, "\\|")[[1]]
    if (length(parts) >= 2) {
      return(clean_species(parts[2]))
    } else {
      return("Unknown")
    }
  })
  
  tip_info$seq_type <- ifelse(tip_info$is_known, "Known", "Predicted")
  
  tip_info$known_id <- ifelse(
    tip_info$is_known,
    sub("^\\*?Known\\|([^|]+)\\|.*$", "\\1", tip_info$label),
    NA
  )
  
  tip_info$display_label <- tip_info$label
  
  species_levels <- tip_info %>%
    filter(!is_known) %>%
    pull(species) %>%
    unique() %>%
    sort()
  
  species_colors <- make_many_species_colors(species_levels)
  
  tree$node.label <- get_support_labels(
    node_labels_raw = tree$node.label,
    mode = SUPPORT_MODE,
    min_value = SUPPORT_MIN
  )
  
  domain_df <- read_domain_table(files$domain_file)
  motif_df <- read_meme_xml(files$meme_file, zero_based = MEME_POSITION_ZERO_BASED)
  
  domain_df$label <- clean_label(domain_df$label)
  motif_df$label <- clean_label(motif_df$label)
  
  domain_df <- domain_df %>% filter(label %in% tree$tip.label)
  motif_df <- motif_df %>% filter(label %in% tree$tip.label)
  
  domain_df <- resolve_domain_overlap(domain_df)
  
  write.csv(domain_df, files$out_domain_csv, row.names = FALSE)
  write.csv(motif_df, files$out_motif_csv, row.names = FALSE)
  
  seq_len_df <- read_fasta_lengths(files$protein_fasta)
  
  if (is.null(seq_len_df)) {
    seq_len_df <- bind_rows(
      domain_df %>% transmute(label, pos = end),
      motif_df %>% transmute(label, pos = end)
    ) %>%
      group_by(label) %>%
      summarise(seq_len = max(pos, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(seq_len_df) || nrow(seq_len_df) == 0) {
    seq_len_df <- data.frame(
      label = tree$tip.label,
      seq_len = 100,
      stringsAsFactors = FALSE
    )
  } else {
    max_len_all <- max(seq_len_df$seq_len, na.rm = TRUE)
    if (!is.finite(max_len_all)) max_len_all <- 100
    
    seq_len_df <- full_join(
      data.frame(label = tree$tip.label, stringsAsFactors = FALSE),
      seq_len_df,
      by = "label"
    )
    
    seq_len_df$seq_len[is.na(seq_len_df$seq_len)] <- max_len_all
  }
  
  if (USE_TOPOLOGY_ONLY) {
    p0 <- ggtree(tree, linewidth = TREE_LINEWIDTH, branch.length = "none")
  } else {
    p0 <- ggtree(tree, linewidth = TREE_LINEWIDTH)
  }
  
  tree_data_raw <- p0$data
  
  tip_pos_raw <- tree_data_raw %>%
    filter(isTip) %>%
    select(label, y, x_tip = x) %>%
    left_join(tip_info, by = "label")
  
  seq_len_plot_df_raw <- seq_len_df %>%
    left_join(tip_pos_raw %>% select(label, y), by = "label") %>%
    filter(!is.na(y))
  
  layout <- auto_layout_params(
    tree_data_raw = tree_data_raw,
    tip_info = tip_info,
    seq_len_plot_df = seq_len_plot_df_raw,
    domain_df = domain_df,
    motif_df = motif_df
  )
  
  TREE_X_SCALE <- layout$TREE_X_SCALE
  DOMAIN_PANEL_WIDTH <- layout$DOMAIN_PANEL_WIDTH
  MOTIF_PANEL_WIDTH <- layout$MOTIF_PANEL_WIDTH
  ROW_HEIGHT <- layout$ROW_HEIGHT
  TIP_FONT_SIZE <- layout$TIP_FONT_SIZE
  BASE_PLOT_WIDTH <- layout$BASE_PLOT_WIDTH
  
  p0$data$x <- p0$data$x * TREE_X_SCALE
  tree_data <- p0$data
  
  tip_pos <- tree_data %>%
    filter(isTip) %>%
    select(label, y, x_tip = x) %>%
    left_join(tip_info, by = "label")
  
  domain_plot_df <- domain_df %>%
    left_join(tip_pos %>% select(label, y), by = "label") %>%
    filter(!is.na(y))
  
  motif_plot_df <- motif_df %>%
    left_join(tip_pos %>% select(label, y), by = "label") %>%
    filter(!is.na(y))
  
  seq_len_plot_df <- seq_len_df %>%
    left_join(tip_pos %>% select(label, y), by = "label") %>%
    filter(!is.na(y))
  
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
  
  tree_max_x <- max(tree_data$x, na.rm = TRUE)
  seq_max_len <- max(seq_len_plot_df$seq_len, na.rm = TRUE)
  if (!is.finite(seq_max_len) || seq_max_len <= 0) seq_max_len <- 100
  
  LABEL_OFFSET <- 0.05
  
  max_tip_x <- max(tip_pos$x_tip, na.rm = TRUE)
  max_label_chars <- max(nchar(tip_pos$display_label), na.rm = TRUE)
  
  if (!is.finite(max_tip_x)) max_tip_x <- tree_max_x
  if (!is.finite(max_label_chars)) max_label_chars <- 20
  
  label_px_width <- max_label_chars * TIP_FONT_SIZE * 0.58
  initial_coord_width <- max_tip_x + DOMAIN_PANEL_WIDTH + MOTIF_PANEL_WIDTH + 3
  label_coord_width <- (label_px_width / BASE_PLOT_WIDTH) * initial_coord_width
  
  LABEL_RIGHT_PADDING_COORD <- 1.2
  
  DOMAIN_PANEL_START <- max_tip_x +
    LABEL_OFFSET +
    label_coord_width +
    LABEL_RIGHT_PADDING_COORD
  
  PANEL_GAP <- max(0.35, DOMAIN_PANEL_WIDTH * 0.12)
  
  MOTIF_PANEL_START <- DOMAIN_PANEL_START + DOMAIN_PANEL_WIDTH + PANEL_GAP
  PLOT_X_MAX <- MOTIF_PANEL_START + MOTIF_PANEL_WIDTH
  
  aa_to_domain_x <- function(pos) {
    DOMAIN_PANEL_START + (pos / seq_max_len) * DOMAIN_PANEL_WIDTH
  }
  
  aa_to_motif_x <- function(pos) {
    MOTIF_PANEL_START + (pos / seq_max_len) * MOTIF_PANEL_WIDTH
  }
  
  domain_backbone_df <- seq_len_plot_df %>%
    mutate(
      x = aa_to_domain_x(1),
      xend = aa_to_domain_x(seq_len)
    )
  
  motif_backbone_df <- seq_len_plot_df %>%
    mutate(
      x = aa_to_motif_x(1),
      xend = aa_to_motif_x(seq_len)
    )
  
  domain_plot_df <- domain_plot_df %>%
    mutate(
      x = aa_to_domain_x(start),
      xend = aa_to_domain_x(end),
      href = vapply(accession, make_domain_href, character(1))
    )
  
  motif_plot_df <- motif_plot_df %>%
    mutate(
      x = aa_to_motif_x(start),
      xend = aa_to_motif_x(end)
    )
  
  edge_df <- build_edge_df(tree_data)
  
  node_support_df <- tree_data %>%
    filter(!isTip, !is.na(label), label != "")
  
  PLOT_LEFT <- 40
  TOP_PAD <- 90
  BOTTOM_PAD <- 40
  
  PLOT_RIGHT <- PLOT_LEFT + BASE_PLOT_WIDTH
  SVG_WIDTH <- PLOT_RIGHT + 40
  
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
  
  result_list <- list()
  
  for (link_step_code in step_codes_for_links) {
    
    message("Generating SVG for link step: ", link_step_code)
    
    motif_plot_df_current <- motif_plot_df %>%
      mutate(
        href = mapply(
          make_motif_href,
          label,
          motif,
          MoreArgs = list(link_step_code = link_step_code)
        )
      )
    out_svg <- file.path(
      files$folder_path,
      paste0(link_step_code, ".svg")
    )
    
    tip_label_df <- tip_pos %>%
      mutate(
        x_label = x_tip + LABEL_OFFSET,
        x_guide_end = DOMAIN_PANEL_START - 0.25,
        href = mapply(
          make_tip_href,
          label,
          is_known,
          MoreArgs = list(link_step_code = link_step_code)
        ),
        species_clean = clean_species(species)
      )
    
    tip_label_df$tip_color <- ifelse(
      tip_label_df$is_known,
      known_color_hex,
      species_colors_hex[tip_label_df$species_clean]
    )
    
    tip_label_df$tip_color[is.na(tip_label_df$tip_color) | tip_label_df$tip_color == ""] <- "#666666"
    
    species_legend_hex <- species_colors_hex
    
    if (any(tip_label_df$tip_color == "#666666" & !tip_label_df$is_known)) {
      species_legend_hex <- c(species_legend_hex, "Unmatched" = "#666666")
    }
    
    if (length(species_legend_hex) > MAX_SPECIES_IN_LEGEND) {
      shown_names <- names(species_legend_hex)[1:MAX_SPECIES_IN_LEGEND]
      species_legend_hex <- species_legend_hex[shown_names]
      species_legend_hex <- c(species_legend_hex, "..." = "#999999")
    }
    
    legend_h1 <- estimate_legend_row_height("Sequence type", names(known_legend_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
    legend_h2 <- estimate_legend_row_height("Motif", names(motif_colors_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
    legend_h3 <- estimate_legend_row_height("Domain", names(domain_colors_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
    legend_h4 <- estimate_legend_row_height("Species", names(species_legend_hex), LEGEND_MAX_WIDTH, LEGEND_LEFT)
    
    legend_total_height <- legend_h1 + legend_h2 + legend_h3 + legend_h4 + 40
    svg_height <- tree_height + legend_total_height
    
    x_to_px <- function(x) {
      PLOT_LEFT + (x / PLOT_X_MAX) * (PLOT_RIGHT - PLOT_LEFT)
    }
    
    y_to_px <- function(y) {
      TOP_PAD + (y_max - y) * ROW_HEIGHT
    }
    
    edge_df_px <- edge_df %>%
      mutate(
        px_parent = x_to_px(x_parent),
        py_parent = y_to_px(y_parent),
        px_child = x_to_px(x_child),
        py_child = y_to_px(y_child)
      )
    
    tip_label_df <- tip_label_df %>%
      mutate(
        px_tip = x_to_px(x_tip),
        py = y_to_px(y),
        px_label = x_to_px(x_label),
        px_guide_end = x_to_px(x_guide_end)
      )
    
    node_support_df_px <- node_support_df %>%
      mutate(
        px = x_to_px(x),
        py = y_to_px(y)
      )
    
    domain_backbone_df_px <- domain_backbone_df %>%
      mutate(
        px = x_to_px(x),
        pxend = x_to_px(xend),
        py = y_to_px(y)
      )
    
    motif_backbone_df_px <- motif_backbone_df %>%
      mutate(
        px = x_to_px(x),
        pxend = x_to_px(xend),
        py = y_to_px(y)
      )
    
    domain_plot_df_px <- domain_plot_df %>%
      mutate(
        px = x_to_px(x),
        pxend = x_to_px(xend),
        py = y_to_px(y)
      )
    
    motif_plot_df_px <- motif_plot_df_current %>%
      mutate(
        px = x_to_px(x),
        pxend = x_to_px(xend),
        py = y_to_px(y)
      )
    
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
    
    push_svg(paste0('
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
  font-size:', TIP_FONT_SIZE, 'px;
  dominant-baseline:middle;
}
.tip-label-known{
  font-family:Arial, Helvetica, sans-serif;
  font-size:', TIP_FONT_SIZE + 1, 'px;
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
'))
    
    push_svg(
      '<text class="main-title" x="', SVG_WIDTH / 2, '" y="', TITLE_Y, '">',
      svg_escape(paste0("Phylogenetic tree with domain and motif architecture - ", link_step_code)),
      '</text>'
    )
    
    push_svg(
      '<text class="panel-title" x="', x_to_px(DOMAIN_PANEL_START + DOMAIN_PANEL_WIDTH / 2),
      '" y="', SUBTITLE_Y, '">', svg_escape("Domain architecture"), '</text>'
    )
    
    push_svg(
      '<text class="panel-title" x="', x_to_px(MOTIF_PANEL_START + MOTIF_PANEL_WIDTH / 2),
      '" y="', SUBTITLE_Y, '">', svg_escape("Motif architecture"), '</text>'
    )
    
    for (i in seq_len(nrow(edge_df_px))) {
      r <- edge_df_px[i, ]
      
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
    
    if (SHOW_NODE_SUPPORT && nrow(node_support_df_px) > 0) {
      for (i in seq_len(nrow(node_support_df_px))) {
        r <- node_support_df_px[i, ]
        push_svg(
          '<text class="node-support" x="', sprintf("%.2f", r$px - 8),
          '" y="', sprintf("%.2f", r$py - 4), '">',
          svg_escape(r$label),
          '</text>'
        )
      }
    }
    
    for (i in seq_len(nrow(tip_label_df))) {
      r <- tip_label_df[i, ]
      push_svg(
        '<line class="guide-line" x1="', sprintf("%.2f", r$px_tip),
        '" y1="', sprintf("%.2f", r$py),
        '" x2="', sprintf("%.2f", r$px_guide_end),
        '" y2="', sprintf("%.2f", r$py), '"/>'
      )
    }
    
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
    
    for (i in seq_len(nrow(domain_backbone_df_px))) {
      r <- domain_backbone_df_px[i, ]
      push_svg(
        '<line class="backbone" x1="', sprintf("%.2f", r$px),
        '" y1="', sprintf("%.2f", r$py),
        '" x2="', sprintf("%.2f", r$pxend),
        '" y2="', sprintf("%.2f", r$py), '"/>'
      )
    }
    
    DOMAIN_BOX_H <- 10
    
    if (nrow(domain_plot_df_px) > 0) {
      for (i in seq_len(nrow(domain_plot_df_px))) {
        r <- domain_plot_df_px[i, ]
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
    
    for (i in seq_len(nrow(motif_backbone_df_px))) {
      r <- motif_backbone_df_px[i, ]
      push_svg(
        '<line class="backbone" x1="', sprintf("%.2f", r$px),
        '" y1="', sprintf("%.2f", r$py),
        '" x2="', sprintf("%.2f", r$pxend),
        '" y2="', sprintf("%.2f", r$py), '"/>'
      )
    }
    
    MOTIF_BOX_H <- 10
    
    if (nrow(motif_plot_df_px) > 0) {
      for (i in seq_len(nrow(motif_plot_df_px))) {
        r <- motif_plot_df_px[i, ]
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
    
    draw_legend_row <- function(title, items, colors, x_start, y,
                                shape = "rect",
                                max_width = SVG_WIDTH - 40,
                                line_gap = 22) {
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
    
    message("SVG output: ", out_svg)
    
    result_list[[link_step_code]] <- data.frame(
      step_folder = step_name,
      link_step_code = link_step_code,
      status = "success",
      svg = out_svg,
      domain_csv = files$out_domain_csv,
      motif_csv = files$out_motif_csv,
      domain_n = nrow(domain_df),
      motif_n = nrow(motif_df),
      known_n = sum(tip_info$is_known),
      predicted_n = sum(!tip_info$is_known),
      max_seq_len = max(seq_len_df$seq_len, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

message("Domain CSV output: ", files$out_domain_csv)
message("Motif CSV output: ", files$out_motif_csv)

dplyr::bind_rows(result_list)
}
