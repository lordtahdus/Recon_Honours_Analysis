# ---------- 1) Encoder: adds family, variant, labels, and aesthetics ----------
encode_recon_models <- function(df, model_col = ".model", drop_large_pc = FALSE) {
  df %>%
    mutate(
      method_raw = .data[[model_col]],
      # Families
      family = case_when(
        str_starts(method_raw, "base") ~ "Base",
        str_starts(method_raw, "ols") ~ "OLS",
        str_starts(method_raw, "mint_shr") ~ "MinT-S",
        str_starts(method_raw, "mint_n") ~ "MinT-N",
        TRUE ~ "Other"
      ),
      # Variants
      variant_type = case_when(
        str_detect(method_raw, "_pc_") ~ "pc",
        str_ends(method_raw, "_sv") ~ "sv",
        str_ends(method_raw, "_hcov") ~ "hcov",
        TRUE ~ "vanilla"
      ),
      pc_k = if_else(str_detect(method_raw, "_pc_"),
                     as.integer(str_extract(method_raw, "(?<=_pc_K)\\d+")),
                     NA_integer_),
      # model labels (clean legend text)
      model_label = case_when(
        method_raw == "base" ~ "Base",
        method_raw == "ols" ~ "OLS",
        method_raw == "mint_shr" ~ "MinT-S",
        method_raw == "mint_n" ~ "MinT-N",
        # shrinkage of OLS and Base (for probabilistic forecasts)
        method_raw == "ols_shr" ~ "OLS-S",
        method_raw == "base_shr" ~ "Base-S",
        # PC-adjusted (explicit PC1/PC2; general fallback for other K)
        str_detect(method_raw, "mint_shr_pc_K1") ~ "MinT-S(PC1)",
        str_detect(method_raw, "mint_shr_pc_K2") ~ "MinT-S(PC2)",
        str_detect(method_raw, "mint_n_pc_K1")   ~ "MinT-N(PC1)",
        str_detect(method_raw, "mint_n_pc_K2")   ~ "MinT-N(PC2)",
        str_detect(method_raw, "_pc_K") ~ paste0(if_else(str_detect(method_raw, "mint_shr"), "MinT-S", "MinT-N"),
                                                 "(PC", pc_k, ")"),
        # sv / hcov
        str_ends(method_raw, "_sv") ~ paste0(if_else(str_detect(method_raw, "mint_shr"), "MinT-S", "MinT-N"), "(sv)"),
        str_ends(method_raw, "_hcov") ~ paste0(if_else(str_detect(method_raw, "mint_shr"), "MinT-S", "MinT-N"), "(hcov)"),
        
        
        TRUE ~ method_raw
      ),
      # Shapes for scatter plots:
      # 0 PC (vanilla, sv, hcov) = circle; PC1 = triangle; PC2 = square; other PCs = diamond
      shape_lab = case_when(
        variant_type == "pc" & pc_k == 1 ~ "PC1",
        variant_type == "pc" & pc_k == 2 ~ "PC2",
        variant_type == "pc" & !is.na(pc_k) ~ paste0("PC", pc_k),
        TRUE ~ "0 PC"
      ),
      # Linetypes for line plots
      linetype_lab = case_when(
        family %in% c("Base", "OLS") ~ "solid",
        variant_type == "pc"  ~ "longdash",
        variant_type == "sv"  ~ "dotted",
        variant_type == "hcov" ~ "dotdash",
        TRUE ~ "solid"
      )
    ) %>%
    {
      out <- .
      if (drop_large_pc) {
        out <- out %>% filter(!(variant_type == "pc" & !pc_k %in% c(1, 2)))
      }
      out
    } %>%
    mutate(
      # Factor orders for consistent legends
      family = factor(family, levels = c("Base", "OLS", "MinT-S", "MinT-N")),
      variant_type = factor(variant_type, levels = c("vanilla", "pc", "sv", "hcov")),
      model_label = factor(model_label, levels = c(
        "Base","Base-S",
        "OLS","OLS-S",
        "MinT-S","MinT-S(PC1)","MinT-S(PC2)",
        "MinT-N","MinT-N(PC1)","MinT-N(PC2)",
        "MinT-S(sv)","MinT-N(sv)",
        "MinT-S(hcov)","MinT-N(hcov)"
      ))
    ) |> 
    # turn pc_k to character for better legend handling in ggplot
    mutate(pc_k = ifelse(is.na(pc_k), "0", as.character(pc_k)))
}

# ---------- 2) Reusable ggplot scales ----------
recon_scales <- list(
  # Color by FAMILY (4 colors total)
  scale_color_manual(
    name   = "Family",
    values = c(
      "Base"   = "#7A7A7A",  # grey
      "OLS"    = "#d1d400ff",  # black
      "MinT-S" = "#00c086ff",  # teal (colorblind-friendly)
      "MinT-N" = "#9500c2ff"   # purple
    ),
    drop = FALSE
  ),
  # Linetype by VARIANT
  scale_linetype_manual(
    name   = "Variant",
    values = c(
      "vanilla" = "solid",
      "pc"      = "longdash",
      "sv"      = "dotted",
      "hcov"    = "dotdash"
    ),
    drop = FALSE
  ),
  # Shape for scatter points (0 PC, PC1, PC2; generic PCk shown as diamond)
  scale_shape_manual(
    name   = "PC",
    values = c(
      "0 PC" = NA,  # no point
      "PC1"  = 16,  # triangle
      "PC2"  = 7   # square
      # If you keep K5/K10/K20, add e.g. "PC5"=18 (diamond), etc.
    ),
    drop = FALSE
  )
)

# ---------- 3) Color fill for probabilistic plotting ----------
recon_fill_scale <- scale_fill_manual(
  name = "Family",
  values = c(
    "Base"   = "#9490b7ff",  # grey
    "OLS"    = "#d1d400ff",  # black
    "MinT-S" = "#00c086ff",  # teal (colorblind-friendly)
    "MinT-N" = "#9500c2ff"   # purple
  )
)
