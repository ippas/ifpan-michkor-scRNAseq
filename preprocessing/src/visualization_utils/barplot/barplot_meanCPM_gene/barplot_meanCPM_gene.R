barplot_meanCPM_gene <- function(
    summary_list,
    cluster_name,
    gene_symbol,
    bar_color = "black",
    bar_width = 0.7,
    bar_size = 1,
    add_points = TRUE,
    add_se = FALSE,
    se_width = 0.15,
    point_color = "black",
    point_size = 10,
    log2_transform = FALSE,
    y_limits = NULL,
    show_p = TRUE,
    text_size_values = 12,
    text_size_p      = 12,
    text_size_title  = 18,
    text_size_axis   = 12,
    xaxis_text_size  = 14,
    yaxis_text_size  = 14,
    custom_title     = NULL,
    stats_bracket_mult = 0.75,
    stats_text_mult    = 0.28
) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # convert sizes from pt to mm
  text_size_values_mm <- text_size_values / 2.83
  text_size_p_mm      <- text_size_p      / 2.83
  point_size_mm       <- point_size       / 2.83
  
  # extract row for gene
  df_raw <- summary_list[[cluster_name]] %>%
    filter(gene_symbol == !!gene_symbol)
  
  # prepare summary df for bars
  df <- df_raw %>%
    select(meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX) %>%
    pivot_longer(cols = everything(),
                 names_to = "condition",
                 values_to = "mean_CPM") %>%
    mutate(
      condition = recode(condition,
                         "meanCPM_CTRL" = "CTRL",
                         "meanCPM_CORT" = "CORT",
                         "meanCPM_DEX"  = "DEX"),
      condition = factor(condition, levels=c("CTRL","CORT","DEX"))
    )
  
  # prepare points for jitter
  df_points <- df_raw %>%
    select(starts_with("A")) %>%
    pivot_longer(cols = everything(), names_to="sample", values_to="expr") %>%
    mutate(
      condition = case_when(
        sample %in% paste0("A",1:4)  ~ "CTRL",
        sample %in% paste0("A",5:8)  ~ "DEX",
        sample %in% paste0("A",9:12) ~ "CORT",
        TRUE ~ NA_character_
      ),
      condition = factor(condition, levels=c("CTRL","CORT","DEX"))
    )
  
  # optional log2 transform
  if (log2_transform) {
    df        <- df        %>% mutate(mean_CPM = log2(mean_CPM + 1))
    df_points <- df_points %>% mutate(expr     = log2(expr     + 1))
  }
  
  # compute SE if needed
  if (add_se) {
    df_se <- df_points %>%
      group_by(condition) %>%
      summarise(
        SE = sd(expr, na.rm=TRUE) / sqrt(n()),
        mean_CPM = mean(expr, na.rm=TRUE),
        .groups="drop"
      )
  }
  
  # determine ymax based on points if available
  if (add_points) {
    ymax <- max(df_points$expr, df$mean_CPM, na.rm = TRUE)
  } else {
    ymax <- max(df$mean_CPM, na.rm = TRUE)
  }
  
  pad <- ymax * 0.18
  
  # default title
  auto_title  <- paste0(gene_symbol," mean CPM in ", cluster_name,
                        ifelse(log2_transform," (log2)",""))
  final_title <- ifelse(is.null(custom_title), auto_title, custom_title)
  
  # base plot
  p <- ggplot(df, aes(x=condition, y=mean_CPM)) +
    geom_col(color=bar_color, fill=NA, width=bar_width, linewidth=bar_size)
  
  # errorbars
  if (add_se) {
    p <- p + geom_errorbar(data=df_se,
                           aes(ymin=mean_CPM-SE, ymax=mean_CPM+SE),
                           width=se_width, linewidth=.8)
  }
  
  # text above bars
  p <- p + geom_text(aes(label=round(mean_CPM,2)),
                     vjust=-2, size=text_size_values_mm)
  
  # points
  if (add_points) {
    p <- p + geom_point(data=df_points, aes(y=expr),
                        position=position_jitter(width=0.1),
                        size=point_size_mm, color=point_color)
  }
  
  # theme
  p <- p +
    theme_classic(base_size=text_size_axis) +
    scale_x_discrete(limits=c("CTRL","CORT","DEX")) +
    labs(title=final_title,
         x="", y=ifelse(log2_transform,"log2(Mean CPM+1)","Mean CPM")) +
    theme(
      plot.title = element_text(size=text_size_title, face="plain", hjust=.5),
      axis.text.x = element_text(size=xaxis_text_size, face="plain"),
      axis.text.y = element_text(size=yaxis_text_size, face="plain")
    )
  
  # p-value formatting helper
  format_pval <- function(p) {
    if (is.na(p)) return("NA")
    if (p >= 0.0001) return(format(round(p,4), nsmall=4))
    return(format(p, scientific=TRUE, digits=3))
  }
  
  # add statistics brackets
  if (show_p) {
    pval_cort <- format_pval(as.numeric(df_raw$PValue_edgerCORT))
    pval_dex  <- format_pval(as.numeric(df_raw$PValue_edgerDEX))
    
    label_cort <- paste0("p = ", pval_cort,", log2FC = ",
                         round(df_raw$log2FC_CORT_CPM,2))
    label_dex <- paste0("p = ", pval_dex,", log2FC = ",
                        round(df_raw$log2FC_DEX_CPM,2))
    
    top1 <- ymax + pad * stats_bracket_mult
    top2 <- ymax + pad * (stats_bracket_mult * 2.0)
    tick <- pad * .10
    
    p <- p +
      geom_segment(aes(x=1,xend=2,y=top1,yend=top1), linewidth=.9) +
      geom_segment(aes(x=1,xend=1,y=top1,yend=top1-tick), linewidth=.9) +
      geom_segment(aes(x=2,xend=2,y=top1,yend=top1-tick), linewidth=.9) +
      annotate("text", x=1.5, y=top1 + pad*stats_text_mult,
               label=label_cort, size=text_size_p_mm, lineheight=.95) +
      
      geom_segment(aes(x=1,xend=3,y=top2,yend=top2), linewidth=.9) +
      geom_segment(aes(x=1,xend=1,y=top2,yend=top2-tick), linewidth=.9) +
      geom_segment(aes(x=3,xend=3,y=top2,yend=top2-tick), linewidth=.9) +
      annotate("text", x=2, y=top2 + pad*stats_text_mult,
               label=label_dex, size=text_size_p_mm, lineheight=.95)
  }
  
  return(p)
}
