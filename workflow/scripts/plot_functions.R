### set colors
risk_colors = c("1" = "#FD8008", "0" = "#009193", 
               "Asymptomatic" = "darkblue", "Symptomatic" = "darkred")

### plot peptide boxplots

plot_peptides_boxplots <- function(df, peptide_list, outpath, save_index, gensini_col){
    print(paste0("plotting ", length(peptide_list), " peptides"))
    g <- ggplot(df[df$coordinate_unique %in% peptide_list, ], 
       aes(x = sequence_label, y = value, color = as.factor(.data[[gensini_col]]))) + 
    geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge()) + 
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)) + scale_x_discrete(labels = scales::label_wrap(15)) + 
    scale_color_manual(values = risk_colors) + labs(color = "Risk", x = "", y = "peptide intensity (QN-norm.)")  + 
    stat_compare_means(aes(group = as.factor(.data[[gensini_col]]), 
                            method = 't.test',
                            label = sprintf("%5.2f", as.numeric(..p.format..)))
                      )

    ggsave(file.path(outpath, paste0(save_index, ".pdf")), g, 
           width = min(30, length(peptide_list) + 2), height = 8)  
    
    g1 <- ggplot(df[df$coordinate_unique %in% peptide_list, ], 
       aes(x = sequence_label_short, y = value, color = as.factor(.data[[gensini_col]]))) + 
    geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge()) + 
    theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)) + 
    scale_color_manual(values = risk_colors) + 
    labs(color = "Risk", x = "", y = "peptide intensity (QN-norm.)")  + 
    stat_compare_means(aes(group = as.factor(.data[[gensini_col]]),
                           method = 't.test',
                          label = sprintf("%5.2f", as.numeric(..p.format..))))
    ggsave(file.path(outdir, paste0(save_index, "_sequence.pdf")), g1, 
           width = max(3, min(30, length(peptide_list))), height = 8)    
}