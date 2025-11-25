###quantile normalisation per group for two groups (ideally to do it for more!)
get_QN_norm <- function(df, ids_group1, ids_group2){
    mtx_chow <- as.matrix(df[, ids_group1])
    norm_chow = normalize.quantiles(as.matrix(mtx_chow))
    colnames(norm_chow) <- colnames(mtx_chow)
    mtx_hcd <- as.matrix(df[, ids_group2])
    norm_hcd = normalize.quantiles(as.matrix(mtx_hcd))
    colnames(norm_hcd) <- colnames(mtx_hcd)

    array_norm <- merge(as.data.frame(norm_hcd, row.names = row.names(df)),
           as.data.frame(norm_chow, row.names = row.names(df)),
          by = 'row.names', all = TRUE)
    row.names(array_norm) <- as.character(array_norm$Row.names)
    array_norm <- array_norm[, -c(1)]
    return(array_norm)
}

##QN if we only have one group
get_QN_norm_onegroup <- function(df, ids_group1){
    mtx_chow <- as.matrix(df[, ids_group1])
    norm_chow = normalize.quantiles(as.matrix(mtx_chow))
    colnames(norm_chow) <- colnames(mtx_chow)
    array_norm <- as.data.frame(norm_chow, row.names= row.names(df))
    return(array_norm)
}


### run PCA
run_pca <- function(df, outdir){
    full_pca <- prcomp(df, scale = FALSE)

    var_expl <- data.table(sdev = full_pca$sdev, comp = paste0("PC", seq(1, length(full_pca$sdev))))
    var_expl$comp <- factor(var_expl$comp, levels = var_expl$comp)
    
    ggplot(var_expl, aes(x = comp, y = sdev)) + geom_bar(stat = "identity")
    ggsave(file.path(outdir,  "var_expl.pdf"), width = 5, height = 5)
    df_rot <- as.data.frame(full_pca$rotation)
    saveRDS(full_pca,file.path(outdir,  "PCA.RDS"))
    return(df_rot)
    }