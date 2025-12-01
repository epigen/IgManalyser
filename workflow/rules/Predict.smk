rule VolcanoPlots:
    input:
        array_subset = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_samples_subsetted.csv'),
        annot_df = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_metadata_subsetted.csv'),
    output:
        volcano_df = os.path.join(config["project_path"], "volcano", 'By_{split_group}_{subset}', "Gensini", "logFC_df.csv")
    params:
        metadata = config['sample_metadata'],
        gensini_col = config['gensini_score'],
        gensini_thr = '6,32',
        project_dir = os.path.join(config["project_path"], "volcano", 'By_{split_group}_{subset}', "Gensini")
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","quant_split_{split_group}_{subset}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/volcano_plots.R"


rule VizTopPeptides_Gensini:
    input: 
        array_subset = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_samples_subsetted.csv'),
        annot_df = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}',  '{subset}_metadata_subsetted.csv'),
        volcano_df = os.path.join(config["project_path"], "volcano", 'By_{split_group}_{subset}', "Gensini", "logFC_df.csv"),
        corr_df = os.path.join(config["project_path"], "correlation", 'By_{split_group}_{subset}', '{subset}_corr_df_imp.csv'),
    output:
        pept_hits = os.path.join(config["project_path"], "top_hit_viz", 'By_{split_group}_{subset}', "petides_hits.csv")
    params:
        metadata = config['sample_metadata'],
        gensini_col = config['gensini_score'],
        gensini_thr = '6,32',
        CORR_THR = 0.45,
        LFC_THR = 0.2,
        project_dir = os.path.join(config["project_path"], "top_hit_viz", 'By_{split_group}_{subset}'),
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","quant_split_{split_group}_{subset}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/peptide_plots.R"
