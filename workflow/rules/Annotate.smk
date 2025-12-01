rule CFA_annotation:
    input:
        metadata_file = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_metadata_subsetted.csv'),
    output:
        cfa_df = os.path.join(config["project_path"], "CFA", 'By_{split_group}_{subset}', '{subset}_cfa_annot.csv'),
    params:
        metadataPath = config['sample_metadata'],
        split_group = '{split_group}',
        split_subset = '{subset}', 
        param_columns = config['core_params_to_plot'],
        project_dir = os.path.join(config["project_path"], "CFA", 'By_{split_group}_{subset}')
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","cfa_annot_{split_group}_{subset}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/CFA_annot.R"

rule CFA_PCA:
    input:
        metadata_file = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_metadata_subsetted.csv'),
        pca_file =  os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_pca.csv'),
    output:
        cfa_df = os.path.join(config["project_path"], "CFA", 'By_{split_group}_{subset}', '{subset}_cfa_pca.csv'),
    params:
        metadataPath = config['sample_metadata'],
        split_group = '{split_group}',
        split_subset = '{subset}', 
        param_columns = config['core_params_to_plot'],
        project_dir = os.path.join(config["project_path"], "CFA", 'By_{split_group}_{subset}')
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","cfa_annot_{split_group}_{subset}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/CFA_PCA.R"


rule CorrelateAnnotation:
    input:
        array_subset = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_samples_subsetted.csv'),
        metadata_file = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', '{subset}', '{subset}_metadata_subsetted.csv'),
    output:
        corr_df = os.path.join(config["project_path"], "correlation", 'By_{split_group}_{subset}', '{subset}_corr_df.csv'),
        corr_df_imp = os.path.join(config["project_path"], "correlation", 'By_{split_group}_{subset}', '{subset}_corr_df_imp.csv'),
    params:
        metadataPath = config['sample_metadata'],
        split_group = '{split_group}',
        split_subset = '{subset}', 
        param_columns = config['core_params_to_plot'],
        project_dir = os.path.join(config["project_path"], "correlation", 'By_{split_group}_{subset}')
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
        "../scripts/correlation.R"


rule plotPerPeptide:
    input:
        corr_df_imp = [os.path.join(config["project_path"], "correlation", 'By_{split_group}_'+subset, subset + '_corr_df_imp.csv') for subset in split_subsets],
        mtd_df = [os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', subset, subset + '_metadata_subsetted.csv') for subset in split_subsets],
        array_subsets= [os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', subset, subset + '_samples_subsetted.csv') for subset in split_subsets],
    output:
        corr_dist = os.path.join(config["project_path"], "correlation", 'By_{split_group}_corr_dist.pdf')
    params:
        metadataPath = config['sample_metadata'],
        split_group = '{split_group}',
        split_subsets = ','.join(split_subsets), 
        param_columns = config['core_params_to_plot'],
        project_dir = os.path.join(config["project_path"], "correlation")
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","perPeptide_{split_group}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/correlation_per_peptide.R"