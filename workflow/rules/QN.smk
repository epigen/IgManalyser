## aggregate
rule quantNorm:
    input: 
        array_file=os.path.join(config["project_path"], "aggregated", "full_samples_combined.csv")
    output:
        array_file_QN=os.path.join(config["project_path"], "normalised_quantiles", "full_samples_normalised.csv"),
        PCA=os.path.join(config["project_path"], "normalised_quantiles", "PCA.RDS")
    params:
    # cluster parameters
        group_by = config['group_by'],
        metadataPath = config['sample_annotation'],
        normPath=os.path.join(config["project_path"], "normalised_quantiles"),
        array_anotation_path = os.path.join(config['array_anotation_path'], "human_array_annotation.tsv"),
        CORTHR=config['CORTHR']
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","quantNorm.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/quantile_normalisation.R"

rule splitByCategory:
    input:
        array_file_QN=os.path.join(config["project_path"], "normalised_quantiles", "full_samples_normalised.csv") if not config.get("use_custom_QN", False) else config['custom_QN_file']
    output:
        array_subsets= [os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', subset, subset + '_samples_subsetted.csv') for subset in split_subsets],
        pca_df = [os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', subset, subset + '_pca.csv') for subset in split_subsets],
        metadata_file = [os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}', subset, subset + '_metadata_subsetted.csv') for subset in split_subsets],
    params:
        metadataPath = config['sample_metadata'],
        split_group = '{split_group}',
        split_subsets = ';'.join(split_subsets), 
        project_dir = os.path.join(config["project_path"], "normalised_quantiles", 'By_{split_group}')
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
        slurm_extra="--qos="+partition
    log:
        os.path.join("logs","rules","quant_split_{split_group}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/split_dataset.R"