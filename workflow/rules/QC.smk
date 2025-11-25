## perform per-sample QC
rule QC:
    input:
        get_raw_bams,
    output:
        QC_file=os.path.join(results_dir, "{sample}", "replica_compared.csv"),
        corr_file=os.path.join(results_dir, "{sample}", "corr_summary.csv")
    params:
    # cluster parameters
        QC_path=os.path.join(results_dir, "{sample}"),
        sample_id="{sample}",
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
    log:
        os.path.join("logs","rules","QC_{sample}.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/array_preprocessing.R"
        
        
## aggregate
rule aggregate:
    input: 
        array_files=expand(os.path.join(results_dir, "{sample_name}", "replica_compared.csv"), sample_name=samples.keys()),
        corr_files=expand(os.path.join(results_dir, "{sample_name}", "corr_summary.csv"), sample_name=samples.keys())
    output:
        array_file=os.path.join(config["project_path"], "aggregated", "full_samples_combined.csv")
    params:
    # cluster parameters
        group_by = config['group_by'],
        metadataPath = config['sample_annotation'],
        aggregatePath=os.path.join(config["project_path"], "aggregated"),
        array_anotation_path = config['array_anotation_path'],
        CORTHR=config['CORTHR']
    threads: threads
    resources:
        mem=mem,
        slurm_partition=partition,
    log:
        os.path.join("logs","rules","aggregate.log"),
    conda:
        "../../envs/peptide_array_r.yaml",
    script:
        "../scripts/aggregate.R"