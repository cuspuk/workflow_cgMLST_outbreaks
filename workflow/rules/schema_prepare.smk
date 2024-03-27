rule ridom_schema_download:
    output:
        db_dir=directory("{cgMLST_schema_dir}"),
    params:
        url=infer_url_for_schema_download,
    localrule: True
    conda:
        "../envs/python_downloader.yaml"
    log:
        "{cgMLST_schema_dir}/download_schema.log",
    script:
        "../scripts/ridom_download_db.py"


rule prodigal_download_training_file:
    output:
        "{training_file_full_path}",
    params:
        url=infer_url_for_training_file,
    localrule: True
    conda:
        "../envs/curl.yaml"
    log:
        "{training_file_full_path}.download.log",
    shell:
        "curl -SL {params.url} -o {output} > {log} 2>&1"


rule chewbacca_prepare_schema:
    input:
        db_dir=infer_cgMLST_schema_dir_for_taxa_label,
        trn=infer_training_file_for_taxa_label,
    output:
        summary="{chewbacca_schemas_dir}/{taxa_label}_summary_stats.tsv",
        out_dir=directory("{chewbacca_schemas_dir}/{taxa_label}"),
    conda:
        "../envs/chewbacca.yaml"
    threads: min(config["threads"]["chewbacca"], config["max_threads"])
    log:
        "{chewbacca_schemas_dir}/logs/{taxa_label}.log",
    shell:
        "chewBBACA.py PrepExternalSchema --schema-directory {input.db_dir} --output-directory {output.out_dir}"
        " --ptf {input.trn} --cpu {threads} > {log} 2>&1"
