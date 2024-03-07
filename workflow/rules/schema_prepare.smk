rule ridom_schema_download:
    output:
        db_dir=directory("{schemas_root_dir}/ridom/{ridom_schema_name}"),
    params:
        url=lambda wildcards, output: f"https://www.cgmlst.org/ncs/schema/{wildcards.ridom_schema_name}/alleles/",
    localrule: True
    conda:
        "../envs/python_downloader.yaml"
    log:
        "{schemas_root_dir}/logs/ridom/{ridom_schema_name}.log",
    script:
        "../scripts/ridom_download_db.py"


rule prodigal_download_training_file:
    output:
        "{schemas_root_dir}/chewbacca/{schema_name}/{training_file}",
    params:
        url=lambda wildcards: f"https://github.com/B-UMMI/chewBBACA/raw/master/CHEWBBACA/prodigal_training_files/{wildcards.training_file}",
    localrule: True
    conda:
        "../envs/curl.yaml"
    log:
        "{schemas_root_dir}/logs/chewbacca/download_training_file/{schema_name}_{training_file}.log",
    shell:
        "curl -SL {params.url} -o {output} > {log} 2>&1"


rule chewbacca_prepare_schema:
    input:
        db_dir="{schemas_root_dir}/{schema_name}",
        trn=infer_training_file_for_taxa_label,
    output:
        summary="{schemas_root_dir}/chewbacca/{schema_name}/{taxa_label}_summary_stats.tsv",
        out_dir=directory("{schemas_root_dir}/chewbacca/{schema_name}/{taxa_label}"),
    conda:
        "../envs/chewbacca.yaml"
    threads: min(config["threads"]["chewbacca"], config["max_threads"])
    log:
        "{schemas_root_dir}/logs/chewbacca/prepare/{schema_name}_{taxa_label}.log",
    shell:
        "chewBBACA.py PrepExternalSchema --schema-directory {input.db_dir} --output-directory {output.out_dir}"
        " --ptf {input.trn} --cpu {threads} > {log} 2>&1"
