
rule ridom_schema_download:
    output:
        db_dir=directory("{schemas_dir}/ridom/{schema_name}"),
    params:
        url=lambda wildcards, output: f"https://www.cgmlst.org/ncs/schema/{wildcards.schema_name}/alleles/",
    localrule: True
    conda:
        "../envs/python_downloader.yaml"
    log:
        "{schemas_dir}/logs/ridom/{schema_name}.log",
    script:
        "../scripts/ridom_download_db.py"


rule prodigal_download_training_file:
    output:
        "{schemas_dir}/chewbacca/{schema_name}/{training_file}",
    params:
        url=lambda wildcards: f"https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/prodigal_training_files/{wildcards.training_file}",
    localrule: True
    conda:
        "../envs/curl.yaml"
    log:
        "{schemas_dir}/logs/chewbacca/download_training_file/{schema_name}_{training_file}.log",
    shell:
        "curl -SL {params.url} -o {output} > {log} 2>&1"


rule chewbacca_prepare_schema:
    input:
        db_dir="{schemas_dir}/ridom/{schema_name}",
        trn=infer_training_file_for_taxa_label,
    output:
        summary="{schemas_dir}/chewbacca/{schema_name}/{taxa_label}_summary_stats.tsv",
        out_dir=directory("{schemas_dir}/chewbacca/{schema_name}/{taxa_label}"),
    conda:
        "../envs/chewbacca.yaml"
    threads: min(config["threads"]["chewbacca"], config["max_threads"])
    log:
        "{schemas_dir}/logs/chewbacca/prepare/{schema_name}_{taxa_label}.log",
    shell:
        "chewBBACA.py PrepExternalSchema --schema-directory {input.db_dir} --output-directory {output.out_dir}"
        " --ptf {input.trn} --cpu {threads} > {log} 2>&1"


rule preprocess_assembly_for_chewbacca:
    input:
        assembly=infer_assembly_fasta,
    output:
        assembly=temp("results/assembly/cleaned/{sample}.fasta"),
    conda:
        "../envs/sed.yaml"
    localrule: True
    log:
        "logs/preprocess_assembly_for_chewbacca/{sample}.log",
    shell:
        "sed 's/=/_/g' {input:q} > {output:q} 2> {log}"


rule prepare_assemblies_paths:
    input:
        assemblies=infer_cleaned_assemblies_for_taxa_label,
    output:
        assembly="results/assembly/{taxa_label}.infile",
    params:
        file_lst=lambda wildcards, input: "\n".join([str(x) for x in input.assemblies]),
    localrule: True
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/prepare_assemblies_paths/{taxa_label}.log",
    shell:
        "echo -e {params.file_lst} > {output} 2> {log}"


rule chewbacca_allele_call:
    input:
        schema_dir=infer_chewbacca_schema_for_taxa_label,
        assembly_file="results/assembly/{taxa_label}.infile",
        assemblies=infer_cleaned_assemblies_for_taxa_label,
    output:
        multiext(
            "results/cgMLST/{taxa_label}/",
            "cds_coordinates.tsv",
            "invalid_cds.txt",
            "loci_summary_stats.tsv",
            "results_statistics.tsv",
            "results_contigsInfo.tsv",
            "results_alleles.tsv",
            "paralogous_counts.tsv",
            "paralogous_loci.tsv",
        ),
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/chewbacca.yaml"
    threads: min(config["threads"]["chewbacca"], config["max_threads"])
    log:
        "logs/chewbacca/allele_call/{taxa_label}.log",
    shell:
        "chewBBACA.py AlleleCall --input-files {input.assembly_file} --schema-directory {input.schema_dir}"
        " --output-directory {params.out_dir} --cpu {threads} > {log} 2>&1"


rule chewbacca_remove_genes:
    input:
        tsv="results/cgMLST/{taxa_label}/results_alleles.tsv",
        paralogs="results/cgMLST/{taxa_label}/paralogous_counts.tsv",
    output:
        tsv="results/cgMLST/{taxa_label}/results_alleles_NoParalogs.tsv",
    conda:
        "../envs/chewbacca.yaml"
    log:
        "logs/chewbacca/remove_genes/{taxa_label}.log",
    shell:
        "chewBBACA.py RemoveGenes -i {input.tsv} -g {input.paralogs} -o {output.tsv} > {log} 2>&1"


rule chewbacca_extract_cgMLST:
    input:
        tsv="results/cgMLST/{taxa_label}/results_alleles_NoParalogs.tsv",
    output:
        html="results/cgMLST/{taxa_label}/extracted_genes/cgMLST.html",
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.html),
    conda:
        "../envs/chewbacca.yaml"
    log:
        "logs/chewbacca/extract_cgMLST/{taxa_label}.log",
    shell:
        "chewBBACA.py ExtractCgMLST -i {input.tsv} -o {params.out_dir} > {log} 2>&1"
