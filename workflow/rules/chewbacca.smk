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
        "(rm -rf {params.out_dir} && chewBBACA.py AlleleCall --input-files {input.assembly_file} --schema-directory {input.schema_dir}"
        " --output-directory {params.out_dir} --cpu {threads}) > {log} 2>&1"


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
