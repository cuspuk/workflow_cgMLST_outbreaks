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
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95.tsv",
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.html),
    conda:
        "../envs/chewbacca.yaml"
    log:
        "logs/chewbacca/extract_cgMLST/{taxa_label}.log",
    shell:
        "chewBBACA.py ExtractCgMLST -i {input.tsv} -o {params.out_dir} > {log} 2>&1"


rule chewbacca_cleanup_names:
    input:
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95.tsv",
    output:
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_cleaned.tsv",
    params:
        suffix="_cleaned",
    localrule: True
    conda:
        "../envs/python.yaml"
    log:
        "logs/cgMLST_distances/{taxa_label}_cleanup.log",
    script:
        "../scripts/chewbacca_cleanup_names.py"


rule cgMLST_distances:
    input:
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_cleaned.tsv",
    output:
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_distances.tsv",
    conda:
        "../envs/cgmlst_dists.yaml"
    log:
        "logs/cgMLST_distances/{taxa_label}.log",
    shell:
        "cgmlst-dists {input.tsv} > {output.tsv} 2> {log}"


rule newick_tree:
    input:
        tsv="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_distances.tsv",
    output:
        newick="results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_tree.newick",
    conda:
        "../envs/newick.yaml"
    log:
        "logs/cgMLST_distances/{taxa_label}_tree.log",
    script:
        "../scripts/newick.R"
