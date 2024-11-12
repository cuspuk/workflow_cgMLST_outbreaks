rule preprocess_assembly_for_chewbacca:
    input:
        assembly=infer_assembly_fasta,
    output:
        assembly=temp("results/assembly/cleaned/{sample}_cleaned.fasta"),
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
        file_lst=lambda wildcards, input: "\n".join([os.path.abspath(x) for x in input.assemblies]),
    localrule: True
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/prepare_assemblies_paths/{taxa_label}.log",
    shell:
        "echo -e {params.file_lst:q} > {output} 2> {log}"


rule dump_samples:
    output:
        "results/_metadata/{taxa_label}.ini",
    params:
        samples=lambda wildcards: "\n".join(get_sample_names_for_taxa_label(wildcards.taxa_label)),
    localrule: True
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/dump_samples/{taxa_label}.log",
    shell:
        "echo -e {params.samples:q} > {output} 2> {log}"
