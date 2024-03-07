from snakemake.utils import validate
from dataclasses import dataclass


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


pepfile: config.get("pepfile", "config/pep/config.yaml")


validate(pep.sample_table, "../schemas/samples.schema.yaml")


### Layer for adapting other workflows  ###############################################################################


### Data input handling independent of wildcards ######################################################################


@dataclass
class SchemaMapping:
    GTDBtk_taxa: list[str]
    cgMLST_ridom_schema: str
    training_file: str
    taxa_label: str


MAPPINGS: list[SchemaMapping] = [SchemaMapping(**cfg) for cfg in config["organism_schemas_mapping"]]

SCHEMA_DIR = os.path.abspath(config["schemas_dir"])


def get_constraints():
    return {
        "schemas_dir": SCHEMA_DIR,
        "schema_name": "|".join([mapping.cgMLST_ridom_schema for mapping in MAPPINGS]),
        "training_file": "|".join([mapping.training_file for mapping in MAPPINGS]),
        "taxa_label": "|".join([mapping.taxa_label for mapping in MAPPINGS]),
    }


# def get_sample_names() -> list[str]:
#     return list(pep.sample_table["sample_name"].values)


def get_fasta_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["fasta"]][0]


def get_mapping_for_taxa_label(taxa_label: str) -> SchemaMapping:
    return next(mapping for mapping in MAPPINGS if mapping.taxa_label == taxa_label)


def get_sample_names_for_taxa_label(taxa_label: str) -> list[str]:
    mapping = get_mapping_for_taxa_label(taxa_label)
    return list(pep.sample_table[pep.sample_table["GTDBtk_taxa"].isin(mapping.GTDBtk_taxa)]["sample_name"].values)


def get_taxa_labels():
    return [mapping.taxa_label for mapping in MAPPINGS]


# def get_all_assemblies() -> list[str]:
#     return list(pep.sample_table["fasta"].values)


### Global rule-set stuff #############################################################################################


def infer_assembly_fasta(wildcards) -> str:
    return get_fasta_for_sample_from_pep(wildcards.sample)


def infer_chewbacca_schema_for_taxa_label(wildcards) -> str:
    schema_name = get_mapping_for_taxa_label(wildcards.taxa_label).cgMLST_ridom_schema
    return os.path.join(SCHEMA_DIR, "chewbacca", schema_name, wildcards.taxa_label)


def infer_training_file_for_taxa_label(wildcards) -> str:
    mapping = get_mapping_for_taxa_label(wildcards.taxa_label)
    return os.path.join(SCHEMA_DIR, "chewbacca", mapping.cgMLST_ridom_schema, mapping.training_file)


def infer_cleaned_assemblies_for_taxa_label(wildcards):
    return expand(
        "results/assembly/cleaned/{sample}.fasta", sample=get_sample_names_for_taxa_label(wildcards.taxa_label)
    )


def get_outputs():
    return {
        "distances": expand(
            "results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_distances.tsv", taxa_label=get_taxa_labels()
        ),
    }


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


### Resource handling #################################################################################################


def get_mem_mb_for_XY(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["XY_mem_mb"] * attempt)
