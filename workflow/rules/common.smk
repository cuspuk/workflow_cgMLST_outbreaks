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
    taxa_label: str
    GTDBtk_taxa: list[str]
    cgMLST_schema_dir: str
    cgMLST_schema_download_url: str | None
    training_file_full_path: str
    training_file_download_url: str | None

    def __post_init__(self):
        if not self.cgMLST_schema_download_url:
            if not os.path.exists(self.cgMLST_schema_dir):
                raise ValueError(
                    f"cgMLST schema directory {self.cgMLST_schema_dir} does not exist and no download URL was provided."
                )
        if not self.training_file_download_url:
            if not os.path.exists(self.training_file_full_path):
                raise ValueError(
                    f"Training file {self.training_file_full_path} does not exist and no download URL was provided."
                )
        if config["chewbacca_schemas_dir"].startswith(self.cgMLST_schema_dir):
            raise ValueError(
                f"Config invalid. cgMLST schema directory for taxa_label={self.taxa_label} is a subdirectory of chewbacca_schemas_dir."
            )
        if self.cgMLST_schema_dir.startswith(config["chewbacca_schemas_dir"]):
            raise ValueError(
                f"Config invalid. chewbacca_schemas_dir is a subdirectory of cgMLST_schema_dir for taxa_label={self.taxa_label}."
            )


MAPPINGS: list[SchemaMapping] = [SchemaMapping(**cfg) for cfg in config["organism_schemas_mapping"]]


def get_constraints():
    return {
        "taxa_label": "|".join([mapping.taxa_label for mapping in MAPPINGS]),
        "cgMLST_schema_dir": "|".join(
            [mapping.cgMLST_schema_dir for mapping in MAPPINGS if mapping.cgMLST_schema_download_url]
        ),
        "training_file_full_path": "|".join(
            [mapping.training_file_full_path for mapping in MAPPINGS if mapping.training_file_download_url]
        ),
    }


def get_fasta_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["fasta"]][0]


def get_mapping_for_taxa_label(taxa_label: str) -> SchemaMapping:
    return next(mapping for mapping in MAPPINGS if mapping.taxa_label == taxa_label)


def get_sample_names_for_taxa_label(taxa_label: str) -> list[str]:
    mapping = get_mapping_for_taxa_label(taxa_label)
    return list(pep.sample_table[pep.sample_table["GTDBtk_taxa"].isin(mapping.GTDBtk_taxa)]["sample_name"].values)


def get_valid_taxa_labels():
    return [mapping.taxa_label for mapping in MAPPINGS if get_sample_names_for_taxa_label(mapping.taxa_label)]


### Global rule-set stuff #############################################################################################


def infer_assembly_fasta(wildcards) -> str:
    return get_fasta_for_sample_from_pep(wildcards.sample)


def infer_chewbacca_schema_for_taxa_label(wildcards) -> str:
    return os.path.join(config["chewbacca_schemas_dir"], wildcards.taxa_label)


def infer_training_file_for_taxa_label(wildcards) -> str:
    return get_mapping_for_taxa_label(wildcards.taxa_label).training_file_full_path


def infer_cgMLST_schema_dir_for_taxa_label(wildcards) -> str:
    return get_mapping_for_taxa_label(wildcards.taxa_label).cgMLST_schema_dir


def infer_cleaned_assemblies_for_taxa_label(wildcards):
    return expand(
        "results/assembly/cleaned/{sample}_cleaned.fasta", sample=get_sample_names_for_taxa_label(wildcards.taxa_label)
    )


def infer_url_for_training_file(wildcards):
    for mapping in MAPPINGS:
        if mapping.training_file_full_path != wildcards.training_file_full_path:
            continue

        if not mapping.training_file_download_url or mapping.training_file_download_url == "null":
            return ""
        return mapping.training_file_download_url


def infer_url_for_schema_download(wildcards):
    for mapping in MAPPINGS:
        if mapping.cgMLST_schema_dir != wildcards.cgMLST_schema_dir:
            continue

        if not mapping.cgMLST_schema_download_url or mapping.cgMLST_schema_download_url == "null":
            return ""
        return mapping.cgMLST_schema_download_url


def get_outputs():
    taxa_labels = get_valid_taxa_labels()
    if config.get("without_cgmlst_dist", False):
        outs = {
            "cgmlst": expand("results/cgMLST/{taxa_label}/extracted_genes/cgMLST95.tsv", taxa_label=taxa_labels),
        }
    else:
        taxa_labels_with_enough_samples = []
        for label in taxa_labels:
            if len(get_sample_names_for_taxa_label(label)) > 2:
                taxa_labels_with_enough_samples.append(label)
        outs = {
            "distances": expand(
                "results/cgMLST/{taxa_label}/extracted_genes/cgMLST95_tree.newick",
                taxa_label=taxa_labels_with_enough_samples,
            ),
        }
    samples = {"samples": expand("results/_metadata/{taxa_label}.ini", taxa_label=taxa_labels)}
    return {**outs, **samples}


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


### Resource handling #################################################################################################


def get_mem_mb_for_XY(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["XY_mem_mb"] * attempt)
