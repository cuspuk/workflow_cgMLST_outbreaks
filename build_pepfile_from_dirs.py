import argparse
import os
import sys
from dataclasses import dataclass

DEFAULT_RELATIVE_TAXA_PATH = "results/taxonomy/{sample}/parsed_taxa.txt"
DEFAULT_RELATIVE_FASTA_PATH = "results/assembly/{sample}/assembly.fasta"
DEFAULT_PEPFILE = "config/pep/samples.csv"


class SampleWithoutOutputException(Exception):
    sample_name: str
    analysis_dir: str
    files: list[str]

    def __init__(self, sample_name: str, analysis_dir: str, files: list[str]) -> None:
        self.sample_name = sample_name
        self.analysis_dir = analysis_dir
        message = f"{sample_name=} was found to be processed in {analysis_dir=} but there were missing {files=}"
        super().__init__(message)


@dataclass
class Sample:
    sample_name: str
    fasta: str
    GTDBtk_taxa: str

    def __str__(self):
        return f"{self.sample_name},{self.fasta},{self.GTDBtk_taxa}"

    @classmethod
    def header(cls):
        return "sample_name,fasta,GTDBtk_taxa"


def parse_pepfile(root: str) -> list[str]:
    pepfile = os.path.join(root, DEFAULT_PEPFILE)
    with open(pepfile, "r") as f:
        header = f.readline().strip().split(",")
        if "sample_name" in header:
            sample_idx = header.index("sample_name")
        elif "sample" in header:
            sample_idx = header.index("sample")
        else:
            raise ValueError(f"Could not find sample_name in the header of {pepfile}")
        samples = [x.rstrip().split(",")[sample_idx] for x in f.readlines() if len(x) > 1]
    return samples


def parse_sample(root_dir: str, sample: str) -> Sample:
    taxa_file = os.path.join(root_dir, DEFAULT_RELATIVE_TAXA_PATH.format(sample=sample))
    assembly_file = os.path.join(root_dir, DEFAULT_RELATIVE_FASTA_PATH.format(sample=sample))
    if not os.path.exists(taxa_file) or not os.path.exists(assembly_file):
        if os.path.exists(assembly_file):
            raise SampleWithoutOutputException(sample, root_dir, [taxa_file])
        else:
            raise SampleWithoutOutputException(sample, root_dir, [assembly_file, taxa_file])
    with open(taxa_file, "r") as f:
        taxa = f.readline().strip()

    return Sample(sample, assembly_file, taxa)


def parse_analysis_dir(root_dir: str):
    sample_names = parse_pepfile(root_dir)

    success_samples: list[Sample] = []
    failed_samples: list[SampleWithoutOutputException] = []
    for name in sample_names:
        try:
            success_samples.append(parse_sample(root_dir, name))
        except SampleWithoutOutputException as err:
            failed_samples.append(err)
    return success_samples, failed_samples


def process_analysis_dirs(analysis_dirs: list[str]):
    success_samples_all: dict[str, list[Sample]] = {}
    failed_samples_all: dict[str, list[SampleWithoutOutputException]] = {}
    for analysis_dir in analysis_dirs:
        success_samples, failed_samples = parse_analysis_dir(analysis_dir)
        success_samples_all[analysis_dir] = success_samples
        failed_samples_all[analysis_dir] = failed_samples

    print(Sample.header())
    for analysis_dir, success_samples in success_samples_all.items():
        for sample in success_samples:
            print(sample)

    print("-" * 80, file=sys.stderr)
    print("Error log:", file=sys.stderr)
    for analysis_dir, failed_samples in failed_samples_all.items():
        for sample in failed_samples:
            print(sample, file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a pepfile from the results of the analysis")
    parser.add_argument(
        "analysis_dir",
        nargs="+",
        help="Analysis directory having results/ directory in it",
    )
    args = parser.parse_args()

    process_analysis_dirs(args.analysis_dir)
