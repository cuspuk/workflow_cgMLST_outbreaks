import os


def chewbacca_cleanup_names(in_file: str, out_file: str, suffix: str, replacement_map: dict):
    with open(in_file, "r") as f_in:
        with open(out_file, "w") as f_out:
            header = f_in.readline()
            f_out.write(header)
            for line in f_in:
                line = line.strip().split("\t")
                line[0] = os.path.basename(line[0]).replace(suffix, "")
                for old, new in replacement_map.items():
                    line[0] = line[0].replace(old, new)
                f_out.write("\t".join(line) + "\n")


if __name__ == "__main__":
    chewbacca_cleanup_names(
        snakemake.input.tsv, snakemake.output.tsv, snakemake.params.suffix, snakemake.params.replacement_map
    )
