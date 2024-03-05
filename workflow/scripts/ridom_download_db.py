import os
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not os.path.exists(snakemake.output.db_dir):
    os.makedirs(snakemake.output.db_dir, exist_ok=True)

with TemporaryDirectory() as tempdir:
    tmp_zip = os.path.join(tempdir, "test.zip")
    shell("(curl -SL {snakemake.params.url} -o {tmp_zip} &&" " unzip {tmp_zip} -d {snakemake.output.db_dir} ) {log}")
