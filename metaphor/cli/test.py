#!/usr/bin/env python

__doc__ = """
    Tests Metaphor.

    1. Checks for test data.
        1a. Downloads if doesn't match the checksum.
    2. Create input file.
    3. Run workflow with default config.
    4. Generate report.
    5. Delete output.
    """


from argparse import Namespace
from pathlib import Path
from hashlib import md5
from textwrap import dedent
from subprocess import check_call, CalledProcessError

import requests
from tqdm import tqdm

from metaphor import wrapper_prefix
from metaphor.workflow import snakefile
from metaphor.config import test_config
from metaphor.utils import confirm_message, get_successful_completion

from .create_input_table import main as create_input_table


def get_md5(file):
    with open(file, "rb") as f:
        filehash = md5(f.read())

    return filehash.hexdigest()


def download_file(url, filename):
    """
    Downloads a large file from a URL.
    From: https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    """

    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        # print(f"\tWriting to '{filename}'.")
        with open(filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                f.write(chunk)
    return filename


def check_test_directory(directory):
    directory = Path(directory)
    if not (directory.exists() and directory.is_dir()):
        print(f"Creating test directory: '{directory}'.")
        directory.mkdir(exist_ok=True)
    else:
        print("Directory exists, skipping.")
    return directory


def download_data(directory, test_files, download_url):
    test_directory = check_test_directory(directory)
    downloaded, skipped = 0, 0
    test_files = {
        test_directory.joinpath(filename): hexdigest
        for filename, hexdigest in test_files.items()
    }

    # Skip file download if they already exist
    if all(
        local_file.exists() and (get_md5(local_file) == hexdigest)
        for local_file, hexdigest in test_files.items()
    ):
        skipped += len(test_files)
    else:
        print("Starting data download.")
        for filename in tqdm(test_files.keys()):
            url = download_url + Path(filename).name
            download_file(url, filename)
            downloaded += 1

    for v in "downloaded", "skipped":
        if eval(v):
            msg = v.capitalize() + f" {eval(v)} files."
            print(msg)

    return


test_files = {
    "RL1_S001__insert_270_short_R1.fq.gz": "09a0e7403a9cea7c2d3c355db8c649cc",
    "RL1_S001__insert_270_short_R2.fq.gz": "1c1178a5686422d9f5228bdc3a87f875",
    "RL1_S002__insert_270_short_R1.fq.gz": "2ebfa1a25b3942bac2bda7f20dea97d7",
    "RL1_S002__insert_270_short_R2.fq.gz": "6d1dfa872190e9bf3586328a81f307e3",
    "RL2_S001__insert_270_short_R1.fq.gz": "9528b0e89ee769e236e720d2e8e71c33",
    "RL2_S001__insert_270_short_R2.fq.gz": "66684a2b0b0a869a6016a8a615dd7932",
    "RL2_S002__insert_270_short_R1.fq.gz": "0d85c5f8eb8fedc428e3804c3f8a0134",
    "RL2_S002__insert_270_short_R2.fq.gz": "c46563a7b9d0e9be7cba475c5514f8b7",
}
download_url = "https://github.com/vinisalazar/mg-example-data/raw/main/data/"


def main(args):
    directory = args.directory
    coassembly = args.coassembly
    cores = args.cores
    max_mb = args.max_mb
    join = args.join_units
    samples_file = "samples.csv"
    create_input_table_args = Namespace()
    create_input_table_args.input_dir = directory
    create_input_table_args.output_file = samples_file
    create_input_table_args.join_units = join
    confirm = args.confirm
    dry_run = args.dry_run
    extras = args.extras
    profile = args.profile

    # Start execution
    download_data(directory, test_files, download_url)
    print("\nCreating input table for test files.\n")
    create_input_table(create_input_table_args)
    print()
    print("Starting Snakemake.")
    print(
        "This may require the installation of conda environments which should take a while.\n"
    )
    if not confirm and not dry_run:
        confirm_message(cores, max_mb)

    cmd = f"""
    snakemake --snakefile {snakefile}           \
              --configfile {test_config}        \
              --cores {cores}                   \
              -p -r                             \
              --use-conda                       \
              --wrapper-prefix {wrapper_prefix}
    """

    for arg in ("max_mb", "coassembly"):
        if value := eval(arg):
            cmd += f" --config {arg}={value} "

    if profile:
        cmd += f" --profile {profile} "

    if dry_run:
        cmd += f" --dry-run --conda-frontend conda "

    cmd += f" {extras}"

    test_complete_message = (
        """
    Test complete!

    This means that you can use this directory to run your actual analysis.
    All necessary software is in the .snakemake/conda/ directory, and databases are in the data/ directory.
    Simply delete the output/ directory and you're good to go.
    """
        if not dry_run
        else """
    Dry run complete. You can run the full testing suite by repeating without the `--dry-run` flag.
    """
    )

    try:
        retcode = check_call(cmd.split())
    except CalledProcessError as e:
        retcode = e.returncode
        raise
    except Exception as e:
        retcode = 1
        raise

    get_successful_completion(retcode, dedent(test_complete_message))
