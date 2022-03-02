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


import sys
from argparse import Namespace
from pathlib import Path
from hashlib import md5

import requests
from tqdm import tqdm
from snakemake import snakemake

from metaphor import (
    snakefile,
    test_config,
    ascii_art,
    get_successful_completion,
    wrapper_prefix,
)
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
    print("Starting data download.")
    downloaded, skipped = 0, 0
    for filename, hexdigest in tqdm(test_files.items()):
        local_file = test_directory.joinpath(filename)
        if (not local_file.exists()) or (get_md5(local_file) != hexdigest):
            url = download_url + filename
            download_file(url, local_file)
            downloaded += 1
        else:
            # print(f"File '{filename}' exists and hash is correct. Skipping download.")
            skipped += 1

    for v in "downloaded", "skipped":
        if eval(v):
            msg = v.capitalize() + f" {eval(v)} files."
            print(msg)


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
    cores = int(args.cores)
    mem_mb = int(args.mem_mb)
    samples_file = str(Path(directory).joinpath("samples.csv"))
    create_input_table_args = Namespace()
    create_input_table_args.input_dir = directory
    create_input_table_args.output_file = samples_file
    create_input_table_args.join_units = False
    conda_prefix = directory if args.remove_conda else None

    # Start execution
    download_data(directory, test_files, download_url)
    print("\nCreating input table for test files.\n")
    create_input_table(create_input_table_args)
    print()
    print("Starting Snakemake.")
    print(
        "This may require the installation of conda environments which should take a while.\n"
    )
    if not args.confirm:
        confirm = input(
            f"Snakemake will start with {cores} cores and {mem_mb} MB RAM. Ok to continue? [y/N]"
        )
        if confirm.lower() != "y":
            print("Metaphor test cancelled.")
            sys.exit()
    print(ascii_art)
    sm_exit = snakemake(
        snakefile=snakefile,
        configfiles=[
            test_config,
        ],
        config={
            "samples": samples_file,
            "coassembly": coassembly,
        },
        cores=cores,
        resources={"mem_mb": mem_mb},
        use_conda=True,
        conda_prefix=conda_prefix,
        printshellcmds=True,
        wrapper_prefix=wrapper_prefix,
    )
    get_successful_completion(sm_exit, "Test complete.")
