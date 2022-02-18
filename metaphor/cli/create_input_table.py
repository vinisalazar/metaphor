#!/usr/bin/env python

__doc__ = """
Create Metaphor input table from directory of fastq files.
"""


import re
import difflib
import argparse
import pandas as pd
from collections import namedtuple
from pathlib import Path

fastq_dir = "fastq/"
s_as_units = False
unit_expr = "S[0-9][0-9][0-9]"


def is_valid_ext(filepath):
    if not isinstance(filepath, str):
        filepath = str(filepath)
    valid_exts = [".fq", ".fastq", ".fna", ".fa", "fasta"]
    valid_exts += [i + ".gz" for i in valid_exts]
    for ext in valid_exts:
        if filepath.endswith(ext):
            return True

    return False


def get_sample_unit(s):
    unit_split = re.split("S[0-9][0-9][0-9]", s)
    if len(unit_split) == 1:
        return "single"
    else:
        return re.search(unit_expr, s).group()


def get_files(fastq_dir):
    files = [file.name for file in Path(fastq_dir).iterdir() if is_valid_ext(file)]
    n_files = len(files)
    print(f"{n_files} files detected:")
    print("\t" + "\n\t".join(files) + "\n")
    return files


def create_file_lists(files):
    R1, R2, unpaired = [], [], []

    for file in files:
        if re.search("_R1", file):
            R1.append(file)
        elif re.search("_R2", file):
            R2.append(file)
        else:
            unpaired.append(file)
    R1, R2, unpaired = sorted(R1), sorted(R2), sorted(unpaired)

    return R1, R2, unpaired


def create_paired_files_dict(R1, R2, n_files, s_as_units):
    PairedSample = namedtuple("PairedSample", ["unit", "R1", "R2"])
    paired_samples_dict = {}

    for f, r in zip(R1, R2):
        match = difflib.SequenceMatcher(None, f, r)
        size = match.find_longest_match().size
        sample_name = f[:size].replace("_R", "")
        unit = get_sample_unit(sample_name) if s_as_units else "single"
        sample_name = (
            sample_name.replace(unit, "").replace("__", "")
            if s_as_units
            else sample_name
        )

        if s_as_units:
            if sample_name in paired_samples_dict.keys():
                paired_samples_dict[sample_name].unit.append(unit)
                paired_samples_dict[sample_name].R1.append(f)
                paired_samples_dict[sample_name].R2.append(r)
            else:
                paired_samples_dict[sample_name] = PairedSample(
                    # black, go home, you're drunk
                    [
                        unit,
                    ],
                    [
                        f,
                    ],
                    [
                        r,
                    ],
                )
        else:
            paired_samples_dict[sample_name] = PairedSample(
                unit,
                f,
                r,
            )

    n_samples = len(paired_samples_dict)
    files_per_sample = n_files / n_samples
    files_per_sample = (
        int(files_per_sample) if files_per_sample.is_integer() else files_per_sample
    )
    print(f"{n_samples} samples detected, {files_per_sample} files per sample.")
    for k, v in paired_samples_dict.items():
        fmt_file = (
            lambda file: "\n\t".join(file) if isinstance(file, list) else "\t" + file
        )
        print(k + ":" + "\n\t" + fmt_file(v.R1) + "\n\t" + fmt_file(v.R2))

    return paired_samples_dict


def create_paired_samples_df(paired_samples_dict, s_as_units):
    df = pd.DataFrame(paired_samples_dict).T
    df.columns = ["unit_name", "R1", "R2"]
    df.index.name = "sample_name"

    if s_as_units:
        df = df.explode(df.columns.to_list())

    return df


def process_final_df(dfs, fastq_dir):
    final_df = pd.concat(dfs)
    extra_columns = [
        "metaquast_reference",
    ]
    for col in extra_columns:
        final_df[col] = ""

    # Add absolute path to files
    for col in ["R1", "R2"]:
        final_df[col] = final_df[col].apply(
            lambda file: str(Path(fastq_dir).absolute().joinpath(file))
        )

    return final_df


def main(args):
    input_dir, join_units = args.input_dir, args.join_units
    files = get_files(input_dir)
    R1, R2, unpaired = create_file_lists(files)
    paired_files_dict = create_paired_files_dict(
        R1, R2, len(files), s_as_units=join_units
    )
    paired_files_df = create_paired_samples_df(paired_files_dict, s_as_units=join_units)

    # TODO: add single-end, interleaved dfs here.
    dfs = [
        paired_files_df,
    ]
    final_df = process_final_df(dfs, input_dir)

    outfile = "./samples.csv" if args.output_file is None else args.output_file
    final_df.to_csv(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input-dir", help="Input directory containing FASTQ files."
    )
    parser.add_argument(
        "-j",
        "--join-units",
        help="If this option is on, files with the same preffix but with "
        "S001, S002, S00N distinctions in the filenames will be treated as different units of the same sample, "
        "i.e. they will be joined into a single file.",
        action="store_true",
    )
    parser.add_argument("-o", "--output-file", help="Path to output file.")
    args = parser.parse_args()
    main(args)
