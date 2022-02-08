"""
Helper functions that are common across scripts.
"""
import logging
import argparse
import traceback


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for directive in "input", "output", "params":
        try:
            for k, v in getattr(snakemake, directive).items():
                args_dict[k] = v
        except AttributeError:
            pass

    return args


def driver(main_fn, snakemake_obj, file, parse_args_fn=None):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    if snakemake_obj:
        fh = logging.FileHandler(str(snakemake_obj.log), encoding="utf-8")
        logging.getLogger().addHandler(fh)
        if not parse_args_fn:
            parse_args_fn = parse_snakemake_args
        args = parse_args_fn(snakemake_obj)
    else:
        args = parse_args_fn()
    try:
        logging.info(f"Starting script '{file.split('/')[-1]}'.")
        logging.info(f"Full script path: '{file}'.\n")
        main_fn(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())


cog_csv_names = [
    "Gene ID (GenBank or ad hoc)",
    "NCBI Assembly ID",
    "Protein ID",
    "Protein length",
    "COG footprint coordinates",
    "Length of the COG footprint on the proteins",
    "COG ID",
    "reserved",
    "COG membership class",
    "PSI-BLAST bit score",
    "PSI-BLAST e-value",
    "COG profile length",
    "Protein footprint coordinates on the COG profile",
]
