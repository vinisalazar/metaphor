"""
Helper functions that are common across scripts.
"""
import logging
import argparse
import traceback


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for attr in "input", "output", "params", "wildcards":
        try:
            for k, v in getattr(snakemake, attr).items():
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
        logging.info("Done.\n")
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())


def write_dfs(absolute, absolute_out, relative_out, relative=None):
    if relative is None:
        relative = round(absolute / absolute.sum(), 4)
    for kind in "absolute", "relative":
        df, outfile = eval(kind), eval(f"{kind}_out")
        df.columns = [i.replace("-to-genes.sorted.bam", "") for i in df.columns]
        df.to_csv(outfile, sep="\t")
        logging.info(f"Wrote {len(df)} rows to '{outfile}'.")


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
