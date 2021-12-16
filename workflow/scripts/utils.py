"""
Helper functions that are common across scripts.
"""
import logging
import argparse


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


def driver(main_fn, snakemake_obj, parse_args_fn=parse_snakemake_args):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    if snakemake_obj:
        fh = logging.FileHandler(str(snakemake_obj.log), encoding="utf-8")
        logging.getLogger().addHandler(fh)
        args = parse_args_fn(snakemake_obj)
    else:
        args = parse_args_fn()
    try:
        logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
        logging.debug(f"Full script path: '{__file__}'.")
        main_fn(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())
