"""
Helper functions that are common across scripts.
"""


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


def driver(parse_snakemake_args_fn=parse_snakemake_args):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    if "snakemake" in locals():
        fh = logging.FileHandler(str(snakemake.log), encoding="utf-8")
        logging.getLogger().addHandler(fh)
        args = parse_snakemake_args(snakemake)
    else:
        args = parse_args()
    try:
        logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
        logging.debug(f"Full script path: '{__file__}'.")
        main(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())
