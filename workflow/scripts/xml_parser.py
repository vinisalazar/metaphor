import re
import sys
import shelve
import logging
import argparse
import xml.etree.ElementTree as ET

import pandas as pd


def main(args):

	tax = shelve.open(args.db)
	df = pd.DataFrame(columns=("code", "sample_name"))

	sample_names = []
	for xml in args.xmls:
		sample_name = xml.split("/")[-1].split(".")[0]
		sample_names.append(sample_name)

		root = ET.parse(xml).getroot()
		for hit in root.findall("./BlastOutput_iterations/Iteration/Iteration_hits/Hit"):
			hit_text = hit.find("Hit_id").text
			code = hit_text.split(":")[0]
			df = df.append( {"code":code, "sample_name":sample_name}, ignore_index=True)

	logging.info(f"Writing to {args.outfile}")
	with open(args.outfile, "w") as f:
		columns = ["Name", "Kegg Code", "Lineage"] + sample_names
		# import pdb; pdb.set_trace()
		codes = sorted(df["code"].unique())

		if codes:
			for code in codes:
				if code in tax.keys():
					lineage = tax[code]
					name = lineage.split(";")[-1]
					name = name[ len(code)+2:].strip()
				else:
					lineage = "unknown"
					name = "Unknown"

				f.write( "\t".join( [name, code, lineage] )  )
				for sample_name in sample_names:
					filtered = df[ (df["sample_name"] == sample_name) & (df["code"] == code) ]
					count = len(filtered.index)
					f.write( f"\t{count}" )

				f.write("\n")
		else:
			input_files = "\n\t\t".join(args.xmls)
			logging.debug(
				f"The output appears to be empty. Please check the input files:\n{input_files}."
			)


def run(snakemake):
	args = argparse.Namespace()
	args_dict = vars(args)

	for rule_input in ("xmls",):
		args_dict[rule_input] = snakemake.input[rule_input]

	for rule_param in ("db",):
		args_dict[rule_param] = snakemake.params[rule_param]

	for rule_output in ("outfile",):
		args_dict[rule_output] = snakemake.output[rule_output]

	main(args)


# Logging config must be at the beginning
logging.basicConfig(
    filename=str(snakemake.log),
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
)
logging.info(f"Starting script {__file__.split('/')[-1]}.")
logging.debug(f"Full script path: {__file__}")
run(snakemake)
logging.info(f"Done.")