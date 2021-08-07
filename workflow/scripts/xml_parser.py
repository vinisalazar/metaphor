import xml.etree.ElementTree as ET
import sys
import json
import re
import argparse
import pandas as pd


def search( node, query, parent_lineage=[] ):
	code = node["name"].split()[0]

	# Append this node to the lineage
	my_linage = list(parent_lineage) + [node["name"]]

	if code == query:
		return my_linage

	# if it doesn't have the code, then search children
	if "children" in node:
		for child in node["children"]:
			result = search( child, query, my_linage )
			if result != False:
				return result

	# if no children then return false
	return False


parser = argparse.ArgumentParser()
parser.add_argument('lineage', help="The KEGG linage file.")
parser.add_argument('xmls', help="The XML files to search for the hits.",  nargs='+')
parser.add_argument('--outfile', help="The name of the output file. Default: OTU.tsv", default="OTU.tsv")

args = parser.parse_args()


with open(args.lineage) as f:
	tax = json.load(f)


df = pd.DataFrame()

sample_names = []
for xml in args.xmls:
	sample_name = xml.split(".")[0]
	sample_names.append(sample_name)

	root = ET.parse(xml).getroot()
	for hit in root.findall("./BlastOutput_iterations/Iteration/Iteration_hits/Hit"):
		hit_text = hit.find("Hit_id").text
		code = hit_text.split(":")[0]
		df = df.append( {'code':code, "sample_name":sample_name}, ignore_index=True)

print(f"Writing to {args.outfile}")
with open(args.outfile, "w") as f:
	columns = ["Name", "Kegg Code", "Lineage"] + sample_names
	print("\t".join(columns), file=f)

	codes = sorted(df.code.unique())

	for code in codes:
		lineage = search( tax, code )
		lineage.pop(0) # Get rid of root
		lineage_text = "; ".join(lineage)
		name = lineage[-1]
		name = name[ len(code):].strip()

		f.write( "\t".join( [name, code, lineage_text] )  )
		for sample_name in sample_names:
			filtered = df[ (df.sample_name == sample_name) & (df.code == code) ]
			count = len(filtered.index)
			f.write( f"\t{count}" )

		f.write("\n")
