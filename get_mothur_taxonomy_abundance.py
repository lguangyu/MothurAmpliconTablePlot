#!/usr/bin/env python3

import argparse
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser("calculate relative abundances of taxons at "
		"given taxon rank")
	ap.add_argument("-s", "--shared", type = str, required = True,
		metavar = "shared",
		help = "the otu count table in mothur SOP output (required)")
	ap.add_argument("-t", "--taxonomy", type = str, required = True,
		metavar = "tax",
		help = "the otu taxonomy table in mothur SOP output (required)")
	ap.add_argument("-r", "--tax-rank", type = str, default = "genus",
		choices = pylib.taxonomy.MothurLineage.RANK,
		help = "taxonomic rank to extrack (default: genus)")
	ap.add_argument("--min-bootstrap", type = int, default = 0,
		metavar = "int",
		help = "minimum bootstrap value (included) to account for taxonomic "
			"classification; 0 means accepting everything (default: 0)")
	ap.add_argument("-o", "--output", type = str, default = "-",
		metavar = "tsv",
		help = "output taxonomic abundance table (default: <stdout>)")
	ap.add_argument("--transpose-output", action = "store_true",
		help = "transpose output abundance table; by default orientation, "
			"each row is a sample and each column is a taxon; row/col will be "
			"transposed if this flag is set (default: off)")
	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout
	return args


def main():
	args = get_args()
	# load data
	otu_count	= pylib.shared_table.MothurSharedTable.from_file(args.shared)
	taxonomy	= pylib.taxonomy.MothurOtuTaxonomyDict.from_file(args.taxonomy)
	# calculate taxonomic abundances
	otu_abund	= otu_count.normalize_rows()
	tax_abund	= otu_abund.group_by_taxonomy(taxonomy,
		tax_rank = args.tax_rank, min_bootstrap = args.min_bootstrap)
	# save
	tax_abund.save(args.output, transpose = args.transpose_output)
	return


if __name__ == "__main__":
	main()
