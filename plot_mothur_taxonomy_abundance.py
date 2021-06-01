#!/usr/bin/env python3

import argparse
import itertools
import matplotlib
import matplotlib.pyplot
import numpy
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser("plot relative abundances table after abudance"
		"analysis")
	ap.add_argument("input", type = str, nargs = "?", default = "-",
		help = "input taxonomic abundance table (default: <stdin>)")
	ap.add_argument("--input-transposed", action = "store_true",
		help = "read input table as if transposed (default: off)")
	ap.add_argument("--taxon-list", type = str,
		metavar = "txt",
		help = "if provided, only plot taxons provided by this list")
	ap.add_argument("--max-n-taxons", type = int, default = 20,
		metavar = "int",
		help = "plot at most this many taxons (default: 20); ignored if "
			"--taxon-list is used")
	ap.add_argument("--plot-percent", action = "store_true",
		help = "plot in percent instead of fraction (default: no)")
	ap.add_argument("--plot-title", type = str, default = "",
		metavar = "str",
		help = "set to include a plot title (default: <empty>)")
	ap.add_argument("-p", "--plot", type = str, default = "-",
		metavar = "png",
		help = "plot image file (default: <stdout>)")
	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.plot == "-":
		args.plot = sys.stdout.buffer
	return args


def get_plot_abund_table(tax_abund, *, taxon_list = None, max_n_taxons = 20)\
		-> pylib.shared_table.TaggedTable:
	return _filter_taxons_by_list(tax_abund, taxon_list) if taxon_list\
		else _filter_taxons_by_max_n(tax_abund, max_n_taxons)


def _filter_taxons_by_list(tax_abund, taxon_list):
	with pylib.file_util.get_fp(taxon_list, "r") as fp:
		user_tax = fp.read().splitlines()
	user_tax_order = {k:i for (i, k) in enumerate(user_tax)}
	# filter taxons
	data = numpy.zeros((tax_abund.nrow, len(user_tax)), dtype = tax_abund.dtype)
	for i, t in enumerate(tax_abund.col_tag):
		if t in user_tax_order:
			data[:, user_tax_order[t]] = tax_abund.data[:, i]
	ret = pylib.shared_table.TaggedTable(row_tag = tax_abund.row_tag,
		col_tag = numpy.asarray(user_tax, dtype = object), data = data)
	return ret


def _filter_taxons_by_max_n(tax_abund, max_n_taxons):
	nrow, ncol = tax_abund.shape
	# if too many taxons (> max_n_taxons), combine everthing ranked lower
	# into '[all others]' category
	# rows: samples; cols: taxons
	if max_n_taxons < ncol:
		ncol	= max_n_taxons + 1
		col_tag	= numpy.empty(ncol, dtype = object)
		col_tag[:max_n_taxons]	= tax_abund.col_tag[:max_n_taxons]
		col_tag[max_n_taxons]	= "[all others classified]"
		data	= numpy.hstack([tax_abund.data[:, :max_n_taxons],
			tax_abund.data[:, max_n_taxons:].sum(axis = 1, keepdims = True)])
	else:
		col_tag	= tax_abund.col_tag
		data	= tax_abund.data
	ret = pylib.shared_table.TaggedTable(row_tag = tax_abund.row_tag,
		col_tag = col_tag, data = data)
	return ret


def setup_layout(figure, nrow, ncol):
	layout = dict(figure = figure)

	# margins
	left_margin_inch	= 3.0
	right_margin_inch	= 0.5
	top_margin_inch		= 0.6
	bottom_margin_inch	= 1.6

	# heatmap dims
	cell_width_inch		= 0.4
	cell_height_inch	= 0.2
	heatmap_width_inch	= cell_width_inch * ncol
	heatmap_height_inch	= cell_height_inch * nrow

	# figure size
	figure_width_inch	= left_margin_inch + heatmap_width_inch\
		+ right_margin_inch
	figure_height_inch	= bottom_margin_inch + heatmap_height_inch\
		+ top_margin_inch
	figure.set_size_inches(figure_width_inch, figure_height_inch)

	# heatmap axes
	heatmap_left	= left_margin_inch / figure_width_inch
	heatmap_bottom	= bottom_margin_inch / figure_height_inch
	heatmap_width	= heatmap_width_inch / figure_width_inch
	heatmap_height	= heatmap_height_inch / figure_height_inch
	heatmap_axes	= figure.add_axes([heatmap_left, heatmap_bottom,
		heatmap_width, heatmap_height])
	layout["heatmap"] = heatmap_axes

	# style axes
	for sp in heatmap_axes.spines.values():
		sp.set_visible(False)
	heatmap_axes.tick_params(
		left = False, labelleft = True,
		right = False, labelright = False,
		top = False, labeltop = False,
		bottom = False, labelbottom = True)

	# colorbar
	colorbar_width_inch		= 0.75 * left_margin_inch
	colorbar_height_inch	= 0.20 * bottom_margin_inch
	colorbar_width	= colorbar_width_inch / figure_width_inch
	colorbar_height	= colorbar_height_inch / figure_height_inch
	colorbar_left	= (left_margin_inch / figure_width_inch\
		- colorbar_width) / 2
	colorbar_bottom	= (bottom_margin_inch / figure_height_inch\
		- colorbar_height) / 2
	colorbar_axes = figure.add_axes([colorbar_left, colorbar_bottom,
		colorbar_width, colorbar_height])
	layout["colorbar"] = colorbar_axes

	return layout


def plot(png, abund_table, *, title = "", percent = False):
	# gather info
	n_sample, n_taxon = abund_table.shape
	sample_list	= abund_table.row_tag
	taxon_list	= abund_table.col_tag
	data		= abund_table.data.T
	if percent:
		vmax	= 100
		vformat	= "%.1f"
		min_val	= 0.1
		min_txt	= "min."
	else:
		vmax	= 1.0
		vformat	= "%.2f"
		min_val	= 0.01
		min_txt	= "min."

	# layout
	figure = matplotlib.pyplot.figure()
	layout = setup_layout(figure, n_taxon, n_sample)

	# plot table heatmap
	axes = layout["heatmap"]
	cmap = matplotlib.pyplot.get_cmap("Spectral_r")
	heatmap = axes.pcolor(data * vmax, cmap = cmap, vmin = 0, vmax = vmax)
	# add text on heatmap
	for r, c in itertools.product(range(n_taxon), range(n_sample)):
		val = data[r, c] * vmax
		if val == 0:
			text = ""
		elif val < min_val:
			text = min_txt
		else:
			text = vformat % val
		color = ("#000000" if (val <= 0.7 * vmax) and (val >= 0.2 * vmax)\
			else "#FFFFFF")
		axes.text(x = c + 0.5, y = r + 0.5, s = text, fontsize = 8,
			color = color,
			horizontalalignment = "center", verticalalignment = "center")
	# add lines
	grid_props = dict(linestyle = "-", linewidth = 0.5, color = "#FFFFFF",
		clip_on = False)
	for i in numpy.arange(n_taxon + 1):
		axes.axhline(i, **grid_props)
	for i in numpy.arange(n_sample + 1):
		axes.axvline(i, **grid_props)
	# misc
	xticks = numpy.arange(n_sample) + 0.5
	axes.set_xticks(xticks)
	axes.set_xticklabels(sample_list, family = "sans-serif",
		rotation = 90, fontsize = 12, fontweight = "medium",
		horizontalalignment = "center", verticalalignment = "top")
	axes.set_xlim(0, n_sample)
	yticks = numpy.arange(n_taxon) + 0.5
	axes.set_yticks(yticks)
	axes.set_yticklabels(taxon_list, family = "sans-serif",
		fontsize = 10,
		horizontalalignment = "right", verticalalignment = "center")
	axes.set_ylim(n_taxon, 0)

	# colorbar
	axes = layout["colorbar"]
	colorbar = figure.colorbar(heatmap, cax = axes, orientation = "horizontal")
	colorbar.set_label("Relative abundance (%s)" % ("%" if percent else "frac"))
	colorbar.outline.set_visible(False)

	# title
	figure.suptitle(title, fontsize = 18, fontweight = "medium",
		horizontalalignment = "center", verticalalignment = "top")
	# save
	matplotlib.pyplot.savefig(png, dpi = 300)
	matplotlib.pyplot.close()
	return


def main():
	args = get_args()
	# load data
	tax_abund	= pylib.shared_table.TaggedTable.from_file(args.input,
		transposed = args.input_transposed)
	# plot
	plot_abund	= get_plot_abund_table(tax_abund, taxon_list = args.taxon_list,
		max_n_taxons = args.max_n_taxons)
	plot(args.plot, plot_abund, title = args.plot_title,
		percent = args.plot_percent)
	return


if __name__ == "__main__":
	main()
