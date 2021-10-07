#!/usr/bin/env python3

import argparse
import functools
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import sklearn
import sklearn.decomposition
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser("PCA plot based on relative abundances "
		"analysis")
	ap.add_argument("input", type = str, nargs = "?", default = "-",
		help = "input taxonomic abundance table (default: <stdin>)")
	ap.add_argument("--input-transposed", action = "store_true",
		help = "read input table as if transposed (default: off)")
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


class PcaAnalyzer(sklearn.decomposition.PCA):
	def analyze(self, abund_table: pylib.shared_table.TaggedTable):
		if not isinstance(abund_table, pylib.shared_table.TaggedTable):
			raise TypeError("abund_table must be TaggedTable, not '%s'"\
				% type(abund_table).__name__)
		self.abund_table = abund_table
		self.transformed_abund_table = self.fit_transform(abund_table.data)
		return

	@property
	def var_ratio(self):
		return self.explained_variance_ratio_


class AspectInvariantFancyArrowPatch(matplotlib.patches.FancyArrowPatch):
	@functools.wraps(matplotlib.patches.FancyArrowPatch.draw)
	def draw(self, renderer, *ka, **kw):
		#self.set_mutation_aspect(numpy.sqrt(1 / get_data_aspect_ratio(self.axes)))
		return super().draw(renderer, *ka, **kw)


def setup_layout(figure):
	layout = dict(figure = figure)

	# margins
	left_margin_inch	= 1.0
	right_margin_inch	= 0.3
	top_margin_inch		= 0.6
	bottom_margin_inch	= 0.7

	# pca dims
	pca_width_inch		= 6.0
	pca_height_inch		= 6.0

	# figure size
	figure_width_inch	= left_margin_inch + pca_width_inch\
		+ right_margin_inch
	figure_height_inch	= bottom_margin_inch + pca_height_inch\
		+ top_margin_inch
	figure.set_size_inches(figure_width_inch, figure_height_inch)

	# pca axes
	pca_left	= left_margin_inch / figure_width_inch
	pca_bottom	= bottom_margin_inch / figure_height_inch
	pca_width	= pca_width_inch / figure_width_inch
	pca_height	= pca_height_inch / figure_height_inch
	pca_axes	= figure.add_axes([pca_left, pca_bottom,
		pca_width, pca_height])
	layout["pca"] = pca_axes

	# style axes
	for sp in pca_axes.spines.values():
		sp.set_visible(False)
	pca_axes.tick_params(
		left = False, labelleft = True,
		right = False, labelright = False,
		top = False, labeltop = False,
		bottom = False, labelbottom = True)
	pca_axes.set_facecolor("#F0F0F8")
	pca_axes.grid(linestyle = "-", linewidth = 1.5, color = "#FFFFFF")
	pca_axes.axhline(0.0, linestyle = "-", linewidth = 1.0, color = "#808080")
	pca_axes.axvline(0.0, linestyle = "-", linewidth = 1.0, color = "#808080")

	return layout


def max_norm(mat, *, ord = 2, axis = 1, **kw):
	return max(numpy.linalg.norm(mat, ord = ord, axis = axis, **kw))


def plot(png, analyzed_pca: PcaAnalyzer, *, biplot_n_feat = 5, title = ""):
	#assert isinstance(analyzed_pca, PcaAnalyzer)
	# layout
	figure = matplotlib.pyplot.figure()
	layout = setup_layout(figure)

	# plot pca
	axes = layout["pca"]
	# rescale pca scatter to 5
	pca_scatter_scale = 5.0 / max_norm(analyzed_pca.transformed_abund_table)
	x = analyzed_pca.transformed_abund_table[:, 0] * pca_scatter_scale
	y = analyzed_pca.transformed_abund_table[:, 1] * pca_scatter_scale
	axes.scatter(x, y, marker = "o", s = 40, edgecolors = "#4040FF",
		facecolors = "#FFFFFF80", linewidths = 1.0, zorder = 3)

	# add text
	for tx, ty, s in zip(x, y, analyzed_pca.abund_table.row_tag):
		axes.text(tx, ty, s, fontsize = 10, color = "#000000")

	# add biplot
	biplot_n_feat = 5
	biplot_scale = 3.0 / max_norm(analyzed_pca.components_[:, :biplot_n_feat],
		axis = 0)
	for i in range(biplot_n_feat):
		xy = analyzed_pca.components_[:, i] * biplot_scale
		arrow = AspectInvariantFancyArrowPatch(posA = (0, 0), posB = xy,
			arrowstyle = "-|>, head_length = 5.0, head_width=2.5",
			linewidth = 1.0, color = "#000000", zorder = 4)
		axes.add_patch(arrow)
		axes.text(*xy, analyzed_pca.abund_table.col_tag[i],	fontsize = 14,
			color = "#FF4040", zorder = 3,
			horizontalalignment = "center", verticalalignment = "center")

	# misc
	axes.set_xlabel("PC1 (%.1f%%)" % (analyzed_pca.var_ratio[0] * 100),
		fontsize = 12)
	axes.set_ylabel("PC2 (%.1f%%)" % (analyzed_pca.var_ratio[1] * 100),
		fontsize = 12)

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
	# pca
	pca = PcaAnalyzer(n_components = 2, whiten = False)
	pca.analyze(tax_abund)
	# plot
	plot(args.plot, pca, title = args.plot_title)
	return


if __name__ == "__main__":
	main()
