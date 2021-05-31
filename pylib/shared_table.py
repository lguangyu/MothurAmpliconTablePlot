#!/usr/bin/env python3

import collections
import copy
import numpy
from . import file_util


class TaggedTable(object):
	def __init__(self, *ka, row_tag, col_tag, data, **kw):
		super().__init__(*ka, **kw)
		if not (len(row_tag) == len(data)):
			raise ValueError("row_tag and data (1st dimension) must be of the "
				"same length")
		if not (len(col_tag) == len(data[0])):
			raise ValueError("col_tag and data (2nd dimension) must be of the "
				"same length")
		self.row_tag	= row_tag
		self.col_tag	= col_tag
		self.data		= numpy.asarray(data)
		return

	def normalize_rows(self, inplace = False):
		ret = self if inplace else copy.deepcopy(self)
		ret.data = ret.data / ret.data.sum(axis = 1, keepdims = True)
		return ret

	@property
	def shape(self):
		return self.data.shape
	@property
	def nrow(self):
		return len(self.row_tag)
	@property
	def ncol(self):
		return len(self.col_tag)
	@property
	def dtype(self):
		return self.data.dtype

	@classmethod
	def _get_sorted_index(cls, data, method):
		if method == "by_average":
			sort_data = data.mean(axis = 0)
			idx = numpy.flip(sort_data.argsort())
		elif method == "by_rank":
			nrow, ncol = data.shape
			ranks = numpy.zeros(data.shape, dtype = int)
			for r, d in zip(ranks, data):
				r[d.argsort()] = numpy.arange(len(ncol), 0, -1)
			idx = cls._get_sorted_index(ranks, "by_average")
		else:
			raise ValueError("invalid sort method: '%s'" % method)
		return idx

	def sort_cols_descend(self, method = "by_average", *, inplace = False):
		"""
		sort columns in descending order based on data

		method - by_average:
		sort based on column-wise average of each column

		method - by_rank:
		sort based on average rank of each column; therefore, needs to first
		calcualte the ranks of each column
		"""
		ret = self if inplace else copy.deepcopy(self)
		idx = ret._get_sorted_index(ret.data, method)
		ret.col_tag = ret.col_tag[idx]
		ret.data	= ret.data[:, idx]
		return ret

	@classmethod
	def from_file(cls, file, *, data_type = float, delimiter = "\t",
			transposed = False):
		"""
		parse mothur output shared file as tabular data, assuming 1st row as
		column tags and 1st column as row tags
		"""
		raw = numpy.loadtxt(file, dtype = object, delimiter = delimiter,
			ndmin = 2)
		if transposed:
			raw = raw.T
		new = cls(
			row_tag		= raw[1:, 0],
			col_tag		= raw[0, 1:],
			data		= raw[1:, 1:].astype(data_type),
		)
		return new

	def save(self, file, *, delimiter = "\t", transpose = False):
		sf = numpy.empty((len(self.row_tag) + 1, len(self.col_tag) + 1),
			dtype = object)
		sf[0, 0]	= ""
		sf[1:, 0]	= self.row_tag
		sf[0, 1:]	= self.col_tag
		sf[1:, 1:]	= self.data
		numpy.savetxt(file, sf.T if transpose else sf, fmt = "%s",
			delimiter = delimiter)
		return


class MothurSharedTable(TaggedTable):
	def __init__(self, *ka, label, group, num_otus, otu, data, **kw):
		if not (len(label) == len(group) == len(num_otus) == len(data)):
			raise ValueError("label, group, num_otus and data (1st dimension) "
				"must be of the same length")
		if not (len(otu) == len(data[0])):
			raise ValueError("otu and data (2nd dimension) "
				"must be of the same length")
		super().__init__(*ka, row_tag = group, col_tag = otu, data = data, **kw)
		self.label		= label
		self.num_otus	= num_otus
		return

	@property
	def group(self):
		return self.row_tag
	@property
	def otu(self):
		return self.col_tag

	@classmethod
	def from_file(cls, file, *, delimiter = "\t"):
		"""
		parse mothur output shared file as tabular data
		"""
		raw = numpy.loadtxt(file, dtype = object, delimiter = delimiter)
		new = cls(
			label		= raw[1:, 0],
			group		= raw[1:, 1],
			num_otus	= raw[1:, 2].astype(int),
			otu			= raw[0, 3:],
			data		= raw[1:, 3:].astype(int),
		)
		return new

	def group_by_taxonomy(self, otu_tax_dict: dict, tax_rank) -> TaggedTable:
		"""
		group otus into taxons at <tax_rank>, then sum otu counts into taxon
		counts

		to do this,
		1. create an intermediate dictionary of tax_name->data_per_tax
		2. iterate over all otus, use otu.unique_taxon_name as key, sum data
		   into the dictionary
		3. stack/combine the data in dictionary
		"""
		nrow, ncol = self.data.shape
		tax_data = collections.defaultdict(lambda : numpy.zeros(nrow,
			dtype = self.dtype))
		for i, o in enumerate(self.otu):
			tax = otu_tax_dict[o].taxonomy[tax_rank]
			if not tax:
				continue # if not classified at this level
			tax_data[tax.unique_taxon_name] += self.data[:, i]
		# make return object
		col_tag = numpy.asarray(list(tax_data.keys()), dtype = object)
		data = numpy.empty((nrow, len(col_tag)), dtype = self.dtype)
		# use loop instead of hstack/vstack to avoid empty tax_data error
		# such case will occur when no otu classified at a certain level
		# for example, default mothur SOP classify as deep as genus level
		# therefore error will raise when try to stat at species level
		for i, k in enumerate(col_tag):
			data[:, i] = tax_data[k]
		ret = TaggedTable(row_tag = self.row_tag, col_tag = col_tag,
			data = data)
		ret.sort_cols_descend(inplace = True)
		return ret

