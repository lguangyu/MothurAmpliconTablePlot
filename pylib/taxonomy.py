#!/usr/bin/env python3

import re
from . import file_util


class MothurTaxonomy(object):
	"""
	MothurTaxonomy deals with the most basic taxonomy data structure in mothur
	taxonomy output file;
	"""
	TAX_REGEX = re.compile(r"^(.*)\((\d+)\)$")

	def __init__(self, taxon_name: str, bootstrap: int = 100, *ka, **kw):
		super().__init__(*ka, **kw)
		self.taxon_name	= taxon_name
		self.bootstrap	= bootstrap
		return

	@property
	def unique_taxon_name(self) -> str:
		"""
		some raw taxon names may be ambiguous (e.g. uncultured)
		unique_taxon_name is the way to resolve this naming collision

		unique_taxon_name must be set using self.set_unique_taxon_name()
		"""
		return getattr(self, "_uniq_name", None) or self.taxon_name

	def set_unique_taxon_name(self, s: str):
		self._uniq_name = s
		return

	@property
	def is_empty_taxon(self) -> bool:
		"""
		return true if taxon_name is not set
		"""
		return not self.taxon_name

	def __bool__(self):
		return not self.is_empty_taxon

	@classmethod
	def from_str(cls, s: str):
		"""
		parse taxonomy from the mothur output format: taxonomy_str(bootstrap)

		example:
		Gammaproteobacteria(100)
		"""
		if not s:
			new = cls("", 0)
		else:
			m = cls.TAX_REGEX.match(s)
			if not m:
				raise ValueError("taxon string in invalid format: '%s'" % s)
			new = cls(taxon_name = m.group(1), bootstrap = int(m.group(2)))
		return new

	def to_str(self):
		return "%s(%u)" % (self.taxon_name, self.bootstrap)

	def __str__(self):
		return self.to_str()


class MothurLineage(list):
	"""
	MothurLineage is a list of MothurTaxonomy's, in ranked order:
	kingdom, phylum, class, order, family, genus, species
	"""
	RANK = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
	RANK_ID = {v: i for (i, v) in enumerate(RANK)}

	def _get_by_id(self, rank_id: int) -> "MothurTaxonomy or None":
		"""
		return the taxon if is defined at given rank (kingdom = 0, species = 6)
		otherwise return None
		"""
		return super().__getitem__(rank_id) if rank_id < self.depth else None

	def _get_by_name(self, rank_name: str) -> "MothurTaxonomy or None":
		if rank_name not in self.RANK_ID:
			raise ValueError("invalid rank name: '%s'" % rank_name)
		return self._get_by_id(self.RANK_ID[rank_name])

	def __getitem__(self, arg):
		return self.get_taxon_by_rank(arg)

	def get_taxon_by_rank(self, rank) -> MothurTaxonomy:
		if isinstance(rank, int):
			ret = self._get_by_id(rank)
		elif isinstance(rank, str):
			ret = self._get_by_name(rank)
		else:
			raise IndexError("rank must be int or str, not '%s'"\
				% type(rank).__name__)
		return ret

	@property
	def depth(self) -> int:
		return len(self)

	def set_unique_taxon_names(self, *, _from = 0, _last_known: str = None):
		"""
		set all unique names for underlying taxonomies in this lineage
		the major concern is the label 'uncultured'
		"""
		taxon = self.get_taxon_by_rank(_from)
		if not taxon:
			return # recursion endpoint, should get a None here
		if taxon.taxon_name.startswith("uncultured"):
			uniq_name = (_last_known + "_" + self.RANK[_from]) if _last_known\
				else "unknown"
			taxon.set_unique_taxon_name(uniq_name)
		else:
			_last_known = taxon.taxon_name
		return self.set_unique_taxon_names(_from = _from + 1,
			_last_known = _last_known)

	@classmethod
	def from_str(cls, s: str):
		"""
		parse lineage from mothur output format;
		the lineage is a string of standard taxonomies separated by semi-colon:

		example:
		Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);

		note there is always a final semi-colon at the end
		"""
		sp = s.rstrip(";").split(";")
		new = cls([MothurTaxonomy.from_str(s) for s in sp])
		new.set_unique_taxon_names()
		return new

	def to_str(self):
		return (";").join([str(t) for t in self]) + ";"

	def __str__(self):
		return self.to_str()


class MothurOtuTaxonomy(object):
	"""
	MothurOtuTaxonomy targets at the i/o of each line in mothur taxonomy output;
	the format is 3-column tab-delimited:
	otu_name, size, taxonomy (lineage)

	example:
	Otu000001	412495	Bacteria(100);Bacteroidetes(100);Sphingobacteriia(100);Sphingobacteriales(100);PHOS-HE51(100);PHOS-HE51_ge(100);
	"""
	def __init__(self, otu_name: str, size: int, taxonomy: MothurLineage,
			*ka, **kw):
		super().__init__(*ka, **kw)
		self.otu_name	= otu_name
		self.size		= size
		self.taxonomy	= taxonomy
		return

	@classmethod
	def from_str(cls, s: str):
		sp = s.split("\t")
		return cls(sp[0], int(sp[1]), MothurLineage.from_str(sp[2]))

	def to_str(self):
		return ("\t").join([self.otu_name, str(self.size), str(self.taxonomy)])

	def __str__(self):
		return self.to_str()


class MothurOtuTaxonomyDict(dict):
	"""
	MothurOtuTaxonomyDict is a wrapper to read mothur taxonomy output file;
	dict maps otu_name to taxonomy object
	"""
	@classmethod
	def from_file(cls, file, *ka, **kw):
		"""
		parse mothur taxonomy output file by lines
		note the file contains a header needs to be discarded
		"""
		new = cls(*ka, **kw)
		with file_util.get_fp(file, "r") as fp:
			for line in fp:
				line = line.rstrip()
				if line == "OTU\tSize\tTaxonomy":
					continue # skip header line
				else:
					otu = MothurOtuTaxonomy.from_str(line)
					new[otu.otu_name] = otu
		return new

