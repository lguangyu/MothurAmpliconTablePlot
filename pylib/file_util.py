#!/usr/bin/env python3

import io


def get_fp(f, *ka, factory = open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument must be str or file handle, not '%s'"\
			% type(f).__name__)
	return ret
