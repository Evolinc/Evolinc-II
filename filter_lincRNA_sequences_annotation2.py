#!/usr/bin/env python

import sys

infile = sys.argv[1] # input file for minimap2 SAM output
outfile = sys.argv[2] # output file for filtered SAM output

final = list()
with open(infile, 'rU' ) as fh_in:
	for line in fh_in:
		line = line.strip().split()
	ID = "%s\t%s" % (line[2], line[0])
	empty = "*"
	test = line[2]
	if test == empty:
		pass
	else:
		final.append(ID)

with open(outfile, "w") as fh_out:
	final3 = "\n".join(list(final))
	fh_out.write(final3)
	fh_out.write("\n")
