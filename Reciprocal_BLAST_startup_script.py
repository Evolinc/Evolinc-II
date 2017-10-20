#!/usr/bin/env python

from subprocess import call
import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'hb:v:n:', ['inputfile=', 'e_value=', 'threads=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-b', '--inputfile'):
        inputfile_name = arg
    elif opt in ('-v', '--value'):
        value = arg
    elif opt in ('-n', '--threads'):
        threads = arg
    else:
        usage()
        sys.exit(2)

with open(inputfile_name, 'r') as inpt:
    for line in inpt:
        if line == '\n':
            pass
        else:
            line = line.strip()
            line_split = line.split('\t')

            genome = line_split[0]
            sequences = line_split[1]
            gff = line_split[2]
            query_gff = line_split[3]
            query_species = line_split[4]
            query_genome = line_split[5]


            query = "bash /Reciprocal_BLAST.sh -g %s -s %s -f %s -a %s -b %s -c %s -e %s -n %s" % \
                            (genome, sequences, gff, query_gff, query_species, query_genome, value, threads)
            if gff == query_gff:
                pass
            else:
                call(query, shell=True)
