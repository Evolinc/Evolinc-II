#!/usr/bin/perl
#Author: Andrew Nelson; andrew.d.l.nelson@gmail.com
#Usage: perl <Batch_mlocarna.pl> <listfile>
#Script to create structural alignments in batch mode using mlocarna package (http://www.bioinf.uni-freiburg.de/Software/LocARNA/)
use strict;
use warnings;

my $listFile = $ARGV[0];
my @list;

open (AFILE, $listFile) or die "cannot open $listFile\n";
while (my $line = <AFILE>) {
        chomp $line;
        push @list, $line;
}
close AFILE;
#print "\n@list\n\n"; #to test the elements of the array

for (my $i=0; $i<@list; $i++) {
        my $file = $list[$i];
        system("echo 'Running mlocarna on $file'");
        system("mlocarna --probabilistic --threads=60 -q $file");
        system("mv $file.out Structures_from_MSA/$file.out");
}
