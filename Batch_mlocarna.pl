#!/usr/bin/perl
#Author: Andrew Nelson; andrew.d.l.nelson@gmail.com
#Usage: perl <Batch_mlocarna.pl> <listfile>
#Script to create structural alignments in batch mode using mlocarna package (http://www.bioinf.uni-freiburg.de/Software/LocARNA/)
use strict;
use warnings;

my $listFile = $ARGV[0];
my $threads = $ARGV[1];
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
        system("perl /singleline.pl $file.fasta >$file.sl.fasta");
        system("grep -A 1 'TBH' $file.sl.fasta >$file.fasta");
        system("rm $file.sl.fasta");
        system("mlocarna --threads=$threads -q $file.fasta");
        system("ps2pdf $file.out/results/alirna.ps $file.out/results/alirna.pdf");
        system("ps2pdf $file.out/results/aln.ps $file.out/results/aln.pdf");
        system("mv $file.out Structures_from_MSA/$file.out");
}
