#!/bin/bash
#author: Andrew Nelson; andrew.d.l.nelson@gmail.com
# Script to process cuffcompare output file to generate lincRNA
# Usage: 
# sh Building_Families.sh  -g subjectgenome.fa -s subject_species -q query_species -l lincRNAs.fa -e subject_gff -k known_lincRNAs -v e_value

while getopts ":l:q:s:f:k:e:g:n:" opt; do
  case $opt in
    l)
      lincRNAfasta=$OPTARG
    ;;
    q)
      query_species=$OPTARG
      ;;
    s)
      subject_species=$OPTARG
      ;;
    f)
      subject_gff=$OPTARG
      ;;
    g)
      subject_genome=$OPTARG
      ;;
    k)
      known_lincRNAs=$OPTARG
      ;;
    e)
      value=$OPTARG
      ;;
    n)
      threads=$OPTARG
      ;;  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Building_Families.sh script is reporting this"
echo "Starting! on $subject_species"

###This step is to make sure that query lncRNAs do not contain underscores or equal signs The if/else statement is more to announce whether names look good.
notclean=1
cleanquery=$(grep -c -E '_| ' $lincRNAfasta)
if [ "$cleanquery" -ge "$notclean" ]; then
	echo "Query lincRNA names are messy, attempting to clean them up. Replacing spaces and underscores with periods."
	sed -i 's~_~.~g; s~gene=~~g; s~ ~.~g; s~\.\.~\.~g;' $lincRNAfasta
else 
	echo "Query lincRNA names look good, proceeding"
fi

### Check to see if query FASTA file contains only one lincRNA
count=$(grep -c ">" $lincRNAfasta)
param=2
if [ "$count" -ge "$param" ]; then
   touch $lincRNAfasta
else
   cat $lincRNAfasta $lincRNAfasta >temp && mv temp $lincRNAfasta
   sed -i 's~.>~\n>~g' $lincRNAfasta
fi

#Moving the aligner to Minimap due to speed issues with BLASTn. The algorithm used favors longer read alignments and is used to align nanopore or pacbio reads to genomes.
minimap2 -t $threads -a -w5 --splice -G5k -A2 -B8 -O12,32 -E1,0 -L -o Homology_Search/$subject_species.sam $subject_genome $lincRNAfasta

#convert the alignment to bed format, sorted so that the higher quality alignment is first for each query
samtools sort -@ $threads Homology_Search/$subject_species.sam -o Homology_Search/$subject_species.sorted.sam
samtools view -@ $threads -S -b Homology_Search/$subject_species.sorted.sam -o Homology_Search/$subject_species.sorted.bam 
bamToBed -i Homology_Search/$subject_species.sorted.bam > Homology_Search/$subject_species.sorted.bed 
sort -b -k 4,4 -k 5nr,5 Homology_Search/$subject_species.sorted.bed >Homology_Search/$subject_species.bed

#Convert from bed to gff

Rscript /bed_to_gff.R Homology_Search/$subject_species.bed Homology_Search/$subject_species_exon.gff $subject_species
sed -i 's~sequence_feature~exon~g' Homology_Search/$subject_species_exon.gff
grep -v "#" Homology_Search/$subject_species_exon.gff >Homology_Search/$subject_species.gff
sed -i 's~name=~~g' Homology_Search/$subject_species.gff && awk '{print $2"_"$9 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "gene_id " $2 "_" $9 "; " "transcript_id " $2 "_" $9 ";" "\t" "500"}' Homology_Search/$subject_species.gff >Homology_Search/$subject_species.out.gff

# Merge close hits in the gff file and add TBH to the top ranked return. Each of the other returns get a unique numerical ID appended.
python /merge_close_hits.py Homology_Search/$subject_species.out.gff Homology_Search/$subject_species.out.merged.gff

grep "TBH" Homology_Search/$subject_species.out.merged.gff >Homology_Search/$subject_species.out.TBH.only.gff

#The below lines grab five of the top non-TBH hits (in case of large, multi-hit families we don't want to grab too many). This is informative for examining duplication events.
grep -v "TBH" Homology_Search/$subject_species.out.merged.gff >Homology_Search/$subject_species.out.no_TBH_precursor.only.gff
cut -f 1 Homology_Search/$subject_species.out.no_TBH_precursor.only.gff | grep '._1$\|._2$\|._3$\|._4$\|._5$' | sed 's~$~\t~g' >Homology_Search/$subject_species.out.no_TBH.only.gff.list.txt
grep -f Homology_Search/$subject_species.out.no_TBH.only.gff.list.txt Homology_Search/$subject_species.out.no_TBH_precursor.only.gff >Homology_Search/$subject_species.out.no_TBH.only.gff

# Convert gtf to fasta
# Change the file format since the gffread requires the chromosome to be first column
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.merged.gff > temp && mv temp Homology_Search/$subject_species.out.merged.gff
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.TBH.only.gff > temp && mv temp Homology_Search/$subject_species.out.TBH.only.gff
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Homology_Search/$subject_species.out.no_TBH.only.gff > temp && mv temp Homology_Search/$subject_species.out.no_TBH.only.gff
gffread Homology_Search/$subject_species.out.merged.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.orthologs.fasta
gffread Homology_Search/$subject_species.out.TBH.only.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta
gffread Homology_Search/$subject_species.out.no_TBH.only.gff -g $subject_genome -w Homology_Search/$subject_species.$query_species.no_TBH.orthologs.fasta

# Optional argument - 1 (transcriptome) 
# Bedtools can check for overlap in either strand. To work with downstream scripts, we can run bedtools twice, one checking for sense and making a list, and the other checking for antisense #(and making a list).
if [ ! -z $subject_gff ]; then 
	echo "Testing for overlap with the transcriptome of" $subject_species
    # sense
	intersectBed -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -u -s > Homology_Search/$subject_species.annotation.sense.gff &&
    #the -f of 0.1 indicates a 10% overlap. This may be too much or not enough. We will have to check emperically. -u reports the mere presence of any unique intersecting gene.
    # antisense
    # This is for capturing the GeneID of the overlapping gene
	intersectBed -wb -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -s > Homology_Search/$subject_species.mod.annotation.sense.gff &&
	awk '{print $2}' Homology_Search/$subject_species.annotation.sense.gff > Homology_Search/$subject_species.annotation.sense.list.txt &&
    # Retain only the second column (gene IDs) and save as a list
	intersectBed -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -u -S > Homology_Search/$subject_species.annotation.antisense.gff &&
    # This is the same as above, except for -S searches only for those overlapping candidates that are on the antisense strand.
    # This is for capturing the GeneID of the overlapping gene
	intersectBed -wb -a Homology_Search/$subject_species.out.TBH.only.gff -b $subject_gff -f 0.1 -S > Homology_Search/$subject_species.mod.annotation.antisense.gff &&
	awk '{print $2}' Homology_Search/$subject_species.annotation.antisense.gff > Homology_Search/$subject_species.annotation.antisense.list.txt &&
    # Retain only the second column (gene IDs) and save as a list
    # assign annotation
    # sense
	python /assign_sense_annotation.py Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta Homology_Search/$subject_species.annotation.sense.list.txt Homology_Search/$subject_species.$query_species.orthologs.sense.renamed.fasta &&
    # assign_sense_annotation.py works as .py file_with_sequences_to_rename.fasta list_of_headers_to_rename >renamed_output.fasta
    #sense + antisense
	python /assign_antisense_annotation.py Homology_Search/$subject_species.$query_species.orthologs.sense.renamed.fasta Homology_Search/$subject_species.annotation.antisense.list.txt Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta
else
	mv Homology_Search/$subject_species.$query_species.TBH.orthologs.fasta Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta
fi

# Optional argument - 2 (Known lincRNAs)
if [ ! -z $known_lincRNAs ]; then 
	echo "Testing for similarity with" $known_lincRNAs
	cleanknown=$(grep -c -E '_| ' $known_lincRNAs)
	if [ "$cleanknown" -ge "$notclean" ]; then
	echo "Known subject lincRNA names are messy, attempting to clean them up. Replacing spaces and underscores with periods."
	sed -i 's~_~.~g; s~gene=~~g; s~ ~.~g; s~\.\.~\.~g;' $known_lincRNAs
	else 
	echo "Known subject lincRNA names look good, proceeding"
	fi
    # Comparing to known lincRNAs with Minimap2. Run minimap2 then clean up the file
	minimap2 -t $threads -a -w5 --splice -G5k -A2 -B8 -O12,32 -E1,0 -L -o Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.sam Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta $known_lincRNAs
	grep -v "SQ" Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.sam | grep -v "@PG" > Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.sam.tested.out
	grep "TBH" Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.sam.tested.out | cut -f 3 | sort -u >  Homology_Search/$subject_species.lincRNA_annotation.list.txt
	
	# If list file is empty, remove sam.tested.out file so that it doesn't cause problems towards the end of the evolinc-ii.sh
	if [ -s Homology_Search/$subject_species.lincRNA_annotation.list.txt ]; then
	echo "Overlap found with known lincRNAs"
	else
	echo "no overlap found with known lincRNAs"
        rm Homology_Search/$subject_species.$query_species.orthologs.renamed.lincRNAs_tested.sam.tested.out
        fi
	
	# Assign the annotation of lincRNA to the known lincRNA
	python /assign_annotation_lincRNA.py Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta Homology_Search/$subject_species.lincRNA_annotation.list.txt Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta
else
	mv Homology_Search/$subject_species.$query_species.orthologs.renamed.fasta Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta
fi

# Remove duplicate sequences
perl /Remove_dup_seqs.pl Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta
### This is just to keep continuity with below. I removed a section below and I want to make sure everything continues to work before I rename files.
mv Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta

### Where the old reiterative BLAST section used to start. Not needed anymore but have retained a few steps for file name continuity
grep ">" Homology_Search/$subject_species.$query_species.orthologs.lincRNA_tested.renamed.fasta.dup_removed.fasta |sed 's~^......~~g'| sed 's~_Known_lincRNA~~g'| sed 's~_Known_Gene_Sense~~g' | sed 's~_Known_Gene_Antisense~~g' | sed 's~_TBH_1~~g' | sed 's~_.$~~g' | sed 's~_..$~~g' | sed 's~_...$~~g' | sed 's~>~~g' | sort -u > Homology_Search/List_of_identified_putative_orthologs.txt
grep ">" $lincRNAfasta > Homology_Search/List_of_all_query_lincRNAs.txt
# Search for lincRNA IDs that have not been found with the previous BLAST search. If none are found, still populate the file with the title "List_of_non_identified_lincRNAs" so the next step doesn't fail.
grep -vwf Homology_Search/List_of_identified_putative_orthologs.txt Homology_Search/List_of_all_query_lincRNAs.txt > Homology_Search/List_of_non_identified_query_lincRNAs.txt
sed -i 's~>~~g' Homology_Search/List_of_non_identified_query_lincRNAs.txt
perl /singleline.pl $lincRNAfasta > Homology_Search/query_lincRNAs_singleline.fasta
grep -A 1 -f Homology_Search/List_of_non_identified_query_lincRNAs.txt Homology_Search/query_lincRNAs_singleline.fasta | sed 's~--~~g' > Homology_Search/$query_species.$subject_species.non_identified_from_first_round.fasta

mv Homology_Search/$subject_species.out.TBH.only.gff Homology_Search/$subject_species.200nt_plus_full_length_searches.gff
sort -k 1,1 -k 4,4n Homology_Search/$subject_species.200nt_plus_full_length_searches.gff > Homology_Search/$subject_species.200nt_plus_full_length_searches_sorted.gff


cat Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta Homology_Search/$subject_species.$query_species.no_TBH.orthologs.fasta > Homology_Search/$subject_species.$query_species.200nt_plus_all_put_homologs.fasta

# Move genome and putative ortholog files for reciprocal BLAST to reciprocal BLAST folder. Also, rename putative orthologs file
cp $subject_genome Reciprocal_BLAST/
cp Homology_Search/$subject_species.$query_species.200nt_plus_full_length_searches.fasta Reciprocal_BLAST/$subject_species.$query_species.putative_orthologs.fasta
cp Homology_Search/$subject_species.200nt_plus_full_length_searches_sorted.gff Reciprocal_BLAST/$subject_species.$query_species.coords.gff
cp Homology_Search/$subject_species.$query_species.200nt_plus_all_put_homologs.fasta Orthologs/$subject_species.$query_species.putative_orthologs.fasta
echo "Finished with" $subject_species
