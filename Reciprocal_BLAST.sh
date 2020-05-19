#!/bin/bash
#author: Andrew Nelson; andrew.d.l.nelson@gmail.com
#This script is designed to perform reciprocal BLAST of identified lincRNA orthologs in conjunction with Evolinc II
# Usage: 
# sh Reciprocal_BLAST.sh 

usage() {
      echo ""
      echo "Usage : sh $0 -g genome -s ortho_sequences -f ortho_gff -a query_gff -b query_sp -c query_genome -e evalue -n threads"
      echo ""

cat <<'EOF'
     -g    </path/to/genome in FASTA format>

     -s    </path/to/putative ortholog sequences in FASTA format>

     -f    </path/to/gff file of ortholog sequences>

     -a    <query gff in four letter format. ie., Atha for A. thaliana>

     -b    <query species in four letter format>

     -c    <query genome in fasta format>

     -e    <e-value>

     -n     <Number of threads>

     -h    Show this usage information
EOF
    exit 0
}


while getopts ":g:s:f:a:hb:c:e:n:" opt; do
  case $opt in
    g)
     Genome=$OPTARG
     ;;
    s)
     Put_ortholog=$OPTARG
     ;;
    f)
     Put_ortholog_gff=$OPTARG
     ;;
    a)
     query_gff=$OPTARG
     ;;
    b)
     query_species=$OPTARG
     ;;
    c)
     query_genome=$OPTARG
     ;;
    e)
     value=$OPTARG # e-value
     ;;
    n)
     threads=$OPTARG
     ;;  
    h) 
     usage 
     exit 1
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

echo "Starting on $Put_ortholog"
mkdir -p Reciprocal_BLAST_DB
mkdir -p Reciprocal_BLAST_Return

###Switching over to minimap for the reciprocal search
minimap2 -t $threads -a -w5 --splice -G5k -A2 -B8 -O12,32 -E1,0 -L -o Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sam $query_genome $Put_ortholog

#convert the alignment to bed format, sorted so that the higher quality alignment is first for each query
samtools sort -@ $threads Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sam -o Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.sam
samtools view -@ $threads -S -b Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.sam -o Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.bam 
bamToBed -i Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.bam > Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.bed 
sort -b -k 4,4 -k 5nr,5 Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.sorted.bed >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.bed
sed -i 's~_TBH_1~~g' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.bed

#assign subject_species argument
subject_species=${Put_ortholog:0:4}

# Convert minimap2 result to gff
Rscript /bed_to_gff.R Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.bed Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal_exon.out.gff $subject_species
sed -i 's~sequence_feature~exon~g' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal_exon.out.gff
grep -v "#" Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal_exon.out.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.gff
sed -i 's~name=~~g' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.gff && awk '{print $2"_"$9 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "gene_id " $2 "_" $9 "; " "transcript_id " $2 "_" $9 ";" "\t" "500"}' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.out.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff


# Pick the top hit and merge close hits
python /merge_close_hits.py Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.merged.gff

#Generate a gff file with all hit fragments, instead of merging the whole locus together. This will retain information about regions of high conservation instead of assuming the entire region is conserved
python /merge_close_hits_2.py Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.gff Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.fragments_not_merged.gff

# Change the file format since the gffread requires the chromosome to be first column
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.merged.gff > temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.merged.gff

grep "TBH" Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.merged.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.gff
bedtools sort -i Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.sorted.gff
cut -f 2 Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.sorted.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.almost.txt
grep "TBH" $query_gff >Reciprocal_BLAST_Return/$query_gff.TBH.gff
bedtools sort -i Reciprocal_BLAST_Return/$query_gff.TBH.gff >Reciprocal_BLAST_Return/$query_gff.TBH.sorted.gff
bedtools closest -a Reciprocal_BLAST_Return/$query_gff.TBH.sorted.gff -b Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.sorted.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.gff
cut -f 1,2,4,5,10,11,13,14 Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff

awk -F'\t' '{gsub(/_Known.*/,"",$6); print}' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff > temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff

###
awk -F'\t' 'substr($2,6) == substr($6,11) {print $2}' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.filtered.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.IDs.list.txt

uniq Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.IDs.list.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt

#The line below is being used to generate a list of families (and all their paralogous sequences) to be pulled from the overall FASTA file of sequence homologs in each species.
sed 's~_TBH_1~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt

sed -i 's~^.....~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt
sed -i 's~^~>'$subject_species'_~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt

sleep 5
sed 's~^.....~~g' Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.list.txt > Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.stripped.list.txt
grep -Ef Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.stripped.list.txt Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.almost.txt >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.TBH.filtered.gff
sort -u Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.comparison.TBH.filtered.gff >Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.for.renaming.final.txt

# Prepare the fragment gff file for output
awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 " " $10  " " $11 " " $12}' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.fragments_not_merged.gff > temp && mv temp Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.fragments_not_merged.gff
grep "TBH" Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.stripped.out.fragments_not_merged.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged.gff
grep -Ef Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.for.renaming.final.txt Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged_reciprocal.gff
bedtools sort -i Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged_reciprocal.gff >Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged_reciprocal_sorted.gff
# get rid of TBH nomenclature
sed -i 's~_TBH_.~~g' Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged_reciprocal_sorted.gff
# Move the fragment gff file up one directory so that it doesn't get deleted
cp Reciprocal_BLAST_Return/$Put_ortholog.Reciprocal.TBH.out.fragments_not_merged_reciprocal_sorted.gff .

#find ./Reciprocal_BLAST_Return -type f ! -name '*final*' -delete
python /assign_annotation_ortholog.py ../Orthologs/$Put_ortholog Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.Unique.IDs.for.renaming.final.txt ../Orthologs/$Put_ortholog.orthologs.identified.fasta
rm -f ../Orthologs/$Put_ortholog
perl /singleline.pl ../Orthologs/$Put_ortholog.orthologs.identified.fasta >../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta
grep -A 1 -f Reciprocal_BLAST_Return/$Put_ortholog.reciprocal.final.list.for.pooling.families.txt ../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta >../Orthologs/$Put_ortholog.orthologs.only.fasta

#Clean things up a bit
rm -f ../Orthologs/$Put_ortholog.orthologs.identified.fasta
rm -f ../Orthologs/$Put_ortholog.orthologs.identified.singleline.fasta
#find ./Reciprocal_BLAST_Return -type f ! -name '*renaming*' -delete
echo "Finished with $Put_ortholog"
