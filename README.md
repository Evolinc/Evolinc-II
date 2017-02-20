# EVOLINC-II v1.0
Evolinc-II pipeline is designed to perform a series of comparative genomic and transcriptomic analyses across an evolutionary timescale of the user’s choosing and on any number (1-1000s) of query lincRNAs. 

Evolinc-II minimally requires the following input data

1. A FASTA file of lincRNA sequences of all genomes to be interrogated
2. A single column text file with all species listed in order of phylogenetic relatedness to the query species 
3. Evolinc-II can optionally incorporate genome annotation files (GTF) and known lincRNA datasets from target species in FASTA format. 

 
# Availability
### Using Docker image

Since there are several dependencies (can be seen in Dockerfile) to Evolinc-II to make it run on your linux or MAC OS, we highly recommend to use Docker image for [Evolinc-II](https://hub.docker.com/r/cyverse/evolinc-ii/) or use the [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-ii/~/dockerfile/) to build an image and then use the built image for running Evolinc-II. Use the [sample data](https://github.com/Evolinc/Evolinc-II/releases/download/v1.0/sample_data.zip) to run Evolinc-II.

```
docker pull cyverse/evolinc-ii:1.0
unzip https://github.com/Evolinc/Evolinc-II/releases/download/v1.0/sample_data.zip
docker run --rm -v $(pwd):/working-dir -w /working-dir docker.io/cyverse/evolinc-ii:1.0 -b sample_data/Blasting_list.txt -l sample_data/Evolinc-II Sample_query_lincRNA_data_set_for_webinar.fasta -q Atha -i sample_data -s sample_data/test_species_list.txt -o test_out_Evolinc_II -v 1e-20
```

### Using CyVerse Discovery Environment

The [Evolinc-II app](https://de.cyverse.org/de/?type=apps&app-id=3ef009a2-7b99-11e6-a1d6-008cfa5ae621&system-id=de) is currently integrated in CyVerse’s Discovery Environment and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://wiki.cyverse.org/wiki/display/TUT/Evolinc+in+the+Discovery+Environment). 

# Issues
If you experience any issues with running this Docker image, please contact *upendradevisetty at goooglemail.com* 

# Copyright free
The sources in this github repository are copyright free. Thus you are allowed to use these sources in which ever way you like. Please be aware that other license terms apply for the `Notung.jar` used inside the image. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license. 
