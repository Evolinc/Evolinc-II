[![Docker Pulls](https://img.shields.io/docker/pulls/evolinc/evolinc-ii.svg)](https://hub.docker.com/r/evolinc/evolinc-ii/)
[![Docker Stars](https://img.shields.io/docker/stars/evolinc/evolinc-ii.svg)](https://hub.docker.com/r/evolinc/evolinc-ii/)
[![Docker Build Status](https://img.shields.io/docker/build/evolinc/evolinc-ii.svg)](https://hub.docker.com/r/evolinc/evolinc-ii/)
[![Docker Automated build](https://img.shields.io/docker/automated/evolinc/evolinc-ii.svg)](https://hub.docker.com/r/evolinc/evolinc-ii/)
[![Release](https://shields.beevelop.com/github/release/Evolinc/Evolinc-II.svg?style=flat-square)](https://github.com/Evolinc/Evolinc-II/releases)

# EVOLINC-II 2.0
Evolinc-II pipeline is designed to perform a series of comparative genomic and transcriptomic analyses across an evolutionary timescale of the user’s choosing and on any number (1-1000s) of query lincRNAs. Please see our chapter in Plant Comparative Genomics for a detailed explanation on how to use Evolinc-II in your research (https://link.springer.com/protocol/10.1007/978-1-0716-2429-6_3). 

Evolinc-II minimally requires the following input data:

1. A FASTA file of lincRNA sequences of all genomes to be interrogated.
2. A single column text file with all species listed in order of phylogenetic relatedness to the query species. See Discussion_on_species_list.txt for more info.
3. Evolinc-II can optionally incorporate genome annotation files (GTF) and known lincRNA datasets from target species in FASTA format.
4. Currently, Evolinc-II iterates through all of its analyses using a tab-delimited file that contains all of this information, including file names and species abbreviations. See [Discussion_on_BLASTing_file.md](https://github.com/Evolinc/Evolinc-II/blob/master/Discussion_on_BLASTing_file.md) for more info.

 
# Availability
### Using a Docker image

Since there are several dependencies (can be seen in [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-ii/~/dockerfile/)) for running Evolinc-II on your linux or MAC OS, we highly recommend using the Docker image for [Evolinc-II](https://hub.docker.com/r/cyverse/evolinc-ii/) or use the [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-ii/~/dockerfile/) to build an image and then use the built image.

```
# Pull the image from CyVerse Dockerhub
docker pull evolinc/evolinc-ii:2.0

# See the command line help for the image
docker run evolinc/evolinc-ii:2.0 -h

# Download the test data
wget https://github.com/Evolinc/Evolinc-II/releases/download/v1.0/sample_data.zip
unzip sample_data.zip

# Run Evolinc-II on the test data
docker run --rm -v $(pwd):/working-dir -w /working-dir evolinc/evolinc-ii:2.0 -b sample_data/Blasting_list.txt -l sample_data/Sample_query_lincRNA_data_set_for_webinar.fasta -q Atha -i sample_data -s sample_data/test_species_list.txt -o test_out -v 1e-20 -n 8
```

### Using CyVerse Discovery Environment

The [Evolinc-II app](https://de.cyverse.org/de/?type=apps&app-id=3ef009a2-7b99-11e6-a1d6-008cfa5ae621&system-id=de) is currently integrated in CyVerse’s Discovery Environment and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://cyverse.atlassian.net/wiki/spaces/DEapps/pages/241882272/Evolinc+in+the+Discovery+Environment). If you do not have access to a high performance computing cluster we highly recommend using Evolinc-II within the DE. Limited computing resources can lead to extremely long Evolinc-II run times, especially when searching through large genomes. Not only is the CyVerse DE free to researchers, but many of the genomes are already integrated (and almost all others are available through CoGe; genomevolution.org), thereby alleviating the data storage hurdle common with these types of analyses. 

# Issues
If you experience any issues with running Evolinc-I (DE app or source code or Docker image), please open an issue on this github repo.
 

# Copyright free
The sources in this github repository are copyright free. Thus you are allowed to use these sources in which ever way you like. Please be aware that other license terms apply for the `Notung.jar` used inside the image. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license. 

# Citing Evolinc-II
If you have used Evolinc in your research, please cite the below manuscripts.

*Andrew D. Nelson&ast;, Upendra K. Devisetty&ast;, Kyle Palos, Asher K. Haug-Baltzell, Eric Lyons, Mark A. Beilstein (2017). "Evolinc: a comparative transcriptomics and genomics pipeline for quickly identifying sequence conserved lincRNAs for functional analysis". Frontiers in Genetics. 1(10)*

*Anna C. Nelson Dittrich and Andrew D. L. Nelson (2022). "High-Throughput Evolutionary Comparative Analysis of Long Intergenic Noncoding RNAs in Multiple Organisms". Plant Comparative Genomics (Part of the Methods in Molecular Biology Book Series)*
