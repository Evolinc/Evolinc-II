#!/bin/bash

apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		python \
		perl \
		wget \
		curl \
		python-matplotlib \
		python-numpy \
        python-pandas \
        openjdk-8-jdk

# Bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar zxvf bedtools-2.26.0.tar.gz
cd bedtools2 && make
cd ..

# Cufflinks
wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Mafft
apt-get install -y mafft

# cpan
apt-get install -y cpanminus

# Install BioPerl dependancies, mostly from cpan
apt-get install --yes \
 libpixman-1-0 \
 libpixman-1-dev \
 graphviz \
 libxml-parser-perl \
 libsoap-lite-perl 

cpanm Test::Most \
 Algorithm::Munkres \
 Array::Compare Clone \
 PostScript::TextBlock \
 SVG \
 SVG::Graph \
 Set::Scalar \
 Sort::Naturally \
 Graph \
 GraphViz \
 HTML::TableExtract \
 Convert::Binary::C \
 Math::Random \
 Error \
 Spreadsheet::ParseExcel \
 XML::Parser::PerlSAX \
 XML::SAX::Writer \
 XML::Twig XML::Writer

apt-get install -y \
 libxml-libxml-perl \
 libxml-dom-xpath-perl \
 libxml-libxml-simple-perl \
 libxml-dom-perl

# Install BioPerl last built
cpanm -v  \
 CJFIELDS/BioPerl-1.6.924.tar.gz 

# Biopython
curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
python get-pip.py
pip install biopython

# R libraries
echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
apt-get update
apt-get install -y r-base r-base-dev
Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'
Rscript -e 'install.packages("reshape2", dependencies = TRUE, repos="http://cran.rstudio.com/");'
Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# RAxML
git clone https://github.com/stamatak/standard-RAxML.git
WORKDIR /standard-RAxML
make -f Makefile.SSE3.PTHREADS.gcc
cp raxmlHPC-PTHREADS-SSE3 /usr/bin/
WORKDIR /

# NCBI
wget -O- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zxvf -

# Setting paths to all the softwares
#ENV BINPATH /usr/bin
#ENV PATH /bedtools2/bin/:$PATH
#ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH
#ENV PATH /ncbi-blast-2.6.0+/bin/:$PATH

# Add all the scripts to the root directory Path
chmod +x /Building_Families.sh
chmod +x /evolinc-part-II.sh && cp /evolinc-part-II.sh $BINPATH
