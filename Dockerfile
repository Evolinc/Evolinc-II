FROM ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for evolinc-ii pipeline"

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		libcurl4-openssl-dev \
		libssl-dev \
		python \
		perl \
		wget \
		curl \
		python-matplotlib \
		python-numpy \
        	python-pandas \
		perl \
		bioperl \
		gnupg2 \
        	openjdk-8-jdk
	
# Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN tar zxvf bedtools-2.26.0.tar.gz
RUN cd bedtools2 && make
RUN cd ..

# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Mafft
RUN apt-get install -y mafft

# Biopython
RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
RUN python get-pip.py
RUN pip install biopython

# R libraries
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update
RUN apt-get install -y r-base r-base-dev
RUN Rscript -e 'install.packages("BiocManager");'
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'install.packages("reshape2", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'BiocManager::install("rtracklayer", dependencies = TRUE);'

# RAxML
RUN git clone https://github.com/stamatak/standard-RAxML.git
WORKDIR /standard-RAxML
RUN make -f Makefile.SSE3.PTHREADS.gcc
RUN cp raxmlHPC-PTHREADS-SSE3 /usr/bin/
WORKDIR /

# NCBI
RUN wget -O- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zxvf -

#minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
RUN tar -jxvf minimap2-2.17_x64-linux.tar.bz2

#samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -jxvf samtools-1.10.tar.bz2
RUN cd samtools-1.10 \
	&& ./configure --prefix=/samtools \
	&& make \
	&& make install


# Setting paths to all the softwares
ENV BINPATH /usr/bin
ENV PATH /minimap2-2.17_x64-linux/:$PATH
ENV PATH /bedtools2/bin/:$PATH
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /ncbi-blast-2.6.0+/bin/:$PATH
ENV PATH /samtools-1.10/:$PATH

# Add all the scripts to the root directory Path
ADD *.py *.pl *.R *.sh *.jar /
RUN chmod +x /Building_Families.sh
RUN chmod +x /evolinc-part-II.sh && cp /evolinc-part-II.sh $BINPATH

ENTRYPOINT ["/evolinc-part-II.sh"]
CMD ["-h"]
