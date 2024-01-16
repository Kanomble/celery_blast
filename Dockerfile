#AUTHOR Lukas Becker
# miniconda - django - snakemake - NCBI BLAST+ 2.11.0+ - edirect - Dockerfile
#

#base image; maybe choose another image
FROM ubuntu:focal

# Download and install required software
# netcat for the wait-for script / libdw1 for rpsbproc
RUN apt-get update -y && apt-get upgrade -y && apt-get install curl -y && apt-get install wget bzip2 netcat libdw1 -y 

# Software and packages for the E-Direct Tool
RUN apt-get -y -m update && DEBIAN_FRONTEND="noninteractive" apt-get install -y cpanminus libxml-simple-perl libwww-perl libnet-perl build-essential tzdata

#Set correct timezone
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# set working directory
WORKDIR /blast

# Download and install anaconda; version: 4.9.2
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh
#-b Batch mode with no PATH modifications to ~/.bashrc
RUN bash Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p /blast/miniconda3
RUN rm Miniconda3-py38_4.9.2-Linux-x86_64.sh

# Download & install BLAST
# This RUN cmd might fail due to connection problems with NCBI servers, if you face any issues, repeat the build step
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz | tar -zxvpf-

# Download & install NCBI EDIRECT
#RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz | \
# tar xzf - && \
# cpanm HTML::Entities && \
# edirect/setup.sh
#new command for downloading the e-direct software
RUN sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
RUN mv $HOME/edirect .

# Update path environment variable
ENV PATH="/blast/edirect:$PATH"

# Update path environment variable
ENV PATH /blast/ncbi-blast-2.11.0+/bin:$PATH

# Update path environment variable
ENV PATH /blast/miniconda3/bin:$PATH

# Add bioconda as channel
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Update miniconda
RUN conda update conda
RUN conda update --all

# Download and install snakemake via the bioconda channel
RUN conda install -c bioconda snakemake
# Fixing the error: module 'lib' has no attribute 'OpenSSL_add_all_algorithms'
# Downgrade the cryptography module
RUN python3 -m pip install cryptography==38.0.4
RUN conda install -c conda-forge notebook

# Set-up django and postgres packages
# copy django project into container
RUN mkdir /blast/reciprocal_blast
COPY /celery_blast /blast/reciprocal_blast

# Create default BLASTDB environment variable
ENV BLASTDB /blast/reciprocal_blast/media/databases
# Set appropriate working directory
WORKDIR /blast/reciprocal_blast

# Install mafft and fasttree
RUN apt-get update && (apt-get install -t buster-backports -y mafft fasttree || apt-get install -y mafft fasttree)

# wait-for script to make the web-app and celery worker waiting for the services to finish
RUN wget -qO- https://raw.githubusercontent.com/eficode/wait-for/v2.1.0/wait-for > wait-for
RUN chmod u+x wait-for
RUN mkdir /blast/utilities/
RUN mv wait-for /blast/utilities/
ENV PATH /blast/utilities:$PATH

# Download requirements 
COPY requirements.txt /blast/reciprocal_blast
RUN pip install -r requirements.txt

# Download trimai
RUN wget https://github.com/inab/trimal/archive/refs/heads/trimAl.zip
RUN unzip -d /blast/utilities trimAl.zip
RUN rm trimAl.zip
RUN make -C /blast/utilities/trimal-trimAl/source
ENV PATH /blast/utilities/trimal-trimAl/source:$PATH

# Download mview
RUN wget https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/mview-1.67.tar.gz
RUN tar xvzf mview-1.67.tar.gz -C /blast/utilities/
RUN mkdir /blast/utilities/mview
RUN cd /blast/utilities/mview-1.67/ && perl install.pl /blast/utilities/mview -y
RUN cd /blast/reciprocal_blast
RUN rm mview-1.67.tar.gz
ENV PATH /blast/utilities/mview:$PATH

# Download rpsbproc
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz
RUN tar xvzf RpsbProc-x64-linux.tar.gz -C /blast/utilities/
RUN rm RpsbProc-x64-linux.tar.gz
ENV PATH /blast/utilities/RpsbProc-x64-linux:$PATH

RUN mkdir /blast/utilities//data
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O /blast/utilities//data/cddid.tbl.gz && gzip -d /blast/utilities/data/cddid.tbl.gz
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -O /blast/utilities//data/cdtrack.txt
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -O /blast/utilities//data/family_superfamily_links
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -O /blast/utilities//data/cddannot.dat.gz && gzip -d /blast/utilities/data/cddannot.dat.gz
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -O /blast/utilities//data/cddannot_generic.dat.gz && gzip -d /blast/utilities/data/cddannot_generic.dat.gz
RUN wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -O /blast/utilities//data/bitscore_specific.txt

# Delete not required packages etc..
RUN apt-get autoremove --purge --yes && apt-get clean && rm -rf /var/lib/apt/lists/*

# Optional commands e.g. initiating scripts
CMD ["bash"]
