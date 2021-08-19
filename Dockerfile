#AUTHOR Lukas Becker
# miniconda - django - snakemake - NCBI BLAST+ 2.11.0+ - edirect - Dockerfile
#

#base image; maybe choose another image
FROM ubuntu:focal

# Download and install required software
RUN apt-get update -y && apt-get upgrade -y && apt-get install curl -y && apt-get install wget bzip2 -y

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
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz | tar -zxvpf-

# Download & install NCBI EDIRECT
RUN curl -s ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz | \
 tar xzf - && \
 cpanm HTML::Entities && \
 edirect/setup.sh

# Update path environment variable
ENV PATH="/blast/edirect:$PATH"

# Update path environment variable
ENV PATH /blast/ncbi-blast-2.11.0+/bin:$PATH

# Update path environment variable
ENV PATH /blast/miniconda3/bin:$PATH

# Update miniconda
RUN conda update conda
RUN conda update --all

# Add bioconda as channel
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Download and install snakemake via the bioconda channel; version: 5.32.2
RUN conda install -c bioconda snakemake=5.32.2

# Set-up django and postgres packages


#copy django project into container
RUN mkdir /blast/reciprocal_blast
COPY /celery_blast /blast/reciprocal_blast

#create default BLASTDB environment variable
#RUN mkdir /blast/nr_database
ENV BLASTDB /blast/reciprocal_blast/media/databases
#set appropriate working directory in order to allow development with docker
WORKDIR /blast/reciprocal_blast

#wait-for script to make the web-app and celery worker waiting for the services to finish
RUN apt-get update -y && apt-get install -y netcat
RUN wget -qO- https://raw.githubusercontent.com/eficode/wait-for/v2.1.0/wait-for > wait-for
RUN chmod u+x wait-for

COPY requirements.txt /blast/reciprocal_blast
RUN pip install -r requirements.txt
EXPOSE 8080
# Delete not required packages etc..
RUN apt-get autoremove --purge --yes && apt-get clean && rm -rf /var/lib/apt/lists/*
# Optional commands e.g. initiating scripts
CMD ["bash"]
