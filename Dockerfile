# SpliceLauncher
# https://github.com/raphaelleman/SpliceLauncher
# Commit f2aa2efec29f173ad1dda4edd9c190e6e423cd3d

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV SPLICELAUNCHER_VERSION f2aa2efec29f173ad1dda4edd9c190e6e423cd3d
ENV STAR_VERSION=2.7.1a
ENV SAMTOOLS_VERSION=1.9
ENV BEDTOOLS_VERSION=2.28.0
ENV SEQTK_VERSION=1.3

ENV SPLICELAUNCHER_DIR=/usr/local/splicelauncher
ENV PATH ${SPLICELAUNCHER_DIR}/scripts/:${SPLICELAUNCHER_DIR}:/usr/local/bin/:$PATH

ENV LC_ALL="en_US.UTF-8"
ENV LANG="en_US.UTF-8"
ENV LANGUAGE="en_US.UTF-8"

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

ENV STAR_PACKAGES gcc g++ make wget zlib1g-dev
ENV SAMTOOLS_PACKAGES wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses5-dev bzip2
ENV BEDTOOLS_PACKAGES wget gcc g++ make  zlib1g-dev
ENV R_PACKAGES libcairo2-dev libxt-dev
ENV PACKAGES zlib1g liblzma5 libncurses5 r-base python libcairo2 unzip

## INSTALL ALL DEPS
RUN apt-get update --fix-missing \
    && apt-get install -y ${PACKAGES} ${SAMTOOLS_PACKAGES} ${STAR_PACKAGES} ${BEDTOOLS_PACKAGES} ${R_PACKAGES}

## INSTALL STAR
RUN cd /usr/local \
    && wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip \
    && unzip ${STAR_VERSION}.zip \
    && cd STAR-${STAR_VERSION}/source \
    && make STARstatic \
    && cp STAR /usr/local/bin \
    && cd /usr/local/ \
    && rm -rf STAR-${STAR_VERSION}* \
    && rm /usr/local/${STAR_VERSION}.zip

## INSTALL SAMTOOLS
RUN cd /usr/local/ \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-${SAMTOOLS_VERSION}

## INSTALL BEDTOOLS
RUN cd /usr/local/ \
    && wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && tar -xzvf bedtools.tar.gz \
    && rm -rf bedtools.tar.gz \
    && cd bedtools2 \
    && make \
    && mv bin/* /usr/local/bin \
    && cd /usr/local \
    && rm -rf /usr/local/bedtools2 

## INSTALL R LIBRARIES
RUN R -e "install.packages('WriteXLS',lib=$RLIBPATH,repo=$CRANREPO)"  \
    && R -e "install.packages('Cairo',lib=$RLIBPATH,repo=$CRANREPO)" 

## INSTALL SEQTK
RUN cd /usr/local \
    && wget https://github.com/lh3/seqtk/archive/v${SEQTK_VERSION}.zip \
    && unzip v${SEQTK_VERSION}.zip \
    && rm  v${SEQTK_VERSION}.zip \
    && cd seqtk-${SEQTK_VERSION} \
    && make \
    && mv seqtk /use/local/bin \
    && cd /usr/local \
    && rm -rf seqtk-${SEQTK_VERSION}

## Clean Build packages
RUN apt-get remove -y ${STAR_PACKAGES} ${SAMTOOLS_PACKAGES} ${BEDTOOLS_PACKAGES} ${R_PACKAGES} \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

## COPY FILES INTO DOCKER
COPY . ${SPLICELAUNCHER_DIR}

ENTRYPOINT ${SPLICELAUNCHER_DIR}/scripts/SpliceLauncher.r
