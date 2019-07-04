FROM ubuntu:16.04

LABEL maintainer "luispedro@big-data-biology.org"

RUN apt update && \
    apt install -y \
        bc \
        build-essential \
        cabal-install \
        coreutils \
        curl \
        git \
        libbz2-dev \
        libcairo2-dev \
        libgif-dev \
        libjpeg-dev \
        liblzma-dev \
        libpthread-stubs0-dev \
        m4 \
        openjdk-8-jdk \
        parallel \
        pkg-config \
        python-dev \
        python-pip \
        software-properties-common \
        unzip \
        wget \
        xz-utils \
        zlib1g \
        zlib1g-dev


# Install pandaseq
RUN apt-add-repository ppa:neufeldlab/ppa && apt-get update && apt-get install -y pandaseq

RUN perl -MCPAN -e 'install FAST'

WORKDIR /FACS

# Install pigz
RUN wget https://zlib.net/pigz/pigz-2.4.tar.gz && \
    tar -zxvf pigz-2.4.tar.gz && \
    cd pigz-2.4 && \
    make && \
    mv pigz ../ && \
    mv unpigz ../ && \
    cd ../ && \
    rm -rf pigz-2.4/ pigz-2.4.tar.gz

# Install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.10.1/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    rm -rf seqkit_linux_amd64.tar.gz

# Install ORFM
RUN wget https://github.com/wwood/OrfM/releases/download/v0.7.1/orfm-0.7.1_Linux_x86_64.tar.gz && \
    tar -zxf orfm-0.7.1_Linux_x86_64.tar.gz && \
    mv orfm-0.7.1_Linux_x86_64/orfm-0.7.1_Linux_x86_64 ./orfm && \
    rm -rf orfm-0.7.1_Linux_x86_64 orfm-0.7.1_Linux_x86_64.tar.gz

# Install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    rm -rf Trimmomatic-0.39.zip

