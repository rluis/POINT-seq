FROM ubuntu:20.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ENV PATH="/root/miniconda3/lib:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/lib:${PATH}"

RUN apt-get update
RUN apt-get install -y rename curl libcurl4 wget && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --set rename /usr/bin/file-rename

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

COPY condaPOINTseq.yml .
RUN conda env update -n root -f condaPOINTseq.yml && conda clean -a

