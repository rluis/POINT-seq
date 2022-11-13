FROM continuumio/miniconda3

COPY condaPOINTseq.yml .
RUN conda env update -n root -f condaPOINTseq.yml && conda clean -a

