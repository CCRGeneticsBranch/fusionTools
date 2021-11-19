# syntax=docker/dockerfile:1
FROM python:3
COPY fusionTools /apps/
ENV PATH=/apps:${PATH}
ENV PERL5LIB=/apps/PfamScan:${PERL5LIB}
RUN pip install --upgrade gtfparse pyfaidx dataclasses pysam pyyaml Bio numpy pandas pybedtools
RUN apt-get update && apt-get install -y bedtools
RUN apt-get update && apt-get install -y hmmer
RUN apt-get update && apt-get install -y cpanminus
RUN apt-get update && apt-get install -y libipc-run-perl
RUN cpanm -v Moose
