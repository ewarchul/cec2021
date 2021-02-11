FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
  apt-get -y dist-upgrade && \
  apt-get -y install --no-install-recommends --no-install-suggests \
  apt-utils curl g++ gcc gfortran git make libblas-dev libcurl4-gnutls-dev \
  pkg-config r-base libxml2-dev libssl-dev && \
  apt-get clean 

COPY ./dependencies.R /

RUN Rscript dependencies.R


