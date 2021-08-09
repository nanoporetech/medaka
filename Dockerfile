FROM nvidia/cuda:10.1-cudnn7-runtime-ubuntu18.04

ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Oxford Nanopore Technologies"
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN \
    apt update \
    && apt install -yq --no-install-recommends \
        ca-certificates build-essential cmake curl wget git \
        zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev \
        libssl-dev libffi-dev \
        libreadline7 libreadline-dev sqlite3 libsqlite3-dev file \
        python3-all-dev python3-venv python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY . /tmp/medaka
RUN \
    cd /tmp/medaka \
    && make install_root \
    && cd / \
    && rm -rf /tmp/medaka \
    && medaka --help \
    && samtools --help


