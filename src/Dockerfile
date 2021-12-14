# syntax=docker/dockerfile:1

FROM python:3.9.9-slim

#ARG workdir=/NetProp

WORKDIR /NetProp
COPY requirements.txt .

RUN pip3 install -r requirements.txt

# ENV PATH="${PATH}:${workdir}"

CMD python3 /NetProp/hostfs/networks/randomization.py 500 10 /hostfs/data/H_sapiens_aug_2020.net && \
    /hostfs/data/cov_to_human_ppi.cx /hostfs/data/symbol_to_entrezgene_2021.json && \
    /hostfs/randomized_networks