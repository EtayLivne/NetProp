# syntax=docker/dockerfile:1

FROM python:3.9.9-slim

WORKDIR /NetProp
COPY . .

RUN pip3 install -r requirements.txt


CMD python3 networks/randomization.py 500 10 /hostfs/data/H_sapiens_aug_2020.net && \
    /hostfs/data/cov_to_human_ppi.cx /hostfs/data/symbol_to_entrezgene_2021.json && \
    /hostfs/randomized_networks