#!/bin/bash

python3 networks/randomizations 1 1 \
        /NetProp/data/H_sapiens_aug_2020.net /NetProp/Data/cov_to_human_ppi.cx \
        /NetProp/data/symbol_to_entrezgene_2021.json \
        /NetProp/output
