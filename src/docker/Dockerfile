# syntax=docker/dockerfile:1

FROM python:3.9.9-slim

ARG conf_file="hostfs/conf.json"

WORKDIR /NetProp
COPY . .

RUN pip3 install -r requirements.txt

CMD py main.py ${conf_file}