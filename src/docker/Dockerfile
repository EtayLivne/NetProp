# syntax=docker/dockerfile:1

FROM python:3.9.9-slim

ARG workdir=/NetProp

WORKDIR ${workdir}
COPY . .

RUN pip3 install -r requirements.txt

ENV PYTHONPATH = "${PYTHONPATH}:${workdir}"

CMD ./runscript.sh