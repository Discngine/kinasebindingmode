FROM continuumio/miniconda3:latest
COPY . /opt/kinconform

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

RUN conda config --add channels conda-forge && conda install -y biopython pandas
WORKDIR /opt/kinconform

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD "python getKinaseStrAnnotations.py"
