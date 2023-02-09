# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/cellmaps_hierarchyeval
COPY ./ /tmp/cellmaps_hierarchyeval/
RUN pip install /tmp/cellmaps_hierarchyeval

RUN rm -rf /tmp/cellmaps_hierarchyeval

ENTRYPOINT ["/opt/conda/bin/cellmaps_hierarchyevalcmd.py"]

CMD ["--help"]
