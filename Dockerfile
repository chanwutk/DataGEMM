FROM r-base:4.4.1

WORKDIR /workspace

COPY ./setup.R /workspace/setup.R

RUN Rscript setup.R