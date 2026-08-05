FROM mambaorg/micromamba:latest

LABEL maintainer="matt@datatecnica.com"
LABEL version="1.0.0"
LABEL image="snRNAseq"
LABEL description="Container image for snRNAseq pipeline. Analysis specific container can be built on top of this image to any additional dependencies required for specific analysis not included in this repository."

USER root

# Install system dependencies required for R packages and compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    libgl1 \
    libglib2.0-0 \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libudunits2-dev \
    libgdal-dev \
    libmagick++-dev \
    librsvg2-dev \
    libfftw3-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libdeflate-dev \
    libuv1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy environment definition and create the conda environment
COPY env/environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml \
    && micromamba clean -afy

# Activate the snRNAseq environment for all subsequent RUN steps
SHELL ["micromamba", "run", "-n", "snRNAseq", "/bin/bash", "-c"]

# Install cell_type_mapper from GitHub
RUN git clone https://github.com/AllenInstitute/cell_type_mapper.git \
    && cd cell_type_mapper \
    && pip install .

# Copy and run the R library installation script
COPY env/install_libraries.r /tmp/install_libraries.r
RUN Rscript /tmp/install_libraries.r

# Copy the full project into the container
WORKDIR /snRNAseq
COPY . .

# Ensure the snRNAseq environment is on PATH at runtime
ENV PATH="/opt/conda/envs/snRNAseq/bin:$PATH"
ENV CONDA_DEFAULT_ENV=snRNAseq

# Default: run the full workflow; override CMD at runtime for individual modules
ENTRYPOINT ["micromamba", "run", "-n", "snRNAseq", "snakemake"]
CMD ["--snakefile", "Snakefile", "--cores", "all"]
