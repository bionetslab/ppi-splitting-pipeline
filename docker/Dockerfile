# Image for every pipeline process. .github/workflows/publish-docker.yml builds and pushes it
# to ghcr on each GitHub release; the container profiles in nextflow.config use its latest tag.
FROM mambaorg/micromamba:2.0.5

USER root

# procps: Nextflow's trace polling. ca-certificates: HTTPS for pip and the UniProt/Pfam/Gurobi calls.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# environment.yml is the only package source; do not add packages or a separate torch layer here.
ENV PIP_NO_CACHE_DIR=1
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# Nextflow puts only the root project's bin/ on PATH, so bake it in for embedding pipelines.
COPY bin/ /opt/ppi-splitting/bin/
RUN chmod -R a+rx /opt/ppi-splitting/bin
ENV PATH=/opt/ppi-splitting/bin:$PATH

USER $MAMBA_USER
