# Image for every pipeline process. .github/workflows/publish-docker.yml builds and pushes it
# to ghcr on each GitHub release; the container profiles in nextflow.config use its latest tag.

# pixi.lock is the only package source; do not add packages here.
FROM ghcr.io/prefix-dev/pixi:0.67.1-bookworm-slim AS build
WORKDIR /app
COPY pixi.toml pixi.lock ./
RUN pixi install --locked -e default

# Only the env is copied over, at the same path (conda envs are not relocatable).
FROM debian:bookworm-slim
# procps: Nextflow's trace polling. ca-certificates: HTTPS for the UniProt/Pfam/Gurobi calls.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
ENV CONDA_PREFIX=/app/.pixi/envs/default
ENV PATH=/app/.pixi/envs/default/bin:$PATH

# Nextflow puts only the root project's bin/ on PATH, so bake it in for embedding pipelines.
COPY bin/ /opt/ppi-splitting/bin/
RUN chmod -R a+rx /opt/ppi-splitting/bin
ENV PATH=/opt/ppi-splitting/bin:$PATH
