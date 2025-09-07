# Dockerfile
FROM mambaorg/micromamba:1.5.8

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-lc"]

# Copy the pinned environment
COPY environment.yml /tmp/environment.yml

# Create the env and make it default
RUN micromamba create -y -f /tmp/environment.yml && \
    micromamba clean -a -y
ENV MAMBA_DEFAULT_ENV=exome
ENV CONDA_PREFIX=/opt/conda/envs/${MAMBA_DEFAULT_ENV}
ENV PATH=${CONDA_PREFIX}/bin:$PATH

# Keep workdir clean and writable
WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
