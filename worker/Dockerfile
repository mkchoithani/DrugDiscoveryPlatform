# GPU docking worker for Runpod Serverless.
#
# Local smina runs ~6.6 s per compound on CPU. Uni-Dock on a GPU is the same job
# at roughly a thousand times the throughput, which is what makes screening a
# million compounds cost about $8 rather than a fortnight of laptop time.
#
# Build and push:
#   docker build --platform linux/amd64 -t <registry>/dmta-docking:1 worker/
#   docker push <registry>/dmta-docking:1
# Then register it as a Runpod Serverless endpoint (see worker/README.md).

FROM nvidia/cuda:12.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates bzip2 && \
    rm -rf /var/lib/apt/lists/*

# Miniforge gives us conda-forge, which is where the chemistry packages live.
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/mf.sh && \
    bash /tmp/mf.sh -b -p $CONDA_DIR && rm /tmp/mf.sh

# rdkit for ligand handling; smina as the CPU fallback and reference scorer.
RUN conda install -y -c conda-forge python=3.11 rdkit numpy smina && conda clean -afy

# Uni-Dock: the GPU engine. Pinned by tag so a rebuild cannot silently change
# the scoring function underneath your results.
ARG UNIDOCK_VERSION=v1.1.3
RUN wget -q https://github.com/dptech-corp/Uni-Dock/releases/download/${UNIDOCK_VERSION}/unidock \
         -O /usr/local/bin/unidock && chmod +x /usr/local/bin/unidock || \
    echo "WARNING: Uni-Dock binary not fetched; worker will fall back to smina"

RUN pip install --no-cache-dir runpod

WORKDIR /app
COPY handler.py /app/handler.py

CMD ["python", "-u", "/app/handler.py"]
