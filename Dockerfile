# BlueGen — reproducible runtime image (RELEASE_PLAN.md 2.1.3 / 2.1.3b)
#
#   docker build -t bluegen .          # native on both amd64 and arm64
#   docker run --rm \
#     -v "$PWD/prs_research_pipeline/reference:/app/prs_research_pipeline/reference" \
#     -v "$PWD/prs_research_pipeline/reports:/app/prs_research_pipeline/reports" \
#     -v "/path/to/sample.vcf.gz:/data/sample.vcf.gz:ro" \
#     bluegen run --full --vcf /data/sample.vcf.gz
#
# Multi-arch. cog-genomics publishes no Linux ARM64 PLINK, so PLINK 1.9 and
# bcftools/tabix come from bioconda (same pinned versions on linux-64 and
# linux-aarch64 → identical tooling whatever the host). PLINK 2.0 has no
# ARM64 Linux build anywhere; it is only added on amd64 and the pipeline
# never invokes it (plink 1.9 does everything, plink2 is a `which` fallback).
#
# The image does NOT ship the ~65 GB reference bundle: mount
# prs_research_pipeline/reference/ from the host (see README "Requirements").
FROM python:3.12-slim-bookworm

# Set automatically by BuildKit: amd64 | arm64
ARG TARGETARCH

# Pinned tool versions (verified on bioconda for both arches, 2026-09-14).
ARG PLINK_VERSION=1.90b7.7
ARG BCFTOOLS_VERSION=1.24
ARG HTSLIB_VERSION=1.24
# PLINK 2.0 (amd64 only) — dated build from cog-genomics.
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260914.zip

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MAMBA_ROOT_PREFIX=/opt/micromamba

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates bzip2 unzip bash \
        # WeasyPrint (PDF) runtime libs
        libcairo2 libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf-2.0-0 \
        libffi8 shared-mime-info fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

# PLINK 1.9 + bcftools + tabix/bgzip via micromamba (bioconda), then expose
# them on PATH as plain binaries so nothing else needs to know about conda.
RUN set -eux; \
    case "$TARGETARCH" in \
        amd64) MM_ARCH=linux-64 ;; \
        arm64) MM_ARCH=linux-aarch64 ;; \
        *) echo "Unsupported TARGETARCH=$TARGETARCH (need amd64 or arm64)"; exit 1 ;; \
    esac; \
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${MM_ARCH}/latest" \
        | tar -xj -C /usr/local bin/micromamba; \
    micromamba create -y -p /opt/bio \
        -c conda-forge -c bioconda \
        "plink=${PLINK_VERSION}" "bcftools=${BCFTOOLS_VERSION}" "htslib=${HTSLIB_VERSION}"; \
    micromamba clean -a -y; \
    for b in plink bcftools tabix bgzip; do ln -sf "/opt/bio/bin/$b" "/usr/local/bin/$b"; done; \
    plink --version; bcftools --version | head -1; tabix --version | head -1

# PLINK 2.0: optional, amd64 only (no ARM64 Linux build exists).
RUN set -eux; \
    if [ "$TARGETARCH" = "amd64" ]; then \
        curl -fsSL "$PLINK2_URL" -o /tmp/plink2.zip; \
        unzip -o -q /tmp/plink2.zip -d /tmp/plink2; \
        install -m 0755 /tmp/plink2/plink2 /usr/local/bin/plink2; \
        rm -rf /tmp/plink2 /tmp/plink2.zip; \
        plink2 --version; \
    else \
        echo "plink2 skipped on $TARGETARCH (no Linux ARM64 build; the pipeline only needs plink 1.9)"; \
    fi

WORKDIR /app

# Dependencies first so code edits don't invalidate this layer.
COPY prs_research_pipeline/requirements.lock prs_research_pipeline/requirements.lock
RUN pip install --no-cache-dir -r prs_research_pipeline/requirements.lock

# Code + curated data (reference data and outputs are excluded by .dockerignore).
COPY . /app
# Editable install: prs.py resolves prs_research_pipeline/ relative to itself.
RUN pip install --no-cache-dir --no-deps -e . && bluegen --help >/dev/null

# Output/reference mount points (bind-mount these from the host).
VOLUME ["/app/prs_research_pipeline/reference", "/app/prs_research_pipeline/reports", "/data"]

ENTRYPOINT ["bluegen"]
CMD ["status"]
