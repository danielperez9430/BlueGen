# BlueGen — reproducible runtime image (RELEASE_PLAN.md 2.1.3)
#
#   docker build --platform linux/amd64 -t bluegen .
#   docker run --rm \
#     -v "$PWD/prs_research_pipeline/reference:/app/prs_research_pipeline/reference" \
#     -v "$PWD/prs_research_pipeline/reports:/app/prs_research_pipeline/reports" \
#     -v "/path/to/sample.vcf.gz:/data/sample.vcf.gz:ro" \
#     bluegen run --full --vcf /data/sample.vcf.gz
#
# The image ships Python deps + PLINK 1.9 + PLINK 2.0 + bcftools/tabix + the
# WeasyPrint system libs. It does NOT ship the ~65 GB reference bundle: mount
# prs_research_pipeline/reference/ from the host (see README "Requirements").
#
# Build with --platform linux/amd64: PLINK 1.9 has no Linux ARM64 build, and
# the guard below refuses any other architecture. On Apple Silicon Docker
# Desktop runs the amd64 image under Rosetta/QEMU (slower, but correct).
FROM python:3.12-slim-bookworm

RUN [ "$(uname -m)" = "x86_64" ] || { \
      echo "BlueGen image must be built for linux/amd64 (PLINK 1.9 has no ARM64 Linux build):"; \
      echo "  docker build --platform linux/amd64 -t bluegen ."; exit 1; }

# Dated PLINK builds from cog-genomics (verified 2026-09-14). Override with
# --build-arg if you need a different build; keep README "System Tools" in sync.
ARG PLINK1_URL=https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260913.zip
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260914.zip

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        bcftools tabix \
        curl ca-certificates unzip bash \
        # WeasyPrint (PDF) runtime libs
        libcairo2 libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf-2.0-0 \
        libffi8 shared-mime-info fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    curl -fsSL "$PLINK1_URL" -o /tmp/plink1.zip; \
    unzip -o -q /tmp/plink1.zip -d /tmp/plink1; \
    install -m 0755 /tmp/plink1/plink /usr/local/bin/plink; \
    curl -fsSL "$PLINK2_URL" -o /tmp/plink2.zip; \
    unzip -o -q /tmp/plink2.zip -d /tmp/plink2; \
    install -m 0755 /tmp/plink2/plink2 /usr/local/bin/plink2; \
    rm -rf /tmp/plink1 /tmp/plink2 /tmp/plink1.zip /tmp/plink2.zip; \
    plink --version; plink2 --version

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
