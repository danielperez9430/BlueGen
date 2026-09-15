#!/usr/bin/env python3
"""Download the UCSC liftOver chain files BlueGen uses for genome-build
handling (RELEASE_PLAN 3.0.2):

  reference/liftover/hg38ToHg19.over.chain.gz   GRCh38 input → GRCh37 pipeline (1.2 MB)
  reference/liftover/hg19ToHg38.over.chain.gz   marker probe for headerless VCFs (0.2 MB)

`prs.py run` fetches them on first need; this script is for offline /
Docker setups. Source: https://hgdownload.soe.ucsc.edu/goldenPath/
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from bluegen.genome_build import CHAIN_URLS, ensure_chain  # noqa: E402


def main():
    ref_dir = Path(__file__).resolve().parents[2] / "reference"
    for name in CHAIN_URLS:
        path = ensure_chain(ref_dir, name)
        print(f"✓ {name}: {path} ({path.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
