#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PGS RSID MAPPER — scripts/benchmarking/pgs_rsid_mapper.py                  ║
║                                                                            ║
║   Converts rsID-based variant identifiers to chr:pos format for PLINK      ║
║   scoring compatibility with DeepVariant VCFs (which lack rsIDs).          ║
║                                                                            ║
║   Mapping sources (in priority order):                                      ║
║     1. SNP database (snp_database_annotated.csv) — rsid + chrom + pos      ║
║     2. Local JSON cache                                                     ║
║                                                                            ║
║   Usage:                                                                    ║
║       mapper = RsidMapper()                                                 ║
║       mapper.load_from_snp_db("data/snp_database_annotated.csv")            ║
║       mapper.load_from_cache("prs/pgs_scores/rsid_mapping_cache.json")     ║
║       chr_pos = mapper.lookup("rs7903146")  # → "10:114758349"              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional, List

import pandas as pd

logger = logging.getLogger(__name__)


class RsidMapper:
    """Maps rsIDs to chr:pos format for PGS Catalog scoring compatibility."""

    def __init__(self):
        self._mapping: Dict[str, str] = {}  # rsid → "chr:pos"
        self._stats = {"snp_db": 0, "cache": 0, "miss": 0}

    # ── Data Loading ───────────────────────────────────────────────────────

    def load_from_snp_db(self, csv_path: str) -> int:
        """Index the SNP database for rsID → chr:pos mappings.

        Returns number of entries loaded.
        """
        path = Path(csv_path)
        if not path.exists():
            logger.warning(f"  SNP database not found: {csv_path}")
            return 0

        try:
            db = pd.read_csv(csv_path, dtype=str)
            count = 0
            for _, row in db.iterrows():
                rsid = str(row.get("rsid", "")).strip()
                chrom = str(row.get("chrom", "")).strip()
                pos = str(row.get("pos", "")).strip()
                if rsid and rsid.startswith("rs") and chrom and pos:
                    # Clean chrom prefix
                    if chrom.startswith("chr"):
                        chrom = chrom[3:]
                    self._mapping[rsid] = f"{chrom}:{pos}"
                    count += 1
            logger.info(f"  Loaded {count} rsID mappings from SNP database")
            self._stats["snp_db"] = count
            return count
        except Exception as e:
            logger.warning(f"  SNP database indexing error: {e}")
            return 0

    def load_from_cache(self, json_path: str) -> int:
        """Load cached rsID → chr:pos mappings from JSON.

        Returns number of entries loaded.
        """
        path = Path(json_path)
        if not path.exists():
            return 0

        try:
            data = json.loads(path.read_text())
            mapping = data.get("mapping", {})
            loaded = 0
            for rsid, entry in mapping.items():
                if isinstance(entry, dict):
                    chr_pos = f"{entry.get('chrom', '')}:{entry.get('pos', '')}"
                else:
                    chr_pos = str(entry)
                if ":" in chr_pos:
                    self._mapping[rsid] = chr_pos
                    loaded += 1
            logger.info(f"  Loaded {loaded} entries from cache: {json_path}")
            self._stats["cache"] = loaded
            return loaded
        except Exception as e:
            logger.warning(f"  Cache load error: {e}")
            return 0

    def save_cache(self, json_path: str) -> int:
        """Save current mappings to JSON cache.

        Returns number of entries saved.
        """
        path = Path(json_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        cache_data = {
            "generated_date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M UTC"),
            "n_entries": len(self._mapping),
            "mapping": {
                rsid: {"chrom": chr_pos.split(":")[0], "pos": int(chr_pos.split(":")[1]),
                       "source": "snp_db" if rsid in self._mapping else "unknown"}
                for rsid, chr_pos in self._mapping.items()
            },
        }

        path.write_text(json.dumps(cache_data, indent=2))
        logger.info(f"  Saved {len(self._mapping)} entries to cache: {json_path}")
        return len(self._mapping)

    # ── Lookup API ─────────────────────────────────────────────────────────

    def lookup(self, rsid: str) -> Optional[str]:
        """Look up a single rsID → chr:pos.

        Returns "chr:pos" string or None if not found.
        """
        # Normalize rsID
        rsid = rsid.strip()
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"

        result = self._mapping.get(rsid)
        if result is None:
            self._stats["miss"] += 1
        return result

    def batch_lookup(self, rsids: List[str]) -> Dict[str, Optional[str]]:
        """Convert multiple rsIDs to chr:pos.

        Returns {rsid: "chr:pos" or None} dict.
        """
        return {rsid: self.lookup(rsid) for rsid in rsids}

    # ── Stats ──────────────────────────────────────────────────────────────

    @property
    def n_entries(self) -> int:
        return len(self._mapping)

    @property
    def stats(self) -> Dict[str, int]:
        return dict(self._stats)
