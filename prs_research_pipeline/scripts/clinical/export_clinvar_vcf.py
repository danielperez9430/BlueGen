#!/usr/bin/env python3
"""
Export ClinVar pathogenic variants as annotated VCF.
Output: clinvar/clinvar_pathogenic.vcf.gz — ready for IGV or sharing.
"""

import sys
import os
import json
import logging
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)

HEADER = """##fileformat=VCFv4.2
##source=PRS_Research_Platform_v1.0.0
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance">
##INFO=<ID=CLNDN,Number=1,Type=String,Description="ClinVar disease name">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene symbol(s)">
##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description="ClinVar review status">
##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Confidence tier">
##INFO=<ID=DISEASE_DESC,Number=1,Type=String,Description="Disease description (MedGen)">
"""

def export(user_vcf, clinvar_json, output_vcf):
    logger.info("═══ ClinVar VCF Export ═══")

    # Load ClinVar data
    with open(clinvar_json) as f:
        data = json.load(f)
    variants = data.get("pathogenic_variants", [])
    logger.info(f"  Pathogenic variants: {len(variants)}")

    if not variants:
        logger.warning("  No variants to export")
        return 1

    # Create BED file of positions
    bed_path = Path(tempfile.mkdtemp()) / "regions.bed"
    with open(bed_path, "w") as fh:
        for v in variants:
            # Normalize chromosome name
            chrom = v["chrom"]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            fh.write(f"{chrom}\t{v['pos']-1}\t{v['pos']}\n")

    # Filter user VCF to these positions
    filtered_vcf = str(Path(output_vcf).with_suffix("")) + "_temp.vcf"
    cmd = ["bcftools", "view", "-R", str(bed_path), user_vcf, "-o", filtered_vcf]
    subprocess.run(cmd, check=True, capture_output=True)
    n_filtered = sum(1 for l in open(filtered_vcf) if not l.startswith("#"))
    logger.info(f"  Filtered VCF records: {n_filtered}")

    # Build variant lookup
    lookup = {}
    for v in variants:
        chrom = v["chrom"]
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        key = (chrom, str(v["pos"]), v.get("ref",".").upper(), v.get("alt",".").upper())
        lookup[key] = v

    # Annotate with ClinVar INFO
    output_lines = []
    with open(filtered_vcf) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#"):
                if line.startswith("##fileformat"):
                    output_lines.append(HEADER)
                elif not line.startswith("##"):
                    # Add our INFO headers before the #CHROM line
                    output_lines.append(line)
                else:
                    output_lines.append(line)
            else:
                parts = line.split("\t")
                if len(parts) < 8:
                    output_lines.append(line)
                    continue
                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3].upper(), parts[4].upper()
                key = (chrom, str(pos), ref, alt)
                v = lookup.get(key)
                if not v:
                    # Try without chr prefix
                    key2 = (chrom.replace("chr",""), str(pos), ref, alt)
                    v = lookup.get(key2)
                if v:
                    clnsig = v.get("clinical_significance",".")
                    clndn = (v.get("disease_name",".") or ".").replace(" ","_").replace("|",";")[:100]
                    gene = ",".join(v.get("genes",[])) or v.get("gene_info",".").split(":")[0] or "."
                    revstat = (v.get("review_status",".") or ".").replace(" ","_")
                    tier = v.get("confidence_tier",".")
                    desc = (v.get("disease_description","") or "").replace(" ","_")[:200]
                    info = f"CLNSIG={clnsig};CLNDN={clndn};GENEINFO={gene};CLNREVSTAT={revstat};CONFIDENCE={tier};DISEASE_DESC={desc}"
                    parts[7] = info + ("" if parts[7] == "." else ";" + parts[7])
                output_lines.append("\t".join(parts))

    # Output: filtered VCF + TSV table
    out_base = str(Path(output_vcf).with_suffix(""))

    # 1. Write filtered VCF (original format, no extra annotations)
    import gzip as gz
    vcf_out = out_base + "_filtered.vcf.gz"
    with gz.open(vcf_out, "wt") as fh:
        fh.write("\n".join(output_lines) + "\n")

    # 2. Write annotated TSV (human-readable, IGV-compatible)
    tsv_out = out_base + "_annotated.tsv"
    with open(tsv_out, "w") as fh:
        fh.write("CHROM\tPOS\tREF\tALT\tRSID\tGENE\tCLINICAL_SIGNIFICANCE\tCONFIDENCE_TIER\tDISEASE\tDESCRIPTION\tREVIEW_STATUS\n")
        for v in variants:
            fh.write("\t".join([
                v.get("chrom","."), str(v.get("pos",".")), v.get("ref","."), v.get("alt","."),
                v.get("rsid",".") or ".", ",".join(v.get("genes",[])) or ".",
                v.get("clinical_significance","."),
                v.get("confidence_tier","."),
                (v.get("disease_name",".") or ".").replace("\t"," ")[:120],
                (v.get("disease_description","") or "").replace("\t"," ")[:200],
                (v.get("review_status",".") or ".").replace("_"," "),
            ]) + "\n")

    # Cleanup
    os.remove(filtered_vcf)
    os.remove(str(bed_path))
    os.rmdir(str(bed_path.parent))

    size_kb = Path(vcf_out).stat().st_size / 1024
    tsv_kb = Path(tsv_out).stat().st_size / 1024
    logger.info(f"  ✅ {vcf_out} ({size_kb:.0f} KB) — filtered VCF")
    logger.info(f"  ✅ {tsv_out} ({tsv_kb:.0f} KB) — annotated table")
    logger.info(f"     VCF: open with IGV")
    logger.info(f"     TSV: open with Excel / Numbers")
    return 0

def main():
    import argparse
    p = argparse.ArgumentParser(description="Export ClinVar pathogenic variants as annotated VCF")
    p.add_argument("--vcf", required=True, help="User VCF")
    p.add_argument("--clinvar-json", default="clinvar/clinvar_pathogenic_variants.json")
    p.add_argument("--output", "-o", default="clinvar/clinvar_pathogenic.vcf.gz")
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    pipeline_root = Path(__file__).resolve().parent.parent.parent
    clinvar_json = args.clinvar_json
    if not Path(clinvar_json).is_absolute():
        clinvar_json = str(pipeline_root / clinvar_json)
    output_vcf = args.output
    if not Path(output_vcf).is_absolute():
        output_vcf = str(pipeline_root / output_vcf)

    return export(args.vcf, clinvar_json, output_vcf)

if __name__ == "__main__":
    sys.exit(main())
