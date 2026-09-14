#!/usr/bin/env python3
"""
🧬 BlueGen — Interactive Dashboard
   Run: streamlit run dashboard.py
   Reads all pipeline JSON outputs. No recomputation needed.
"""

import streamlit as st
import json
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from collections import Counter

st.set_page_config(
    page_title="BlueGen",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

REPO_ROOT = Path(__file__).parent
# BLUEGEN_PIPELINE_DIR lets tests point the dashboard at an empty directory
# to prove every page degrades to "no data" instead of a traceback.
PIPELINE = Path(os.environ.get("BLUEGEN_PIPELINE_DIR", REPO_ROOT / "prs_research_pipeline"))

# Single source of truth for the version (IMPROVEMENT_PLAN.md TIER 0.1) —
# never hardcode a "vX.Y.Z" literal here; tests/test_version_consistency.py
# scans this file for exactly that.
import sys
sys.path.insert(0, str(REPO_ROOT / "prs_research_pipeline" / "scripts"))
from utils.constants import PIPELINE_VERSION

# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

# Cache on the *absolute* path: a relative key would survive a change of
# PIPELINE (e.g. the BLUEGEN_PIPELINE_DIR override in tests) and hand back
# another directory's data.
@st.cache_data
def _read_json(abs_path: str):
    p = Path(abs_path)
    if p.exists():
        try:
            with open(p) as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}
    return {}

@st.cache_data
def _read_csv(abs_path: str):
    p = Path(abs_path)
    if p.exists():
        try:
            return pd.read_csv(p)
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()

def load_json(path):
    return _read_json(str(PIPELINE / path))

def load_csv(path):
    return _read_csv(str(PIPELINE / path))

def safe_get(d, *keys, default="N/A"):
    for k in keys:
        if isinstance(d, dict):
            d = d.get(k, default)
        else:
            return default
    return d if d is not None else default

RISK_COLORS = {"HIGH": "#e74c3c", "ELEVATED": "#f39c12", "AVERAGE": "#95a5a6",
               "LOW": "#27ae60", "high": "#e74c3c", "elevated": "#f39c12",
               "medium": "#95a5a6", "average": "#95a5a6", "low": "#27ae60"}

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════════════════════

st.sidebar.markdown("# 🧬 BlueGen")
st.sidebar.markdown("### Navigation")

page = st.sidebar.radio(
    "", [
        "📊 Overview",
        "🧬 PRS Results",
        "🥗 Recommendations",
        "🩺 PGS Catalog",
        "🔬 ClinVar Pathogenic",
        "💊 Pharmacogenomics",
        "🌍 Ancestry",
        "🦴 Archaic DNA",
        "📋 Raw Data",
    ],
    label_visibility="collapsed",
)

st.sidebar.markdown("---")
st.sidebar.caption("Reads from `prs_research_pipeline/` outputs.")
st.sidebar.caption("No recomputation.")

# Load all data once
prs = load_json("prs/PRS_RESULT.json")
anc = load_json("science/ANCESTRY_MODEL.json")
integrity = load_json("FINAL_SCIENTIFIC_SCORE.json")
clinvar = load_json("clinvar/clinvar_pathogenic_variants.json")
pharmgkb = load_json("pharmgkb/pharmgkb_drug_report.json")
deep_anc = load_json("ancestry/deep_ancestry.json")
validation = load_json("science/global_validation_report.json")
pgs_cal = load_csv("prs/pgs_scores/pgs_calibrated.csv")
recommendations = load_json("data/trait_recommendations.json")

# ═══════════════════════════════════════════════════════════════════════════════
# OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════════

if page == "📊 Overview":
    st.title("📊 BlueGen")
    st.caption("Polygenic Risk Score + ClinVar + Pharmacogenomics + Ancestry")

    # First row: 3 metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        traits = len(safe_get(prs, "prs_entries", default=[]))
        st.metric("🧬 PRS Traits", traits if traits != "N/A" else 10)
    with col2:
        pop = safe_get(anc, "assigned_population", default="EUR")
        st.metric("🌍 Ancestry", pop)
    with col3:
        cv_variants = len(safe_get(clinvar, "pathogenic_variants", default=[]))
        st.metric("🔬 Pathogenic Variants", cv_variants)

    # Second row: 3 metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        high_conf = safe_get(clinvar, "pathogenic_variant_summary", "high_confidence_count", default=0)
        st.metric("🏅 High Confidence", high_conf)
    with col2:
        pgx = len(safe_get(pharmgkb, "pharmacogenomic_findings", default=[]))
        st.metric("💊 Drug Findings", pgx)
    with col3:
        score = safe_get(integrity, "scientific_integrity_score", default=0)
        st.metric("📊 Integrity Score", f"{score}/100")

    # Third row: 3 metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("🧬 mtDNA Haplogroup", safe_get(deep_anc, "mt_dna", "haplogroup", default="—"))
    with col2:
        st.metric("📦 Pipeline", f"v{PIPELINE_VERSION}")
    with col3:
        st.metric("🗃️ ClinVar DB", f"{safe_get(clinvar, 'metadata', 'user_vcf_total_variants', default=0):,} variants")

    st.markdown("---")

    # PRS quick view
    st.subheader("🧬 PRS Risk Profile")
    entries = safe_get(prs, "prs_entries", default=[])
    if entries:
        df_prs = pd.DataFrame(entries)
        if "trait" in df_prs.columns and "population_zscore" in df_prs.columns:
            df_prs = df_prs.sort_values("population_zscore")
            colors = ["#27ae60" if z < 0 else "#f39c12" if z < 1 else "#e74c3c" for z in df_prs["population_zscore"]]
            fig = px.bar(df_prs, x="population_zscore", y="trait", color="population_zscore",
                         color_continuous_scale=["#27ae60", "#f39c12", "#e74c3c"],
                         title="PRS Z-Scores by Trait (population-calibrated)")
            fig.add_vline(x=1.0, line_dash="dash", line_color="#f39c12", annotation_text="Elevated")
            fig.add_vline(x=-1.0, line_dash="dash", line_color="#27ae60", annotation_text="Reduced")
            st.plotly_chart(fig, width="stretch")
    else:
        st.info("No PRS data. Run pipeline first.")

    # ClinVar summary
    st.subheader("🔬 ClinVar Pathogenic — Confidence Tiers")
    tiers = safe_get(clinvar, "pathogenic_variant_summary", "by_confidence_tier", default={})
    if tiers:
        tier_df = pd.DataFrame([
            {"Tier": "🏅 High (Expert Panel)", "Count": tiers.get("high", 0), "Color": "#27ae60"},
            {"Tier": "✓ Moderate (Multi-Lab)", "Count": tiers.get("moderate", 0), "Color": "#2e86c1"},
            {"Tier": "⚠️ Low (Single Lab)", "Count": tiers.get("low", 0), "Color": "#f39c12"},
            {"Tier": "❓ Very Low (No Criteria)", "Count": tiers.get("very_low", 0), "Color": "#95a5a6"},
        ])
        fig = px.bar(tier_df, x="Tier", y="Count", color="Tier",
                     color_discrete_map={t["Tier"]: t["Color"] for _, t in tier_df.iterrows()},
                     title="Variant Confidence Distribution")
        st.plotly_chart(fig, width="stretch")

    # Integrity score gauge
    if score:
        st.subheader("🏅 Scientific Integrity")
        fig = go.Figure(go.Indicator(
            mode="gauge+number", value=score,
            title={"text": "Scientific Integrity Score"},
            gauge={"axis": {"range": [0, 100]},
                   "bar": {"color": "#27ae60" if score >= 75 else "#f39c12" if score >= 50 else "#e74c3c"},
                   "steps": [{"range": [0, 50], "color": "#fadbd8"},
                             {"range": [50, 75], "color": "#fdebd0"},
                             {"range": [75, 100], "color": "#d5f5e3"}]}))
        st.plotly_chart(fig, width="stretch")

# ═══════════════════════════════════════════════════════════════════════════════
# PRS RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🧬 PRS Results":
    st.title("🧬 Polygenic Risk Scores")
    st.caption("Population-calibrated PRS with z-scores and percentiles")

    entries = safe_get(prs, "prs_entries", default=[])
    if not entries:
        st.warning("No PRS data. Run `python prs.py run --full --vcf sample.vcf.gz`")
    else:
        df = pd.DataFrame(entries)
        if "population_zscore" in df.columns and "trait" in df.columns:
            df["risk_color"] = df["population_zscore"].apply(
                lambda z: "#e74c3c" if z >= 2 else ("#f39c12" if z >= 1 else "#27ae60"))
            df = df.sort_values("population_zscore", ascending=False)

            cols = st.columns(3)
            with cols[0]:
                st.metric("Traits analyzed", len(df))
            with cols[1]:
                high = len(df[df["population_zscore"] >= 1])
                st.metric("Elevated risk (>1σ)", high)
            with cols[2]:
                low = len(df[df["population_zscore"] <= -1])
                st.metric("Reduced risk (<-1σ)", low)

            # Detailed table
            st.subheader("All Traits")
            display_df = df[["trait", "population_zscore", "population_percentile", "risk_category"]].copy()
            display_df.columns = ["Trait", "Z-Score", "Percentile", "Risk Category"]
            display_df["Z-Score"] = display_df["Z-Score"].round(2)
            display_df["Percentile"] = display_df["Percentile"].round(1)
            st.dataframe(display_df, width="stretch", hide_index=True)

            # Bar chart
            fig = px.bar(df, x="population_zscore", y="trait",
                         title="PRS Z-Scores",
                         color="population_zscore",
                         color_continuous_scale=["#27ae60", "#f5f5f5", "#e74c3c"],
                         range_color=[-2, 2])
            fig.add_vline(x=1.0, line_dash="dash", line_color="#f39c12")
            fig.add_vline(x=-1.0, line_dash="dash", line_color="#27ae60")
            st.plotly_chart(fig, width="stretch")

            # Radar chart
            st.subheader("Risk Radar")
            radar_df = df[["trait", "population_zscore"]].copy()
            radar_df["abs_z"] = radar_df["population_zscore"].abs()
            fig = go.Figure()
            fig.add_trace(go.Scatterpolar(
                r=radar_df["population_zscore"].tolist(),
                theta=radar_df["trait"].tolist(),
                fill="toself",
                name="Z-Score",
            ))
            fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[-2, 2])))
            st.plotly_chart(fig, width="stretch")

# ═══════════════════════════════════════════════════════════════════════════════
# RECOMMENDATIONS (curated, evidence-cited — data/trait_recommendations.json)
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🥗 Recommendations":
    st.title("🥗 Actionable Recommendations")
    st.caption("Curated, evidence-cited guidance per trait (PubMed PMID / NIH fact sheets), "
               "joined with your population-calibrated PRS where the trait was scored.")

    recs = {k: v for k, v in recommendations.items()
            if not k.startswith("_") and isinstance(v, dict)}
    entries = safe_get(prs, "prs_entries", default=[])
    by_trait = {str(e.get("trait", "")).lower(): e for e in entries if isinstance(e, dict)}

    if not recs:
        st.info("No recommendations file found (`data/trait_recommendations.json`).")
    else:
        lang = st.radio("Language / Idioma", ["en", "es"], horizontal=True,
                        format_func=lambda x: {"en": "🇬🇧 English", "es": "🇪🇸 Español"}[x])
        col1, col2, col3 = st.columns(3)
        with col1:
            levels = sorted({str(v.get("evidence_level", "?")) for v in recs.values()})
            level_filter = st.multiselect("Evidence level", levels, default=levels)
        with col2:
            only_scored = st.checkbox("Only traits scored in my PRS", value=bool(by_trait))
        with col3:
            search = st.text_input("Search trait / gene", "")

        rows = []
        for trait, rec in recs.items():
            e = by_trait.get(trait.lower())
            if only_scored and e is None:
                continue
            if str(rec.get("evidence_level", "?")) not in level_filter:
                continue
            blob = f"{trait} {rec.get('recommendation_en', '')} {rec.get('recommendation_es', '')}".lower()
            if search and search.lower() not in blob:
                continue
            z = e.get("population_zscore") if e else None
            rows.append((trait, rec, e, z))

        # Most extreme z-scores first, unscored traits last (alphabetical)
        rows.sort(key=lambda r: (-abs(r[3]) if isinstance(r[3], (int, float)) else 1e9, r[0]))

        m1, m2, m3 = st.columns(3)
        m1.metric("Recommendations", len(rows))
        m2.metric("Curated in panel", len(recs))
        m3.metric("Matched to your PRS", sum(1 for r in rows if r[2] is not None))
        if not rows:
            st.info("Nothing matches the current filters.")

        for trait, rec, e, z in rows:
            if e is not None and isinstance(z, (int, float)):
                badge = "🔴" if z >= 1 else "🟢" if z <= -1 else "⚪"
                head = f"{badge} **{trait}** — z = {z:+.2f}, percentile {e.get('population_percentile', 0):.0f}, {e.get('risk_category', '—')}"
            else:
                head = f"◌ **{trait}** — not scored in this run"
            with st.expander(head):
                st.markdown(rec.get(f"recommendation_{lang}") or rec.get("recommendation_en", "—"))
                c1, c2 = st.columns([1, 3])
                c1.metric("Evidence level", rec.get("evidence_level", "?"))
                c2.caption(f"**Source:** {rec.get('reference', '—')}")
                if e is not None:
                    c2.caption(f"SNPs used: {e.get('n_snps_used', '?')}/{e.get('n_snps_total', '?')} · "
                               f"95% CI [{e.get('ci_95_lower', '?')}, {e.get('ci_95_upper', '?')}] · "
                               f"calibrated vs {e.get('assigned_population', '?')}")

        st.markdown("---")
        st.caption("⚠️ Research use only. Discuss any dietary or lifestyle change with your clinician or dietitian.")

# ═══════════════════════════════════════════════════════════════════════════════
# PGS CATALOG (external validation — prs/pgs_scores/pgs_calibrated.csv)
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🩺 PGS Catalog":
    st.title("🩺 PGS Catalog — External Validation")
    st.caption("Published polygenic scores (EBI PGS Catalog) applied to your genotype and "
               "calibrated against the 1000 Genomes reference.")

    if pgs_cal.empty or "z_score" not in pgs_cal.columns:
        st.info("No PGS Catalog results. Run `python prs.py run --full --vcf sample.vcf.gz`.")
    else:
        df = pgs_cal.copy()
        df["z_score"] = pd.to_numeric(df["z_score"], errors="coerce")
        df["percentile"] = pd.to_numeric(df.get("percentile"), errors="coerce")
        df = df.dropna(subset=["z_score"]).sort_values("z_score", ascending=False)

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Scores calibrated", len(df))
        c2.metric("High (z ≥ 2)", int((df["z_score"] >= 2).sum()))
        c3.metric("Elevated (1 ≤ z < 2)", int(((df["z_score"] >= 1) & (df["z_score"] < 2)).sum()))
        if "reliable" in df.columns:
            c4.metric("Flagged reliable", int(df["reliable"].astype(str).str.lower().eq("true").sum()))

        ref_cols = [c for c in df.columns if c.endswith("_mean")]
        ref_pop = ref_cols[0].replace("_mean", "").upper() if ref_cols else "1000G"
        st.caption(f"Reference distribution: 1000 Genomes **{ref_pop}** for every score in this file.")

        show = [c for c in ["pgs_id", "trait", "n_snps", "z_score", "percentile", "risk_category", "reliable"]
                if c in df.columns]
        table = df[show].copy()
        table["z_score"] = table["z_score"].round(2)
        if "percentile" in table.columns:
            table["percentile"] = table["percentile"].round(1)
        st.dataframe(table, width="stretch", hide_index=True)

        top = df.reindex(df["z_score"].abs().sort_values(ascending=False).index).head(30)
        top = top.sort_values("z_score")
        fig = px.bar(top, x="z_score", y="trait", orientation="h",
                     color="z_score", color_continuous_scale=["#27ae60", "#f5f5f5", "#e74c3c"],
                     range_color=[-3, 3], title="PGS z-scores (30 most extreme)",
                     hover_data=[c for c in ["pgs_id", "n_snps", "percentile"] if c in top.columns])
        fig.add_vline(x=1.0, line_dash="dash", line_color="#f39c12")
        fig.add_vline(x=-1.0, line_dash="dash", line_color="#27ae60")
        fig.update_layout(height=max(400, 22 * len(top)))
        st.plotly_chart(fig, width="stretch")

        st.markdown("---")
        st.markdown("""
        **📖 Reading this page**
        - A PGS is a published score built from a large GWAS; the number of SNPs is what the authors used, not what BlueGen curated.
        - z = how many standard deviations your score sits from the reference population mean; percentile is the share of that population below you.
        - A very large |z| (≫ 3) on a score with few matched SNPs usually means low coverage of that score in your VCF rather than extreme risk — check `n_snps` and `reliable`.
        """)

# ═══════════════════════════════════════════════════════════════════════════════
# CLINVAR
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🔬 ClinVar Pathogenic":
    st.title("🔬 ClinVar Pathogenic Variants")
    st.caption("Genome-wide pathogenic/likely pathogenic variant annotation")

    variants = safe_get(clinvar, "pathogenic_variants", default=[])
    summary = safe_get(clinvar, "pathogenic_variant_summary", default={})
    meta = safe_get(clinvar, "metadata", default={})

    if not variants:
        st.warning("No ClinVar data. Run `python prs.py run --clinvar --vcf sample.vcf.gz`")
    else:
        # Top metrics
        cols = st.columns(4)
        with cols[0]:
            st.metric("Total Pathogenic/Likely", len(variants))
        with cols[1]:
            st.metric("High Confidence", summary.get("high_confidence_count", 0))
        with cols[2]:
            st.metric("With Descriptions", sum(1 for v in variants if v.get("disease_description")))
        with cols[3]:
            st.metric("ClinVar Matches", f"{meta.get('exact_matches', 0):,}")

        # Filters
        st.subheader("Filters")
        col1, col2, col3 = st.columns(3)
        with col1:
            tier_filter = st.multiselect(
                "Confidence Tier",
                ["high", "moderate", "low", "very_low"],
                default=["high", "moderate", "low", "very_low"],
                format_func=lambda x: {"high": "🏅 High", "moderate": "✓ Moderate", "low": "⚠️ Low", "very_low": "❓ Very Low"}[x],
            )
        with col2:
            sig_filter = st.multiselect(
                "Clinical Significance",
                ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Risk_allele"],
                default=["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Risk_allele"],
            )
        with col3:
            search = st.text_input("Search gene or disease", "")

        # Filter
        filtered = [v for v in variants
                    if v.get("confidence_tier") in tier_filter
                    and v.get("clinical_significance") in sig_filter
                    and (search.lower() in str(v).lower() if search else True)]

        st.caption(f"Showing {len(filtered)} of {len(variants)} variants")

        # Table
        rows = []
        for v in filtered[:200]:
            rows.append({
                "rsID": v.get("rsid") or "—",
                "Gene": ", ".join(v.get("genes", []))[:30],
                "Position": f"{v.get('chrom','?')}:{v.get('pos','?')}",
                "Significance": v.get("clinical_significance", "—"),
                "Confidence": v.get("confidence_tier", "—"),
                "Disease": (v.get("disease_name", "—") or "—").replace("_", " ")[:80],
                "Description": (v.get("disease_description", "") or "")[:120],
            })
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True,
                     column_config={"Description": st.column_config.TextColumn(width="large")})

        # Gene chart
        gene_counts = Counter()
        for v in filtered:
            for g in v.get("genes", []):
                gene_counts[g] += 1
        if gene_counts:
            top_genes = dict(gene_counts.most_common(15))
            fig = px.bar(x=list(top_genes.values()), y=list(top_genes.keys()),
                         orientation="h", title="Top Genes with Pathogenic Variants",
                         labels={"x": "Variants", "y": "Gene"})
            st.plotly_chart(fig, width="stretch")

# ═══════════════════════════════════════════════════════════════════════════════
# PHARMACOGENOMICS
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "💊 Pharmacogenomics":
    st.title("💊 Pharmacogenomics — Drug Response")
    st.caption("CPIC/PharmGKB-guided drug-gene associations from your genotype")

    findings = safe_get(pharmgkb, "pharmacogenomic_findings", default=[])
    summary = safe_get(pharmgkb, "summary", default={})

    if not findings:
        st.info("No pharmacogenomic data. Run `python prs.py run --clinvar --vcf sample.vcf.gz`")
    else:
        cols = st.columns(4)
        with cols[0]:
            st.metric("Total Findings", len(findings))
        with cols[1]:
            st.metric("Genes Affected", len(set(f["gene"] for f in findings)))
        with cols[2]:
            st.metric("Drugs Affected", len(set(f["drug"] for f in findings)))
        with cols[3]:
            crit = summary.get("by_actionability", {})
            st.metric("Critical/Important", crit.get("critical", 0) + crit.get("important", 0))

        st.markdown("---")

        # Findings by gene
        for finding in findings:
            icon = {"critical": "🔴", "important": "🟠", "informative": "🟡"}.get(finding.get("actionability", ""), "⚪")
            with st.expander(f"{icon} **{finding['gene']}** → **{finding['drug']}** ({finding.get('drug_class', '')}) — {finding.get('phenotype', '')}"):
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.markdown(f"**Recommendation:** {finding.get('recommendation_en', 'N/A')}")
                    guide_summary = finding.get("guideline_summary", "")
                    if isinstance(guide_summary, dict):
                        guide_summary = guide_summary.get("html", "")
                    if guide_summary:
                        st.markdown("---")
                        # Strip HTML tags for clean display
                        import re as _re
                        clean = _re.sub(r"<[^>]*>", "", str(guide_summary))
                        st.markdown(f"*{clean[:500]}*")
                with col2:
                    st.metric("Copies", finding.get("copies", 0))
                    st.caption(f"Variant: {finding.get('rsid', '?')} {finding.get('star_allele', '')}")
                    st.caption(f"CPIC Level: {finding.get('cpic_level', '?')}")
                    if finding.get("guideline_name"):
                        st.caption(f"Guideline: {finding.get('guideline_source', '')}")

        # Summary chart
        st.subheader("Drugs Affected")
        drug_genes = {}
        for f in findings:
            drug_genes[f["drug"]] = f["gene"]
        df_drugs = pd.DataFrame([
            {"Drug": drug, "Gene": gene, "Actionability": next(
                (f["actionability"] for f in findings if f["drug"] == drug), "informative")}
            for drug, gene in drug_genes.items()
        ])
        st.dataframe(df_drugs, width="stretch", hide_index=True)

        st.markdown("---")
        st.markdown("""
        **📖 About CPIC Guidelines:**
        - **CPIC** = Clinical Pharmacogenetics Implementation Consortium (NIH-funded, U.S.)
        - **DPWG** = Dutch Pharmacogenetics Working Group (European)
        - **Level A/B** = strong evidence → prescribing change recommended
        - **Level C/D** = moderate evidence → consider, not mandatory
        - ⚠️ Do NOT change medications without consulting your doctor.
        """)

# ═══════════════════════════════════════════════════════════════════════════════
# ANCESTRY
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🌍 Ancestry":
    st.title("🌍 Ancestry & Deep Ancestry")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Continental Ancestry")
        pop = safe_get(anc, "assigned_population", default="EUR")
        conf = safe_get(anc, "confidence", default="—")
        st.metric("Assigned Population", f"{pop} ({conf})")

        probs = safe_get(anc, "posterior_probabilities", default={})
        if probs:
            df_prob = pd.DataFrame([
                {"Population": p, "Probability": prob * 100}
                for p, prob in sorted(probs.items(), key=lambda x: -x[1])
            ])
            fig = px.bar(df_prob, x="Population", y="Probability", title="Ancestry Probabilities",
                         color="Population", color_discrete_sequence=px.colors.qualitative.Set2)
            st.plotly_chart(fig, width="stretch")

    with col2:
        st.subheader("Deep Ancestry")
        mtdna = safe_get(deep_anc, "mt_dna", default={})
        ydna = safe_get(deep_anc, "y_dna", default={})
        st.metric("mtDNA Haplogroup", mtdna.get("haplogroup", "—"))
        st.caption(mtdna.get("description", ""))
        st.metric("Y-DNA Haplogroup", ydna.get("haplogroup", "—"))
        st.caption(ydna.get("description", ""))

    # Sub-continental
    subcont = safe_get(deep_anc, "sub_continental", default={})
    if subcont:
        st.subheader("Sub-Continental Reference Populations")
        subs = subcont.get("sub_populations_available", [])
        if subs:
            df_subs = pd.DataFrame(subs)
            st.dataframe(df_subs[["code", "name", "description"]], width="stretch", hide_index=True)

    # PCA plot
    st.subheader("PCA Plot — 1000 Genomes Reference")
    pca_ref = PIPELINE / "pca" / "1000G_pcs.eigenvec"
    pca_target = PIPELINE / "pca" / "target_pcs.eigenvec"
    if pca_ref.exists() and pca_target.exists():
        try:
            ref_df = pd.read_csv(pca_ref, sep=r"\s+")
            target_df = pd.read_csv(pca_target, sep=r"\s+")
            fig = px.scatter(ref_df, x=ref_df.columns[2], y=ref_df.columns[3],
                            opacity=0.3, title="PC1 vs PC2 (1000G reference + target sample)")
            if len(target_df) > 0:
                fig.add_scatter(x=[target_df.iloc[0, 2]], y=[target_df.iloc[0, 3]],
                               mode="markers", marker=dict(size=20, color="red", symbol="star"),
                               name="You")
            st.plotly_chart(fig, width="stretch")
        except Exception:
            st.info("PCA data available but could not render plot.")

# ═══════════════════════════════════════════════════════════════════════════════
# ARCHAIC DNA (ancestry/deep_ancestry.json → "neanderthal")
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "🦴 Archaic DNA":
    st.title("🦴 Archaic DNA — Neanderthal Admixture")
    neand = safe_get(deep_anc, "neanderthal", default={})

    if not isinstance(neand, dict) or not neand:
        st.info("No archaic-admixture results. Run `python prs.py run --full --vcf sample.vcf.gz` "
                "(needs the Vindija/AADR reference bundle, see README).")
    else:
        pct = neand.get("percentage")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Neanderthal ancestry", f"{pct:.2f} %" if isinstance(pct, (int, float)) else "—")
        found, total = neand.get("snps_found"), neand.get("snps_total")
        c2.metric("Archaic SNPs found", f"{found:,} / {total:,}" if isinstance(found, int) and isinstance(total, int) else "—")
        c3.metric("Reliable", "✅ yes" if neand.get("reliable") else "⚠️ no")
        c4.metric("Closest population", neand.get("closest_population", "—"))

        ref = neand.get("reference", {})
        if isinstance(ref, dict) and ref:
            st.caption(f"Method: `{neand.get('method', '?')}` · Reference: {ref.get('source', '?')} "
                       f"({ref.get('coverage', '?')}, {ref.get('chromosomes', '?')} chromosomes)")

        comps = neand.get("population_comparisons", {})
        if isinstance(comps, dict) and comps:
            st.subheader("You vs. 1000 Genomes super-populations")
            rows = []
            for code, c in comps.items():
                if not isinstance(c, dict):
                    continue
                rows.append({"Population": f"{code} — {c.get('label', code)}",
                             "Population mean %": c.get("mean_pct"),
                             "You %": c.get("user_admix_pct", pct),
                             "z-score": c.get("z_score"),
                             "Percentile": c.get("percentile")})
            df_c = pd.DataFrame(rows)
            if not df_c.empty:
                long = df_c.melt(id_vars="Population", value_vars=["Population mean %", "You %"],
                                 var_name="Series", value_name="Neanderthal %")
                fig = px.bar(long, x="Population", y="Neanderthal %", color="Series", barmode="group",
                             color_discrete_map={"Population mean %": "#95a5a6", "You %": "#8e44ad"},
                             title="Neanderthal admixture: you vs. population means")
                st.plotly_chart(fig, width="stretch")
                st.dataframe(df_c, width="stretch", hide_index=True)

        with st.expander("Sharing statistics (raw)"):
            keep = {k: v for k, v in neand.items()
                    if k not in ("population_comparisons", "reference") and not isinstance(v, (dict, list))}
            st.json(keep)

        st.markdown("---")
        st.markdown("""
        **📖 Reading this page**
        - Present-day non-African genomes carry roughly 1.5–2.5 % Neanderthal ancestry; East Asians a little more than Europeans.
        - The estimate compares your genotype directly with the high-coverage Vindija Neanderthal genome; it is an *affinity* measure, not a count of introgressed segments.
        - "Reliable" is false when too few archaic-informative SNPs were found in your VCF (low coverage or a chr22-only reference).
        """)

# ═══════════════════════════════════════════════════════════════════════════════
# RAW DATA
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "📋 Raw Data":
    st.title("📋 Raw JSON Data")
    st.caption("Direct view of pipeline output files")

    files = {
        "PRS Results": "prs/PRS_RESULT.json",
        "PGS Catalog Calibration Report": "prs/pgs_scores/pgs_calibration_report.json",
        "Trait Recommendations (curated)": "data/trait_recommendations.json",
        "Ancestry Model": "science/ANCESTRY_MODEL.json",
        "ClinVar Pathogenic": "clinvar/clinvar_pathogenic_variants.json",
        "PharmGKB Drug Report": "pharmgkb/pharmgkb_drug_report.json",
        "Deep Ancestry": "ancestry/deep_ancestry.json",
        "Final Scientific Score": "FINAL_SCIENTIFIC_SCORE.json",
        "Validation Report": "science/global_validation_report.json",
    }

    selected = st.selectbox("Select file", list(files.keys()))
    if selected:
        data = load_json(files[selected])
        st.json(data)

st.sidebar.markdown("---")
st.sidebar.caption(f"BlueGen v{PIPELINE_VERSION}")
