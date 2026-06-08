#!/usr/bin/env python3
"""
🧬 BlueGen — Interactive Dashboard
   Run: streamlit run dashboard.py
   Reads all pipeline JSON outputs. No recomputation needed.
"""

import streamlit as st
import json
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

PIPELINE = Path(__file__).parent / "prs_research_pipeline"

# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data
def load_json(path):
    p = PIPELINE / path
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}

def safe_get(d, *keys, default="N/A"):
    for k in keys:
        if isinstance(d, dict):
            d = d.get(k, default)
        else:
            return default
    return d if d is not None else default

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════════════════════

st.sidebar.markdown("# 🧬 BlueGen")
st.sidebar.markdown("### Navigation")

page = st.sidebar.radio(
    "", [
        "📊 Overview",
        "🧬 PRS Results",
        "🔬 ClinVar Pathogenic",
        "💊 Pharmacogenomics",
        "🌍 Ancestry",
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
        st.metric("📦 Pipeline", "v1.0.0.0")
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
            st.plotly_chart(fig, use_container_width=True)
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
        st.plotly_chart(fig, use_container_width=True)

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
        st.plotly_chart(fig, use_container_width=True)

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
            st.dataframe(display_df, use_container_width=True, hide_index=True)

            # Bar chart
            fig = px.bar(df, x="population_zscore", y="trait",
                         title="PRS Z-Scores",
                         color="population_zscore",
                         color_continuous_scale=["#27ae60", "#f5f5f5", "#e74c3c"],
                         range_color=[-2, 2])
            fig.add_vline(x=1.0, line_dash="dash", line_color="#f39c12")
            fig.add_vline(x=-1.0, line_dash="dash", line_color="#27ae60")
            st.plotly_chart(fig, use_container_width=True)

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
            st.plotly_chart(fig, use_container_width=True)

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
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True,
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
            st.plotly_chart(fig, use_container_width=True)

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
        st.dataframe(df_drugs, use_container_width=True, hide_index=True)

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
            st.plotly_chart(fig, use_container_width=True)

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
            st.dataframe(df_subs[["code", "name", "description"]], use_container_width=True, hide_index=True)

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
            st.plotly_chart(fig, use_container_width=True)
        except Exception:
            st.info("PCA data available but could not render plot.")

# ═══════════════════════════════════════════════════════════════════════════════
# RAW DATA
# ═══════════════════════════════════════════════════════════════════════════════

elif page == "📋 Raw Data":
    st.title("📋 Raw JSON Data")
    st.caption("Direct view of pipeline output files")

    files = {
        "PRS Results": "prs/PRS_RESULT.json",
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
st.sidebar.caption("BlueGen v1.0.0")
