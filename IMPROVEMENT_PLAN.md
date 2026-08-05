# BlueGen — Plan de mejora (para implementar con Sonnet)

> **Objetivo real del proyecto:** convertir mi genoma (WGS VCF) en un **informe personal accionable** que use con mi nutricionista y para mejorar hábitos de vida. Este documento prioriza mejoras con ese objetivo en mente, no la "publicación científica".
>
> **Cómo usar este plan:** cada tarea tiene `Qué / Por qué / Dónde (archivo:línea) / Esfuerzo / Criterio de aceptación`. Están ordenadas por prioridad. Sonnet puede tomar cualquier bloque de forma independiente.
>
> Fecha: 2026-07-14. Rama: `main`.

---

## Diagnóstico en una frase

El entregable que de verdad importa —el **informe completo** (`comprehensive_report.py`, ~3.000 líneas)— es sólido pero puede **cubrir más rasgos y dar más contexto accionable**, y la lógica científica que lo sustenta (scoring, calibración) casi no tiene tests unitarios. La estrategia es **hacer el informe más rico y más fiable**: más cobertura genética, más interpretación accionable con evidencia citada, sin perder profundidad científica. La maquinaria de "publicación" (fases 7–10) se conserva; solo se ordena para que no estorbe al mantenimiento.

---

## 🔖 HANDOFF — estado actual y trabajo restante inmediato (leer primero)

> Contexto para retomar en una sesión nueva (p. ej. tras cambiar de modelo). H1-H4 del hilo del panel **ya están commiteados en `main`**. H5 está corriendo. El resto del documento (TIER 0–5) es el roadmap de más largo alcance.

**Ya hecho y verificado (2026-07-14/15):**
- ✅ Script de auditoría de posiciones: `prs_research_pipeline/scripts/setup/audit_snp_positions.py`. Reporte: `reports/snp_position_audit.json` (no commiteado, gitignored).
- ✅ **H1** (commit `666f19b`): panel curado commiteado — 108 posiciones GRCh37 corregidas, 10 rsID corruptos resueltos. `.gitignore` extendido (`/reports/`, `*.bak`) para que el HTML con genotipos reales y los backups de curación nunca puedan commitearse.
- ✅ **H2** (commit `7711b13`): QC literatura de las 4 filas marcadas. Bug real encontrado y arreglado: `rs2304672` (PER2) tenía `effect_allele=C` pero `risk_genotype=G/G` — el scoring usa `effect_allele`+`effect_direction` (no `risk_genotype`, que no lo lee ningún script), así que puntuaba el alelo equivocado; corregido a G/C. `rs13217795` (FOXO3) resultó ya estar bien tras encontrar la fuente primaria real (Willcox 2008, PMID 18765803). `rs590787` (RHD) bajado a evidencia D/baja confianza — su único uso documentado es un estudio chino de weak-D/DEL, no tipaje general RhD+/−. Las 4 filas tenían PMIDs falsos (citaban papers no relacionados) heredados de datos placeholder — corregidos.
- ✅ **H3** (commit `f1960d5`): añadidas 7 filas nuevas de la sub-tarea 1.2.B — FUT2 (`rs602662`), HFE H63D (`rs1799945`), ADORA2A (`rs5751876`), HNMT (`rs11558538`), COL1A1 (`rs1800012`), CLOCK (`rs1801260`), ACE (`rs4343`). Cada una verificada independientemente contra dbSNP GRCh37 + cita primaria confirmada en PubMed (no copiadas tal cual de la tabla del plan).
- ✅ **H4** (commit `94c373f`): `tests/test_snp_positions.py` (6 tests) wireado en `.github/workflows/test.yml` como step propio ("SNP position regression"). Escribir el check de duplicados destapó 3 bugs preexistentes — `rs1801131` (MTHFR), `rs1695` (GSTP1), `rs1048943` (CYP1A1) tenían cada uno DOS filas en la misma `trait_category`, duplicando el peso de ese SNP en la suma del PRS. Verificadas contra dbSNP/literatura, eliminada la fila incorrecta de cada par, corregidos 3 PMIDs falsos más.
- ✅ Panel tras H1-H4: 192 filas, `195/195 → 192/192 MATCH` en la auditoría de posiciones (0 errores).

**Hallazgo transversal de H2/H4, confirmado y resuelto a mayor escala (no queda pendiente):** el patrón de `effect_allele`/PMIDs corruptos en filas sueltas resultó ser sistémico. Se auditó el resto del panel y se corrigió en tres commits posteriores al handoff original:
- ✅ **TIER 0.1** (commit `e96c1ab`): versión del pipeline unificada en una sola fuente de verdad (`constants.py::PIPELINE_VERSION`, ahora `2.0.0`); 9 archivos que hardcodeaban su propia versión migrados; test de regresión añadido.
- ✅ **Fix de calibración** (commit `4a41436`): bug z-score/percentile en `population_calibrate_v2.py` forzaba `risk_category="medium"` para varios rasgos sin importar el z real. Corregido; score científico subió de 74.1 (NEEDS_REVISION) a 86.7 (RESEARCH_GRADE).
- ✅ **Fix de strand-orientation** (commit `754a972`): 34 filas con `effect_allele`/`reference_allele` en notación gene-relative (hebra minus) en vez de forward-strand GRCh37 — PLINK las descartaba en silencio. Corregidas 26 mecánicamente + 6 con problema adicional de alelo verificadas contra literatura primaria; 2 filas eliminadas por datos placeholder/gen equivocado. Test de regresión añadido (`tests/test_allele_strand_consistency.py`). Resultado: 23→25 rasgos puntúan sobre el genotipo real del usuario, `well_calibrated` 11→22/23.
- ✅ **Fix GWAS consortium + PGS percentile** (commit `15a8a4e`): alias de sub-rasgo mal mapeados en la validación de consorcios GWAS (4/17→16/17 passed) y un `dir()` sin argumentos que congelaba el percentil de PGS en 50% para las 30 puntuaciones. Ambos corregidos.
- ✅ (commit `d46164a`): `FINAL_SCIENTIFIC_SCORE.json`/`PUBLICATION_LOCK.json` regenerados tras los fixes anteriores — score final actual: **84.4 (RESEARCH_GRADE)** (bajó un poco desde 86.7 porque Population Portability ahora refleja con más precisión la brecha real EUR-vs-otras-ancestrías del panel corregido, no una regresión).
- **Panel actual (verificado 2026-08-04):** `prs_research_pipeline/data/snp_database_annotated.csv` tiene **190 filas / 172 rsID únicos / 59 trait_category** (bajó de 192 por las 2 filas eliminadas en el fix de strand).

- ✅ **H5** (2026-08-04): informe HTML regenerado desde el VCF genome-wide real (`prs.py run --full`, 1511s primera corrida). **25/59 rasgos puntúan** sobre el genotipo real (23→25 según lo previsto por el fix de strand). Confirmado: Lactose intolerance (z=+5.39, GOOD calibración R²=0.882) e Histamine intolerance (z=+1.43, GOOD R²=0.998) muestran los `risk_category` recalculados correctamente. De las 7 filas nuevas de H3, solo ACE (`rs4343`) tiene una llamada de variante en el `qc_filtered.bim` real del usuario — las otras 6 (FUT2, HFE H63D, ADORA2A, HNMT, COL1A1, CLOCK) no están cubiertas en este WGS concreto (misma categoría de limitación que rs671/alcohol flush: posición correcta, sample sin esa región genotipada — no es un bug del panel).

**Tres bugs adicionales encontrados y arreglados verificando el informe post-H5 (2026-08-04), cada uno con test de regresión, verificados con una segunda corrida completa del pipeline (exit 0):**
- ✅ **Referencia huérfana al integrity score** (commit `43837cf`): `41_unified_report_engine.py` y `45_publication_evidence_pack.py` leían por defecto `science/scientific_integrity_score.json`, archivo que solo genera `34_scientific_integrity_score.py` — script registrado en `ROUTES` de `prs.py` pero **nunca invocado** por `run --full` (el score real vive en `FINAL_SCIENTIFIC_SCORE.json`, escrito después por el script 46). Efecto visible: `SCIENTIFIC_MANUSCRIPT_EN/ES.md` decía *"Scientific Integrity Score of 0/100 (Unknown)"* mientras el HTML sí mostraba 84.4 bien (lee el JSON correcto directo). Arreglado apuntando ambos scripts a `FINAL_SCIENTIFIC_SCORE.json` y reordenando `prs.py` para correr `final_score` antes de `report_engine`/`evidence_pack`.
- ✅ **PRS_CORE congelado desde junio, nunca refrescado** (commit `43837cf`): `36_prs_core.py::load_or_create()` reusaba `science/prs_core_definition.json` para siempre una vez creado, sin comparar contra el panel actual. Frozen el 2026-06-03 con 109 SNPs/10 traits, nunca actualizado pese a que H1-H4 llevó el panel a 190/59. Arreglado con un hash SHA-256 del CSV (`source_csv_hash`): si el CSV cambió (o el campo no existe, para archivos legacy), se re-congela. 4 tests nuevos en `tests/test_prs_core_freeze.py`. Verificado: manuscrito ahora dice correctamente "190 SNPs across 59 traits".
- ✅ **Banner "Critical Findings" falso + LOCK-010 con condición imposible** (commit `8faf490`): `43_adversarial_prs_validation.py` metía en `critical_findings` **todos** los tests de severidad CRITICAL sin filtrar por si realmente fallaron (`not r.is_robust`) — el HTML mostraba un banner rojo listando 2 tests que en realidad pasaron (✅ Robust en la misma tabla, justo debajo). Y `45_publication_evidence_pack.py` comparaba el tamaño de un catálogo estático de 18 modos de fallo documentados (`failure_mode_map.json`, siempre 7 CRITICAL, algunos heredados de la arquitectura pre-Phase-6 como "chr22-only bias") contra "expected: 0" — una condición que **nunca podía pasar por diseño**. Arreglado: filtro `not r.is_robust` en el primero; el segundo ahora lee `adversarial.critical_findings` (ya corregido) en vez del catálogo estático. Resultado: `LOCK-010` sigue fallando (9/10), pero ahora por una razón real y accionable — `LD_DISRUPT_SEVERE` (VIF=3.0x bajo disrupción de LD severa), no un fantasma. 2 tests nuevos en `tests/test_adversarial_validation.py` y `tests/test_publication_evidence_pack.py`.

- ✅ **`LD_DISRUPT_SEVERE` investigado (commit `a217980`, 2026-08-04): era un cuarto bug, no una vulnerabilidad real.** `_test_ld_disruption` calculaba `vif = (variance × inflation) / variance` — se simplifica algebraicamente a `vif = inflation` siempre, sin importar los datos reales. Con `robust = vif < 2.0`, MODERATE/SEVERE estaban condenados a "vulnerable" en cualquier corrida, para siempre — no medía nada. Reemplazado por una simulación real (perturbar los z-scores reales con ruido calibrado al nivel de inflación objetivo, medir si el ranking de rasgos sobrevive), misma metodología que los otros 3 stress-tests del archivo. Verificado como no-tautológico con tests que prueban que el resultado depende de los datos de entrada. **Con datos reales, el resultado empeoró** (los 3 niveles fallan ahora, no solo SEVERE): Adversarial Robustness 85→77, score final 84.4→**82.9 (sigue RESEARCH_GRADE)**. La calibración del ruido (cuánta perturbación = "leve" vs "severa") es una decisión de diseño, no un hecho matemático — revisada y aceptada tal cual, no re-calibrada.

**Trabajo restante inmediato:**

- [ ] Pushear los commits locales a `origin/main` (varios acumulados sin pushear esta sesión).
- [ ] Si se quiere refinar más la validación adversarial: revisar si la calibración del ruido en `_test_ld_disruption` (factor idiosincrático `extra_sd*0.3`, escala percentil `×15` heredada de `_test_population_shift`) es la más adecuada, o solo documentar que MILD ahora también falla como hallazgo genuino a explicar en el informe.

**Después del handoff:** seguir con el `Orden sugerido de ataque` (al final del documento) — TIER 0, luego enriquecer el informe (1.1/1.4/1.5), etc.

---

## TIER 0 — Rápidas, alto valor, bajo riesgo (empezar aquí)

### 0.1 Unificar el versionado — ✅ **Hecho (commit `e96c1ab`, verificado 2026-08-05)**
- **Qué:** hay tres versiones distintas: `README.md:11` dice `1.1.0`, `prs_research_pipeline/config.yaml:13` dice `2.0.0`, y `prs.py:345` imprime `v1.0.0` en el header del log.
- **Por qué:** confunde qué versión generó un informe; rompe la trazabilidad que el propio proyecto dice ofrecer.
- **Dónde:** `README.md:11`, `prs_research_pipeline/config.yaml:13`, `prs.py:345`, `prs_research_pipeline/README.md` (tabla de versiones ~L345).
- **Cómo:** definir la versión en **un solo sitio** (p. ej. `prs_research_pipeline/scripts/utils/constants.py`) y que README/config/log la importen o la referencien. Añadir un check en CI que falle si divergen.
- **Esfuerzo:** 1–2 h. **Criterio:** `grep -rn "version" ` muestra una sola fuente de verdad.
- **Estado:** ya estaba implementado de una sesión previa a este handoff — `PIPELINE_VERSION = "2.0.0"` vive en `prs_research_pipeline/scripts/utils/constants.py`, README (badge + tabla), `config.yaml` y `prs.py` lo referencian. `tests/test_version_consistency.py` (5 tests) verifica que las 3 fuentes coincidan y que ningún script tenga un literal `pipeline_version` hardcodeado; corre en CI vía `pytest tests/ -v`. Verificado hoy: 5/5 tests passing.

### 0.2 Arreglar la ruta hardcodeada de referencia chr22 — ✅ **Hecho (2026-08-05)**
- **Qué:** `config.yaml:110-111` y `prs.py:367` fijan `reference/1000G/ALL.chr22...vcf.gz` como fallback. Es una ruta concreta a un único cromosoma.
- **Por qué:** frágil; si falta el genome-wide, el pipeline cae a chr22 silenciosamente y calcula ancestría/calibración con datos parciales sin avisar de forma prominente en el informe.
- **Dónde:** `prs.py:365-380`, `config.yaml:106-121`.
- **Cómo:** que el fallback a chr22 escriba un flag `reference_coverage: "chr22_only"` en `PRS_RESULT.json` y que el informe lo muestre como advertencia visible (no solo en logs).
- **Esfuerzo:** 2–3 h. **Criterio:** informe generado con solo chr22 muestra un banner de "cobertura de referencia parcial".
- **Estado:** `prs.py` ya calculaba `use_full_ref` (L376) pero nunca lo propagaba más allá de un `warn()` en el log. Ahora se pasa como `--reference-coverage genome_wide|chr22_only` a la stage `prs_result` (`sss/37_prs_result_unified.py`, nuevo flag CLI), que lo escribe en `PRS_RESULT.json::metadata.reference_coverage`. `comprehensive_report.py` añade `reference_coverage_banner()` (bilingüe EN/ES) que se muestra al principio de la página — antes de cualquier sección, no enterrado — cuando el valor es `chr22_only`; con `genome_wide` no renderiza nada.
  - **Hallazgo de paso (relacionado con 0.1):** `comprehensive_report.py` tenía `"BlueGen v10.0"` hardcodeado en dos sitios (header + footer), divergiendo de `PIPELINE_VERSION="2.0.0"` — el mismo bug que 0.1 se supone que evitaba, pero el test existente solo vigilaba asignaciones `pipeline_version=...`, no texto libre. Corregido a usar `PIPELINE_VERSION` (import añadido); test nuevo `test_no_hardcoded_bluegen_v_string_anywhere` en `tests/test_version_consistency.py` (ahora 6/6) evita que vuelva a pasar en cualquier script.
  - 151/151 tests pasan (`venv/bin/python -m pytest tests/ -q`).

### 0.3 `.gitignore` / privacidad de datos genómicos — ✅ **Hecho (2026-08-04)**
- **Qué:** confirmar que ningún dato personal identificable (VCF, BAM, FASTQ, informes con genotipos) puede colarse a git. Ya hay reglas, pero conviene un test.
- **Por qué:** el genoma es el PII más sensible que existe; un `git add .` accidental es irreversible una vez pusheado.
- **Dónde:** `.gitignore`, `.github/workflows/test.yml` (añadir un job).
- **Cómo:** añadir a CI un check que falle si el commit incluye `*.vcf.gz`, `*.bam`, `*.fq.gz`, `reports/*.html`, o rutas `/Users/`. Ya existe el check de `/Users/` en `test.yml:30` — extenderlo a extensiones de datos.
- **Esfuerzo:** 1 h. **Criterio:** CI falla si se intenta commitear un VCF de prueba.
- **Estado:** step "Check for accidentally committed genomic/personal data" añadido al job `lint` de `test.yml` — falla si `git ls-files` encuentra `.vcf.gz`/`.vcf.gz.tbi`/`.bam`/`.bam.bai`/`.cram`/`.fq.gz`/`.fastq.gz`/`.bak` o cualquier `reports/*.html` trackeado. Verificado localmente contra el árbol actual (0 offenders).

### 0.4 `requirements.txt` con versiones fijadas (lock) — ✅ **Hecho (2026-08-04)**
- **Qué:** las deps usan `>=` abierto (`prs_research_pipeline/requirements.txt`). `weasyprint>=60` arrastra libs de sistema (cairo/pango) frágiles en macOS.
- **Por qué:** reproducibilidad — el proyecto presume de "deterministic seeds, SHA-256" pero el entorno Python no está pinneado.
- **Dónde:** `prs_research_pipeline/requirements.txt`.
- **Cómo:** generar un `requirements.lock` con `pip freeze` del venv que funciona hoy, y documentar que WeasyPrint es opcional (el HTML es el entregable primario; el PDF, secundario).
- **Esfuerzo:** 1 h. **Criterio:** `pip install -r requirements.lock` reproduce el entorno actual.
- **Estado:** `prs_research_pipeline/requirements.lock` generado desde el venv que corre H5 ahora mismo (Python 3.14.5, macOS). `requirements.txt` apunta a él. Nota de WeasyPrint opcional incluida en el header del lock.

---

## TIER 1 — El informe: hacerlo más completo y más accionable (máxima prioridad de producto)

> **Principio:** el informe se mantiene completo. Todo lo que sigue es **aditivo** — más cobertura, más contexto, más acciones — sin quitar profundidad científica.

### 1.1 Añadir un "Resumen ejecutivo accionable" al principio (sin quitar nada) — ✅ **Hecho (commit `d16fc3a`, 2026-08-04)**
- **Qué:** una primera sección que destile los hallazgos más relevantes en lenguaje claro para el nutricionista/uno mismo, **encima** del informe científico completo que ya existe.
- **Por qué:** el informe es potente pero denso; una portada con "tus 10 hallazgos clave + qué hacer" multiplica su utilidad práctica sin perder el detalle de abajo.
- **Dónde:** `prs_research_pipeline/scripts/publication/comprehensive_report.py` (añadir sección al inicio del render), consume `prs/PRS_RESULT.json`, `clinvar/…`, `pharmgkb/…`, `ancestry/…`.
- **Cómo:** priorizar rasgos por (magnitud del z-score / percentil extremo) × (nivel de evidencia) × (accionabilidad). Cada hallazgo: título en claro → qué significa → **acción concreta** → nivel de evidencia. Con enlaces ancla a la sección científica completa.
- **Esfuerzo:** 2–3 días. **Criterio:** el informe abre con un resumen de hallazgos priorizados y el resto del contenido sigue intacto debajo.
- **Estado:** nueva sección "🔍 Your Top Findings" / "Tus Hallazgos Principales" al inicio del HTML (EN/ES), abierta por defecto. Prioriza por `|z-score| × evidencia × confianza` (no hay campo de "accionabilidad" real todavía — eso es 1.4). Cada card: rasgo, badge de riesgo, badge de evidencia A-D, una frase en lenguaje claro (qué se midió, base estadística), y enlace ancla a la fila exacta en la tabla PRS completa (`id="trait-<slug>"` nuevo en cada fila). **Deliberadamente NO fabrica una recomendación específica por rasgo** — el campo `recommendation_en/es` se lee si existe (lo llenará 1.4), si no hay, muestra un texto honesto de fallback ("consultar con un profesional, ver detalle abajo"). Evita repetir el patrón de PMIDs/notas fabricadas que esta sesión ya encontró y corrigió varias veces en el panel de SNPs. 12 tests nuevos.

### 1.2 Ampliar cobertura y **auditar posiciones** del panel

> **Corrección importante (verificado 2026-07-14):** el panel actual ya es MUCHO más amplio de lo que parecía — cubre **~60 categorías de rasgo**, no ~10. Ya incluye alcohol flush (`ALDH2` rs671), dependencia (`ADH1B` rs1229984), celíaco DQ2.5/DQ8 con **tag-SNPs** (`rs2187668`, `rs7454108` — **no hace falta imputar HLA**), hemocromatosis `HFE` C282Y, histamina (`AOC1`/DAO), presión (`AGT`), cronotipo, vitaminas A/C/D/E, urato, ABO, pigmentación, etc. Por tanto **la mayor palanca NO es añadir genes a lo loco, sino: (a) auditar/corregir posiciones, y (b) dar profundidad a rasgos hoy sostenidos por 1 solo SNP.**

- **Hallazgo concreto que obliga a priorizar la auditoría:** `rs671` (ALDH2, alcohol flush) está en el CSV como `chr12:112240098` pero la posición **GRCh37 correcta es `chr12:112241766`** (verificado en Ensembl GRCh37). Con la posición mal, ese SNP **no casa en PLINK** y el rasgo sale vacío/erróneo. Casi seguro hay más casos así → **la auditoría de posiciones (tarea 3.3) debe ir ANTES de añadir nada.**

**Sub-tarea A — Auditar/curar posiciones existentes (prioritario). ✅ COMPLETADA (2026-07-14, por Sonnet).**

> **Estado final:** panel curado, auditoría `188 MATCH / 0 POS ERROR / 0 CHROM ERROR / 0 NOT FOUND`. 108 posiciones corregidas; 191 → 188 filas (eliminadas 3: duplicado MMP3, duplicado SLC17A1, y `rs41380247` que no existe en dbSNP). Cambios en el working tree, **sin commitear**, `.bak` preservado. Los 10 rsID corruptos se sustituyeron por el SNP real de cada gen (rs5400/SLC2A2, rs10156191/AOC1, rs1165196/SLC17A1, rs590787/RHD, rs12203592/IRF4, rs2304672/PER2, rs2736100/TERT, rs13217795/FOXO3). Cada fila corregida lleva nota `Corrected 2026-07-14:` explicando el cambio.
>
> **Seguimientos pendientes (QC humano antes de confiar en dirección/peso):**
> - `rs590787` (etiqueta ABO pero locus real **RHD**): dirección RhD+/− no re-derivada; revisar `effect_direction`/`risk_genotype`/`weight`, y considerar renombrar `gene` a RHD.
> - `rs2736100` (TERT): asociación bidireccional (telómero largo ↔ riesgo variable de cáncer); revisar convención +/−.
> - `rs2304672` (PER2) y `rs13217795` (FOXO3): dirección de alelo no re-derivada letra a letra; `evidence_level` posiblemente generoso.
> - **FUT2/B12:** la fila `rs602662` se eliminó (era FUT2 mal etiquetado como SLC17A1/urato). Añadir FUT2 correctamente (chr19:49206985) queda para la **sub-tarea 1.2.B**.
> - **Verificación ✅ (2026-07-14):** prueba A/B del scoring (`prs_plink_score.py`) sobre el mismo `qc_filtered` bfile con CSV viejo vs curado → **12 → 22 rasgos, 21 → 40 SNPs casados** (casi ×2). Nuevos rasgos que ahora puntúan: Histamina, Hierro, Pigmentación, Vitamina A, Vitamina D binding, Tendinopatía, Colorrectal, Cerumen, Rendimiento muscular, Dolor, Inflamación CRP. Nota: "Alcohol flush" (rs671) sigue sin puntuar, pero **no por la posición** (ya corregida) sino porque rs671 no está genotipado en este bfile concreto (ausente/filtrado en QC) — limitación del dato del sample, no del panel.


Script entregado: **`scripts/setup/audit_snp_positions.py`** (batch Ensembl GRCh37, separa errores en dos clases, `--fix` solo aplica los seguros, hace `.bak`). Reporte: **`reports/snp_position_audit.json`**.

**Resultado sobre `snp_database_annotated.csv` (191 filas, 170 rsIDs únicos):**

| Clase | Nº | Qué es | Acción |
|---|---|---|---|
| ✅ MATCH | 72 | posición correcta | ninguna |
| ⚠️ POS ERROR | **108** | mismo cromosoma, posición desplazada (el rsID es fiable, la pos está mal) | auto-fix seguro (`--fix`) |
| 🛑 CHROM ERROR | **10** | el rsID mapea a **otro cromosoma** que el gen etiquetado → **el rsID es el campo corrupto** (copy-paste). La posición del CSV sí cae en el gen | **NO** auto-fix; curar rsID por posición |
| ❓ NOT FOUND | 1 | `rs41380247` (rsID fusionado/retirado o indel) | revisar en dbSNP |

- **Impacto real (confirmado):** el scoring casa por `chrom:pos` (`scripts/prs/prs_plink_score.py:119` → `vid = f"{chrom}:{pos}"`, solo puntúa `if vid in bim_ids`). Por tanto **las 118 filas con posición errónea NO casan y se descartan en silencio: ~62% del panel no está entrando en el informe.** Esto explica rasgos vacíos (p. ej. alcohol flush: `rs671` estaba en `chr12:112240098`, correcto `chr12:112241766`).

- **Los 10 CHROM ERROR (curar rsID, NO mover posición):**

| rsID en CSV | etiqueta (gen/rasgo) | pos CSV (cae en el gen) | el rsID realmente está en |
|---|---|---|---|
| rs11605924 | SLC2A2 / Glucosa | chr3:170733088 | chr11 (es CRY2) |
| rs10177833 | AOC1 / Histamina | chr7:150543321 | chr2 |
| rs70991108 | SLC17A1 / Urato | chr6:25821380 | chr5 |
| rs1476413 | MMP3 / Tendinopatía | chr11:102834698 | chr1 (duplicado de la fila MTHFR) |
| rs1053878 | ABO / Grupo Rh | chr1:25592883 | chr9 |
| rs602662 | SLC17A1 / Urato | chr6:25821380 | chr19 (es FUT2/B12) |
| rs1015362 | IRF4 / Pigmentación | chr6:396461 | chr20 (región ASIP) |
| rs228697 | PER2 / Cronotipo | chr2:239164396 | chr1 (es PER3) |
| rs7689424 | TERT / Telómeros | chr5:1254653 | chr4 |
| rs10224002 | FOXO3 / Longevidad | chr6:108960955 | chr7 |

  → Para cada uno: reverse-lookup por posición en dbSNP para hallar el rsID correcto del gen, o corregir la etiqueta. Nota: las dos filas "SLC17A1/Urato" comparten la MISMA posición placeholder (chr6:25821380) — son duplicados a resolver.

- **Cómo aplicar el fix (workflow, para Sonnet):**
  1. `python scripts/setup/audit_snp_positions.py --fix` → corrige las 108 POS ERROR (crea `.bak`), deja intactas las 10 CHROM ERROR.
  2. Curar manualmente las 10 CHROM ERROR + `rs41380247`.
  3. Re-ejecutar el pipeline y **verificar que los rasgos antes vacíos ahora puntúan** (criterio de aceptación).
  4. Añadir un test (TIER 4) que falle si algún rsID diverge de dbSNP/Ensembl.
- **Nota de seguridad (por qué no se auto-corrigió ya):** es un fichero de datos de salud y la auditoría demuestra que el propio rsID a veces está mal; el fix seguro (108) descansa en confiar en el rsID, así que conviene revisar el diff y curar las 10 antes de commitear. El `.bak` y git lo hacen reversible.
- **Criterio:** 0 divergencias tras curar; rasgos afectados (empezando por alcohol flush) aparecen en el informe.

**Sub-tarea B — Rellenar huecos puntuales de alto valor (posiciones ya verificadas).**
Todas verificadas hoy en Ensembl GRCh37. El **alelo de efecto/peso/evidencia se curan en 1.4** (aquí van marcados como *candidato*):

| rsID | gen | rasgo | chr:pos (GRCh37) | alelos | efecto (candidato) | nota |
|---|---|---|---|---|---|---|
| rs1799945 | HFE | Hierro/hemocromatosis (H63D) | chr6:26091179 | C/G/T | G | 2ª variante HFE, **falta** (solo está C282Y) |
| rs5751876 | ADORA2A | Cafeína — ansiedad/sueño | chr22:24837301 | T/C | T | complementa CYP1A2 (velocidad ≠ sensibilidad) |
| rs602662 | FUT2 | Vitamina B12 / estado secretor | chr19:49206985 | G/A | A | B12 sérica; microbiota; **ausente** |
| rs11558538 | HNMT | Intolerancia a histamina | chr2:138759649 | C/T | T | Thr105Ile; complementa AOC1/DAO |
| rs1800012 | COL1A1 | Densidad ósea (Sp1) | chr17:48277749 | C/A | A | complementa el panel óseo |
| rs1801260 | CLOCK | Cronotipo (3111T/C) | chr4:56301369 | A/G | G | vespertino; complementa MTNR1B/PER |
| rs4343 | ACE | Presión arterial / sal (proxy I/D) | chr17:61566031 | G/A | G | tag del alelo D; complementa AGT |

- **Dónde:** añadir filas a `prs_research_pipeline/data/snp_database_annotated.csv`; declarar categorías nuevas si aplica en `config.yaml:78-89`; genes de interés en `config.yaml:140-169`; interpretación bilingüe en `scripts/validation/bilingual_interpretation.py`.
- **Nota HFE C282Y:** de paso, verificar la posición de `rs1800562` ya presente (auditoría A).

**Sub-tarea C — Dar profundidad a rasgos con 1 solo SNP.**
- **Qué:** ~30 rasgos hoy dependen de un único SNP (ver `cut -d, -f3 | uniq -c`): p. ej. "Sleep duration", "Vitamin C", "Vitamin E", "Cognitive function", inflamación. Un solo SNP es ruidoso.
- **Cómo:** para los rasgos accionables en nutrición, añadir 2–4 SNPs bien establecidos por rasgo (con posición verificada) para un score más estable. Priorizar sobre añadir rasgos nuevos.
- **Criterio:** ningún rasgo accionable del informe descansa en un único SNP salvo justificación (variante mayor tipo APOE/HFE).

- **Esfuerzo total 1.2:** 4–6 días (mayormente curación + verificación). **Criterio global:** informe cubre los nuevos rasgos con genotipo+interpretación+acción, posiciones auditadas, y el test E2E (4.2) los detecta.

### 1.3 Integrar más fuentes de datos ya disponibles en el informe
- **Qué:** el pipeline ya calcula/descarga muchas cosas que podrían aparecer con más contexto en el informe: PharmGKB (fármaco→gen), ClinVar+MedGen (variantes patogénicas con definición de enfermedad), ancestría profunda (mtDNA/Y-DNA/Neandertal), y 45 PGS aún sin exponer (ver 3.1).
- **Por qué:** ya se paga el coste de calcularlos; falta presentarlos de forma accionable e integrada.
- **Dónde:** `comprehensive_report.py` (secciones de PharmGKB/ClinVar/ancestría), `clinical/pharmgkb_annotator.py`, `clinical/clinvar_annotator.py`, `clinical/medgen_enrich.py`.
- **Cómo:** para farmacogenética, añadir "qué implica" por fármaco (dosis/precaución según CPIC). Para ClinVar, agrupar por sistema/órgano y enlazar la definición MedGen. Sección de ancestría con contexto divulgativo.
- **Esfuerzo:** 2–3 días. **Criterio:** cada módulo tiene una sub-sección con lectura accionable, no solo tablas crudas.

### 1.4 Traducir cada score a recomendación accionable **con evidencia citada** — 🟡 **Parcial: 10/59 rasgos hecho (2026-08-04)**
- **Qué:** las interpretaciones viven en `config.yaml:140-169` (9 genes) y `bilingual_interpretation.py`. Faltan **acciones** (qué comer/evitar/monitorizar) con su cita y nivel de evidencia.
- **Por qué:** es lo que hace el informe útil para el nutricionista y fiable para decisiones. Añade valor sin quitar contenido.
- **Dónde:** `scripts/validation/bilingual_interpretation.py`, `config.yaml:140-169`, `data/snp_database_annotated.csv`.
- **Cómo:** cada rasgo con campos `recommendation_en/es` + `evidence_level` (A/B/C) + `reference` (PMID/URL). En el informe, mostrar la acción destacada y el nivel de evidencia como badge.
- **Esfuerzo:** 2–3 días (curación). **Criterio:** cada rasgo del informe tiene acción + nivel de evidencia + cita.
- **Estado:** decisión deliberada de **no** apurar los 59 rasgos de una — cada recomendación requiere verificar una fuente real (PubMed/NIH), y esta sesión ya encontró PMIDs fabricados en el panel al apurar curación antes. Se hicieron los **10 rasgos más accionables en nutrición** de los que puntúan sobre datos reales: Lactose intolerance, Caffeine metabolism, Glucose metabolism, Folate & methylation, Vitamin D metabolism, Lipid metabolism, Histamine intolerance, Iron levels, Blood pressure, Inflammation (CRP levels).
  - **Dónde quedó implementado** (distinto de "Dónde" arriba — diseño más aditivo/desacoplado): `prs_research_pipeline/data/trait_recommendations.json` (nuevo archivo, no toca `config.yaml` ni el CSV del panel — cada entrada: `recommendation_en/es`, `evidence_level`, `reference`). `comprehensive_report.py` lo carga y lo pasa a `build_top_findings()` (de 1.1), que solo muestra la recomendación curada cuando `risk_category == "high"` (el texto está escrito para la dirección de riesgo elevado; mostrarlo en un hallazgo protector/promedio sería incorrecto). Sin cobertura curada o sin risk="high" → cae al fallback honesto de 1.1.
  - **Cada recomendación fue verificada contra una fuente real** antes de escribirse (PMID de PubMed o ficha de NIDDK/NIH ODS) — no reescribe conocimiento genérico, cita el hallazgo específico (p. ej. Iron levels queda deliberadamente sin recomendación dietética directa — combina HFE con TMPRSS6 en direcciones opuestas, así que remite a análisis de sangre en vez de arriesgar un consejo dietético incorrecto).
  - Tests: `TestTraitRecommendationsData` (cada trait curado existe en el panel real, cada entrada tiene ambos idiomas + evidencia + referencia) + tests de gating en `TestBuildTopFindings`.
  - **Pendiente:** ~49 rasgos restantes, misma metodología (WebSearch + verificación PMID por rasgo antes de escribir), en sesiones futuras.

### 1.5 Contexto honesto de base (aditivo, refuerza credibilidad) — ✅ **Hecho (mayormente ya existía; brecha cerrada en commit `b589fcb`, 2026-08-04)**
- **Qué:** mostrar en cada rasgo cuántos SNPs lo sustentan y su límite; distinguir "nutrigenética de panel pequeño (bien establecida)" de "PRS de enfermedad compleja (limitado)". El disclaimer de `config.yaml:180-196` es bueno pero está enterrado.
- **Por qué:** un informe completo y **bien calibrado en confianza** es mejor producto que uno que aparenta certeza uniforme; también protege ante decisiones de salud.
- **Dónde:** `config.yaml:180-196`, `comprehensive_report.py`.
- **Cómo:** badge de "n SNPs" y "nivel de evidencia" por rasgo; disclaimer accesible pero no intrusivo (colapsable). Esto acompaña, no reemplaza, la profundidad.
- **Esfuerzo:** 1 día. **Criterio:** cada rasgo muestra su base de evidencia y n de SNPs.
- **Estado:** al auditar el código (2026-08-04) resultó que casi todo esto ya estaba construido en una sesión anterior a este plan (CHANGELOG 1.1.1/1.2.0): barra de cobertura de SNPs, estrellas de confianza, badge de calibración, trust tier T1/T2/T3, badges de limitación — todo por rasgo en la tabla PRS. Lo único genuinamente pendiente era el disclaimer enterrado al final del informe — se agregó un `<details>` colapsable (sin JS) con el mismo texto justo debajo de Top Findings, arriba del todo.

### 1.6 Refactorizar `comprehensive_report.py` (2.998 líneas) para poder crecer
- **Qué:** un solo archivo con HTML, CSS, datos y lógica mezclados. Al añadir secciones (1.1–1.5) crecerá; conviene ordenarlo antes/durante.
- **Por qué:** es EL entregable; para enriquecerlo con seguridad necesita estructura y tests de render.
- **Dónde:** `prs_research_pipeline/scripts/publication/comprehensive_report.py`, plantillas en `prs_research_pipeline/templates/`.
- **Cómo:** extraer a Jinja2 (ya es dependencia), separar (a) carga de datos, (b) interpretación, (c) render. Módulo `report/` con `data_loader.py`, `interpretations.py`, `render.py`. Migrar por secciones.
- **Esfuerzo:** 2–3 días. **Criterio:** HTML equivalente (diff justificado) tras el refactor; tests de render por sección; añadir una sección nueva es trivial.

---

## TIER 2 — Arquitectura y mantenibilidad

### 2.1 Extraer una librería importable del scoring (dejar de shell-out todo)
- **Qué:** `prs.py` orquesta ~40 scripts vía `subprocess` (`prs.py:147-205`, registro en `prs.py:73-138`), pasando estado por archivos JSON en disco. No hay librería compartida de scoring; cada script re-lee/re-parsea.
- **Por qué:** imposible de testear unitariamente, lento (arranque de Python × 40), y frágil (los nombres numerados `01_`, `06_`, `36_` acoplan orden a nombre de archivo).
- **Dónde:** `prs.py:73-138` (ROUTES), scripts en `prs_research_pipeline/scripts/{prs,validation,sss,publication}/`.
- **Cómo:** crear `prs_research_pipeline/bluegen/` como paquete: `io.py` (carga/validación de artefactos), `scoring.py`, `calibration.py`, `ancestry.py`. Los scripts numerados pasan a ser thin wrappers que importan el paquete. Migrar de forma incremental empezando por la ruta crítica (Stages F–H).
- **Esfuerzo:** 3–5 días (incremental). **Criterio:** el scoring PRS se puede llamar como función Python y tiene tests unitarios sin subprocess.

### 2.2 Un esquema/contrato para los artefactos JSON
- **Qué:** decenas de JSON intermedios (`PRS_RESULT.json`, `ANCESTRY_MODEL.json`, etc., ver `prs_research_pipeline/README.md:165-188`) sin esquema formal. El test suite solo comprueba "existe y no hay tipos numpy".
- **Por qué:** cambios de formato rompen consumidores (informe, dashboard) silenciosamente.
- **Dónde:** `prs_research_pipeline/scripts/utils/test_suite.py`, todos los productores de JSON.
- **Cómo:** definir dataclasses o `pydantic`/`jsonschema` para los 4–6 JSON críticos y validar al escribir y al leer.
- **Esfuerzo:** 2 días. **Criterio:** escribir un JSON con forma inválida falla en el productor, no en el consumidor.

### 2.3 Ordenar (no quitar) la maquinaria de validación/publicación
- **Qué:** Fases 7–10 (SSST, adversarial, failure_map, publication_lock, evidence_pack, manuscritos) — ~35 scripts numerados. Se **conservan** (aportan al informe completo y al rigor), pero conviene aislar las que dan valor de las que son coste puro de mantenimiento.
- **Por qué:** el informe sigue siendo completo; el objetivo aquí es solo bajar el coste de mantenimiento y el tiempo de ejecución sin perder capacidades.
- **Dónde:** `scripts/sss/`, `scripts/publication/43-47*`, `scripts/validation/18,20-22,30-34`.
- **Cómo:** clasificar cada fase como "alimenta el informe" vs "artefacto auxiliar". Las auxiliares que no alimentan el informe pueden quedar tras un flag `--research-mode` (siguen disponibles), sin cambiar la salida por defecto del informe.
- **Esfuerzo:** 1 día. **Criterio:** el informe completo se sigue generando por defecto; nada de valor se pierde; el tiempo de `run` baja.

---

## TIER 3 — Rigor científico (tocar con cuidado, afecta decisiones de salud)

### 3.1 Completar la calibración PGS (ya en ROADMAP_NEXT.md #1) — ✅ **Hecho (2026-08-05)**
- **Qué:** solo 9 de 54 PGS se calculan (basados en posición); 45 necesitan mapeo rsID→pos.
- **Dónde:** `scripts/prs/`, `scripts/benchmarking/pgs_catalog_integration.py`, `scripts/utils/pgs_population_calibrate.py`, `scripts/utils/` (existe `pgs_rsid_mapper`, ver `tests/test_pgs_rsid_mapper.py`).
- **Esfuerzo:** 2–3 días. **Criterio:** 54/54 PGS con z-score + percentil.
- **Causa real (no era lo que decía el enunciado):** no era un problema de mapeo rsID→pos sin resolver — `pgs_catalog_integration.py` ya usa los archivos armonizados (`hmPOS_GRCh37`) de PGS Catalog, que traen `hm_chr`/`hm_pos` directo. El problema real eran **tres capas de artefactos desconectados**:
  1. El loop principal solo procesaba los top-3 resultados de una búsqueda por trait en cada corrida — los PGS ya descargados en corridas previas nunca se reintentaban si no volvían a aparecer en el top-3. Fix: `reprocess_downloaded_scores()` reprocesa TODO lo ya descargado en `pgs/`, re-descargando el archivo armonizado por ID directo si falta, y excluyendo los de >500K variantes (poco prácticos para una sola muestra). Resultado: 30→46 PGS puntuando.
  2. `pgs_population_calibrate.py` (el script que calcula z-score/percentil) **nunca estaba conectado a `prs.py`** — nadie lo invocaba, y su input (`prs/pgs_scores/pgs_results.csv`) llevaba congelado desde el 7 de junio. Conectado como paso `pgs_calibrate` después de `pgs_integration`, leyendo del `pgs/pgs_results.csv` real que ahora sí genera `pgs_catalog_integration.py` en cada corrida.
  3. El informe (`build_pgs_calibration_section`) lee `prs/pgs_scores/pgs_calibration_report.json` — un archivo JSON estructurado (distinto del CSV) que **tampoco generaba nada**, congelado desde el 16 de julio. `pgs_population_calibrate.py` ahora también escribe este JSON con nombre real de trait, n_snps y flag de confiabilidad.
  - Bug adicional encontrado corrigiendo esto: `--plink str(PLATFORM_DIR / "tools" / "plink")` en `prs.py` apuntaba a una ruta que no existe (`PLATFORM_DIR` es `prs_research_pipeline/`, el binario vive en la raíz del repo) — corregido a `PROJECT_ROOT`.
  - Optimización de performance: calibrar contra las 84M variantes del 1000G completo por cada uno de los ~50 PGS tomaba 1.5-2h+. Se extrae una sola vez el subconjunto de variantes realmente necesario (`--extract`) antes del loop — bajó a ~15 min.
- **Resultado verificado con el pipeline real (`prs.py run --full`) de punta a punta:** **52/57 PGS descargados puntuando y calibrados** (antes: 30, y esos 30 eran una foto congelada de junio/julio, no datos del run actual). 46 de los 52 son "reliable" (≤500K SNPs). El título de la sección del informe y las tarjetas de "Total Scores" ahora se calculan dinámicamente en vez de tener "30" hardcodeado.

### 3.2 Frecuencias alélicas reales gnomAD (ROADMAP_NEXT.md #2) — ⏸️ **Sin objeto hoy (verificado 2026-08-05)**
- **Qué:** reemplazar estimaciones HWE / MAF 0.25 por AFs reales para los ~30 SNPs que fallan calibración 1000G.
- **Dónde:** nuevo `reference/gnomad_af.json`, consumido en `scripts/prs/population_calibrate_v2.py`.
- **Esfuerzo:** 2 días. **Criterio:** menos fallbacks HWE; informe indica fuente de AF por SNP.
- **Estado:** antes de implementar, se verificó `reference/population_distributions/reference_distributions.json` (el archivo real que usa la calibración) — **0 de 59 rasgos usan fallback HWE hoy** (`"method": "Hardy-Weinberg estimate..."` no aparece en ninguna entrada; los 59 tienen distribuciones empíricas reales con `n_samples=503` de 1000G EUR). El "~30 SNPs" del enunciado es un número de antes de los fixes de posición (H1, commit `666f19b`) y strand-orientation (commit `754a972`) de esta sesión y la anterior — ese trabajo ya resolvió indirectamente el problema que 3.2 buscaba arreglar. No hay nada que implementar hasta que aparezcan nuevos fallbacks (p. ej. si se amplía el panel con SNPs no presentes en 1000G).
- **Nota aparte (no es fallback HWE, pero llamó la atención al verificar):** "Alcohol flush reaction" tiene `n_samples=503` pero **todos** los percentiles/media/mediana son idénticos (0.5), skewness/kurtosis `NaN` — una distribución degenerada, no una empírica real. Consistente con la limitación ya documentada de que `rs671` tiene cobertura pobre (no solo en el sample del usuario, sino aparentemente en gran parte del panel de referencia 1000G también). No investigado a fondo — candidato para revisar si se retoma el tema.

### 3.3 Auditoría/auto-fix de posiciones de SNP — **AHORA CRÍTICO** (ROADMAP_NEXT.md #8)
- **Qué:** verificar y corregir chr:pos GRCh37 de **todos** los rsID del panel contra dbSNP/Ensembl GRCh37.
- **Por qué:** ya se detectó al menos un error real (`rs671` ALDH2: CSV `chr12:112240098` vs correcto `chr12:112241766`) que hace que ese rasgo no puntúe. Posiciones erróneas = rasgos silenciosamente vacíos o mal calculados en un informe que se usa para decisiones. Esto es un **bug de corrección**, no una mejora opcional.
- **Dónde:** `scripts/setup/download_dbsnp.py`, `data/snp_database_annotated.csv`; feeds sub-tarea 1.2.A.
- **Cómo:** batch resolver rsID→pos (Ensembl GRCh37 REST o dbSNP local), diff contra el CSV, reportar y corregir. Añadir un test que falle si algún rsID diverge.
- **Esfuerzo:** 1 día. **Criterio:** 0 divergencias; rasgos afectados (empezando por alcohol flush) aparecen en el informe.

---

## TIER 4 — Testing, CI y calidad

### 4.1 Tests unitarios de la lógica científica (no solo integridad de JSON)
- **Qué:** hoy el "test suite" (`scripts/utils/test_suite.py`) valida que existan archivos y no haya tipos numpy. Los tests reales en `tests/` (~7 archivos) cubren utilidades, no el scoring.
- **Por qué:** un bug en el cálculo de z-score o en la orientación de alelos no lo detecta nada.
- **Dónde:** `tests/` (añadir `test_scoring.py`, `test_calibration.py`, `test_ancestry.py`), depende de 2.1.
- **Cómo:** con un genotipo de juguete y betas conocidas, verificar que PRS = Σ(β·dosis) exacto; verificar flip de strand con SNPs palindrómicos.
- **Esfuerzo:** 2–3 días. **Criterio:** cobertura de la ruta F–H con casos conocidos y de borde (alelo faltante, palindrómico, dosis 0/1/2).

### 4.2 Test end-to-end con VCF sintético en CI
- **Qué:** existe `test_samples_3eur.vcf.gz` y `test_samples`. Ejecutar el pipeline mínimo en CI.
- **Dónde:** `.github/workflows/test.yml`, `prs_research_pipeline/test_samples_3eur.vcf.gz`.
- **Cómo:** job que corre `run` en modo reducido sobre el VCF de prueba y compara salidas clave contra un golden file.
- **Esfuerzo:** 1–2 días. **Criterio:** CI ejecuta el pipeline y detecta regresiones de salida.

### 4.3 Endurecer el manejo de errores del orquestador
- **Qué:** `run_script` captura todo y a menudo continúa aunque falle (muchos `run_script(...)` sin `required=True` en `prs.py:504-560`). Un informe puede generarse con secciones vacías sin fallo claro.
- **Dónde:** `prs.py:147-205`, `prs.py:504-560`.
- **Cómo:** clasificar cada stage como crítico/opcional explícitamente; si un stage crítico falla, abortar con mensaje accionable. El informe final debe listar qué secciones faltan y por qué.
- **Esfuerzo:** 1 día.

---

## TIER 5 — Rendimiento y DX

### 5.1 Docker (ROADMAP_NEXT.md #5)
- **Qué:** `Dockerfile` con Python 3.12 + PLINK + bcftools + tabix + deps.
- **Por qué:** elimina el infierno de libs de sistema (cairo/pango de WeasyPrint, PLINK por OS).
- **Dónde:** nuevo `Dockerfile`, raíz. Referencia de tools en `tools/` y `README.md:70-79`.
- **Esfuerzo:** 1–2 días. **Criterio:** `docker run bluegen run --report nutrition --vcf /data/x.vcf.gz`.

### 5.2 Reducir overhead de subprocess / paralelizar Stage C
- **Qué:** Stage C (LD prune) tarda 45 min en frío (`prs_research_pipeline/README.md:195`). Arranque de Python ×40 por el patrón subprocess.
- **Dónde:** `scripts/stages/03_ld_ancestry_prune.sh`, `prs.py` (orquestación).
- **Cómo:** parte se resuelve con 2.1 (menos procesos). Para Stage C, cachear el resultado por hash del input (ya hay "cached 30s" — verificar que la cache funciona por defecto).
- **Esfuerzo:** 1–2 días.

### 5.3 Logging estructurado en vez de captura+recorte
- **Qué:** `run_script` vuelca stdout/stderr a `pipeline_debug.log` y muestra las últimas 8 líneas (`prs.py:184-199`). Existe `scripts/utils/logging_config.py` infrautilizado.
- **Dónde:** `prs.py:167-205`, `scripts/utils/logging_config.py`.
- **Cómo:** logging estructurado con niveles por stage; el informe enlaza al log del stage que falló.
- **Esfuerzo:** 1 día.

---

## Orden sugerido de ataque

1. **TIER 0 completo** (medio día) — quita fricción y riesgos inmediatos.
2. **3.3 + 1.2.A (auditoría de posiciones)** — **primero**, es un bug de corrección: hay SNPs que no puntúan (rs671 confirmado). Arregla el informe actual antes de ampliarlo.
3. **1.6** (refactor del report) — habilita crecer el informe sin romperlo.
4. **1.1 + 1.4 + 1.5** — resumen ejecutivo accionable + acciones con evidencia + contexto de confianza. Gran salto de valor sin quitar nada.
5. **1.2.B/C + 1.3** — más cobertura (posiciones ya verificadas en la tabla) + profundidad + integración de fuentes ya calculadas (el "completar con más información" que pediste).
6. **2.1** (librería importable) — desbloquea 4.1 (tests reales).
7. **4.1 / 4.2** — red de seguridad antes de tocar más ciencia.
8. **TIER 3** (PGS 45 restantes / gnomAD) — más datos que además alimentan 1.3.
9. **2.2 / 2.3 / TIER 5** — mantenibilidad y DX.

## Notas para Sonnet
- Antes de cada commit, seguir `memory/pre-commit-checklist.md`.
- No commitear nunca datos genómicos reales (VCF/BAM/FASTQ) ni informes con genotipos.
- El HTML (`reports/comprehensive_report_en.html`) es el entregable primario; el PDF es secundario y opcional.
- Cambios en la lógica de scoring/calibración requieren un test que los cubra (TIER 4).
