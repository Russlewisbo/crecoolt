# CRECOOLT — Novel vs. Older Regimens Analysis — Status

_Last updated: 2026-07-22 (redesigned to 2018–2024 cohort; 30-day mortality primary)_

## Project

Statistical re-analysis supporting the manuscript **"Evaluating the impact of novel
regimens in the management of CRE infections among liver transplant recipients"**
(Rinaldi et al., CRECOOLT project). Primary manuscript endpoint = 90-day mortality;
we were asked to add a competing-risks analysis of **CRE recurrence** and **resistance
emergence** accounting for death, and ended up auditing the whole IPTW mortality analysis.

## Data files (in project root)

| File | Rows | What it is |
|---|---|---|
| `CRECOOLT_onlyinfections.sav` | 405 | **THE analysis cohort** (1 row/patient with CRE infection). Use this. |
| `CRECOOLT_overall.sav` / `.csv` | 3170 | Older/partial REDCap export; only 190 infected. Do NOT use for the paper. |

Load: `new <- haven::read_sav("CRECOOLT_onlyinfections.sav")`

### Key variables in `new`
- Treatment: `new_drug` (1 = novel, 0 = old) → regimen (Old/Novel).
- Endpoints: `cre_recurrence`, `resistance` (0/1); `death`, `death_90d`, `death_30d`.
- Times: `CREinf_to_relapse`, `CREinf_to_resistance`, `CREinf_to_death`, `CREinf_deathCOX` (capped at 30d), `cre_infection_date`, `death_date`, `discharge_date`, `clearence_date`, `end_date`.
- PS covariates (9): `study_period`, `age`, `meld_score`, `sofa`, `mec_of_carbapenem_resi___1` (KPC), `infection_source___1` (BSI), `multisite_colonization`, `post_olt_compli___4` (PGNF), `post_olt_compli___2` (CRRT).
- Subgroup extra: `combination_treat`.

## Analysis decisions (CURRENT — 2018–2024 redesign, agreed with Matteo & Maddalena)
- **Cohort: 2018–2025 era ONLY (N=182)** — resolves positivity (novel only post-2018) and the differential ascertainment of 2010–2017. Full-cohort analyses abandoned.
- Time origin = date of CRE infection. Novel = ceftazidime-avibactam, mero/vaborbactam, imipenem/relebactam, cefiderocol.
- **PRIMARY endpoint: 30-day all-cause mortality** (IPTW-weighted Cox).
- **SECONDARY: 90-day CRE recurrence & resistance emergence** (IPTW Fine–Gray, death competing).
- **PS covariates (5): age, SOFA, BSI, multisite colonization, MELD**; stabilized ATE IPTW. Treatment-pathway variables (combination therapy, time-to-treatment) EXCLUDED as mediators.
- Death time = date-derived (`death_date − cre_infection_date`). 54 invalid discharge dates → censoring falls back to latest valid contact date.

## Key findings (current design)
1. **PRIMARY — 30-day mortality:** IPTW Cox **HR 0.49 (0.25–0.95), p=0.035** — ~51% reduction, significant. Crude old 32.9% vs novel 15.6%. (At 30 days "halving mortality" is accurate; the 90-day effect is weaker/borderline: 0.58, 0.34–0.99.)
2. **SECONDARY — recurrence:** IPTW Fine–Gray **sHR 0.98 (0.48–2.01), p=0.97** — comparable.
3. **SECONDARY — resistance:** IPTW Fine–Gray **sHR 1.79 (0.46–6.98), p=0.40** — comparable (11 events, wide CI).
4. **Subgroups (30-day):** effect consistent; all interaction p non-significant.
5. **Fine–Gray is invalid for all-cause mortality** (nothing competes with death) → Cox for mortality, Fine–Gray for recurrence/resistance only.

Rationale for the redesign (background): the full-cohort comparison was confounded by calendar era and by differential ascertainment in 2010–2017 (0 resistance events, ~0 recurrence, no clearance dates, outcomes coded "No"); positivity was violated because novel agents were used only from 2018.

## Figures produced (project root; PNG + PDF + SVG each), JAMA grey/blue, all 2018–2025 era
- `figure1_mortality_iptw.*` — 30-day cumulative mortality, HR 0.49, 95% CI bands, number-at-risk.
- `figure2_recurrence_resistance.*` — 90-day recurrence & resistance CIF, sHRs, death competing.
- `loveplot_era.*` — 5-covariate PS balance (age, SOFA, BSI, multisite, MELD).
- `figure_forest_30d.*` — 30-day mortality subgroup forest.
- STALE (prior 90-day / full-cohort design; delete before submission): `figure_3panel.*`, `figure_forest_90d.*`, `loveplot_iptw.*`, `competing_risks_model_table.*`.

## Website
- Quarto site publishes to **https://crecoolt.org** via `quarto publish netlify` (password: bolognacre).
- New page: `treatment_recurrence_resistance.qmd` ("Recurrence & Resistance"), published, leads with within-era analysis + ascertainment section + "Why death is a competing risk" callout.
- Nav switched to docked sidebar in `_quarto.yml`; project `render:` list scoped to 9 linked pages (bypasses broken `cid_letter_competing_risks.qmd` YAML `reference-doc: null`).

## Outstanding / TODO
- [ ] Manuscript edits (line-referenced list already provided): fix wrong Figure 1 legend (L245–246, it's the prophylaxis legend); remove Fine–Gray for mortality + redundant Table 2 column (L66,131–135,161–165); reframe recurrence/resistance from "survival bias" → ascertainment/era (L72–74,191–198); correct HR 0.37 → 0.58/clarify endpoint (L69–70,157); fix positivity/balance claim + replace Suppl Fig 1 with era version (L128–130,152–156); reconcile 30-vs-90-day endpoint (L116,131–135); complete truncated software sentence (L138).
- [ ] Draft note to authors is written (5 points) — send/adapt.
- [ ] Confirm with authors: retrospective follow-up culture practice (ascertainment); 8 `resistance=0` records that carry a resistance date; source Stata script to reconcile 0.37.
- [ ] Optional: recover data-collection variable to restore forest study-period stratum; add figures to the qmd page; update combined legend already set to 0.58.

## Reproducible script
- **`analysis.R`** — single self-contained script for the 2018–2024 / 30-day design.
  `source("analysis.R")` reloads the data, rebuilds the era cohort, PS/IPTW, and all models,
  prints the key estimates, and regenerates the four figures (30-day mortality Figure 1,
  recurrence/resistance, era Love plot, 30-day forest). Fastest way to resume after closing the IDE.

## Website pages
- `treatment_recurrence_resistance.qmd` → "Recurrence & Resistance" (competing-risks + IPTW + ascertainment).
- `mortality_iptw.qmd` → "Mortality (IPTW)" (finalized mortality figures, forest, Love plots, reproduction notes).
- Both are registered in `_quarto.yml` (sidebar + `render:` list) and published via `quarto publish netlify`.

## To resume in R
Key reproducible objects (regenerate by re-running the executeCode cells, or):
`new` (data), `d2`/`fdat` (mortality analysis frame + 9-cov PS weights `w`), `cr`/`era2`/`e_rec`/`e_res`
(competing-risks frames), `make_cr()` (builds time/status), `era_fg()` (weighted Fine–Gray),
`forest`/`pints` (subgroup results). Palette: `jama2 <- c(Old="#8C8C8C", Novel="#00468B")`.
