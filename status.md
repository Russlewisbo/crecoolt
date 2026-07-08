# CRECOOLT — Novel vs. Older Regimens Analysis — Status

_Last updated: 2026-07-08_

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

## Analysis decisions (agreed)
- Time origin = date of CRE infection.
- Novel = ceftazidime-avibactam, mero/vaborbactam, imipenem/relebactam, cefiderocol.
- Competing risk = death, for recurrence/resistance only (NOT mortality — see findings).
- IPTW: full-cohort for mortality; recurrence/resistance restricted to **2018–2025 era** (positivity).
- 54 discharge dates precede infection date → treated invalid; censoring falls back to latest valid contact date.

## Key findings
1. **Mortality benefit is real** (~40–50% reduction). Reproducible **90-day IPTW HR = 0.58 (0.39–0.86, p≈0.007)**.
2. **Reported HR 0.37 is NOT reproducible** — it's a 30-day estimate; our specs gave 0.49–0.58 (30-day ATT 0.52, p=0.068, non-significant).
3. **Fine–Gray is invalid for all-cause mortality** — nothing competes with death; Cox HR and sHR both 0.50 (identical) confirms it adds nothing. Reserve Fine–Gray for recurrence/resistance.
4. **Recurrence/resistance: crude excess with novel is an era/ascertainment artifact, NOT survival bias.**
   - Full-cohort competing-risks: recurrence sHR 4.58, resistance 7.12 (excess persists → not survival bias).
   - Within-era (2018–2025) IPTW: recurrence sHR **1.03 (0.52–2.05)**, resistance **1.45 (0.38–5.55)** → null.
   - 2010–2017 ascertainment: 0 resistance events, 1.4% recurrence (vs 19.7%), zero clearance dates, outcomes coded "No" not NA.
5. **Supplementary Figure 1 balance requires ATT, not the stated ATE.** Study period cannot be balanced full-cohort (perfect separation, novel only post-2018). ATT drives pre-2018 old weights to ~0 = effectively era restriction. Balance IS achievable within the 2018–2025 era (max |SMD| 0.032).

## Figures produced (project root; PNG + PDF + SVG each), JAMA grey/blue palette
- `figure1_mortality_iptw.*` — cumulative mortality, HR 0.58, 95% CI bands, number-at-risk.
- `figure_3panel.*` — A mortality / B recurrence / C resistance (sHRs, death competing, 0–90d).
- `loveplot_iptw.*` — full cohort (study period unbalanceable).
- `loveplot_era.*` — 2018–2025 era (all balanced, max |SMD| 0.032).
- `figure_forest_90d.*` — 90-day IPTW subgroup forest, overall HR 0.58; subgroups KPC/BSI/SOFA/therapy-type, all interactions n.s. (study-period stratum omitted — data-collection variable absent).
- `competing_risks_model_table.*` — all model variants table.

## Website
- Quarto site publishes to **https://crecoolt.org** via `quarto publish netlify` (password: bolognacre).
- New page: `treatment_recurrence_resistance.qmd` ("Recurrence & Resistance"), published, leads with within-era analysis + ascertainment section + "Why death is a competing risk" callout.
- Nav switched to docked sidebar in `_quarto.yml`; project `render:` list scoped to 9 linked pages (bypasses broken `cid_letter_competing_risks.qmd` YAML `reference-doc: null`).

## Outstanding / TODO
- [ ] Manuscript edits (line-referenced list already provided): fix wrong Figure 1 legend (L245–246, it's the prophylaxis legend); remove Fine–Gray for mortality + redundant Table 2 column (L66,131–135,161–165); reframe recurrence/resistance from "survival bias" → ascertainment/era (L72–74,191–198); correct HR 0.37 → 0.58/clarify endpoint (L69–70,157); fix positivity/balance claim + replace Suppl Fig 1 with era version (L128–130,152–156); reconcile 30-vs-90-day endpoint (L116,131–135); complete truncated software sentence (L138).
- [ ] Draft note to authors is written (5 points) — send/adapt.
- [ ] Confirm with authors: retrospective follow-up culture practice (ascertainment); 8 `resistance=0` records that carry a resistance date; source Stata script to reconcile 0.37.
- [ ] Optional: recover data-collection variable to restore forest study-period stratum; add figures to the qmd page; update combined legend already set to 0.58.

## To resume in R
Key reproducible objects (regenerate by re-running the executeCode cells, or):
`new` (data), `d2`/`fdat` (mortality analysis frame + 9-cov PS weights `w`), `cr`/`era2`/`e_rec`/`e_res`
(competing-risks frames), `make_cr()` (builds time/status), `era_fg()` (weighted Fine–Gray),
`forest`/`pints` (subgroup results). Palette: `jama2 <- c(Old="#8C8C8C", Novel="#00468B")`.
