# Replication Plan: Braithwaite & Lemke (2011), "Unpacking Escalation"

> **Status (2026-10-04): replication complete. See [REPORT.md](REPORT.md) for results.** This is the initial plan, written before the full text was available. Items marked 🔶 were guesses and have since been corrected against the article; Section 2 of REPORT.md has the specification actually used. For example, the replication uses MID v3.02 (not v3.10), the sample period is 1919–1995 (not 1816–2001), and the article's six escalation measures differ from the candidate list Y1–Y7 below.

- **Article:** Alex Braithwaite & Douglas Lemke. 2011. "Unpacking Escalation." *Conflict Management and Peace Science* 28(2): 111–123. DOI: [10.1177/0738894210396631](https://doi.org/10.1177/0738894210396631)
- **Keywords** (journal page): concept measurement, conflict escalation, MIDs, selection effects, territory
- **Plan date:** 2026-10-04

> ⚠️ **Sources of information at the time of writing**
> - The article is closed-access on SAGE. Unpaywall and Semantic Scholar list no open-access copy, and the journal site returned HTTP 403, so **the full text was not available**.
> - Items marked ✅ were verified from the abstract, keywords, Crossref metadata, or downloaded data.
> - Items marked 🔶 were inferred from the literature the article builds on (Reed 2000; Senese 2005; Senese & Vasquez 2008, and other MID escalation studies). **They had to be checked one by one against Sections 3–4 and the tables once the PDF was obtained.**
> - No author replication package was found on Harvard Dataverse, the journal site, or through web search.

---

## 1. Research question ✅

Studies measure "conflict escalation" in very different ways, and theory does not specify how it should be measured. The authors run the same set of explanatory variables against **several measures of escalation** and compare whether the conclusions agree.

**Core findings** (from the abstract):
1. Most relationships are not robust when the measure of escalation changes ("few robust relationships").
2. The positive effect of **territorial issues** on escalation holds under every measure.

The target of the replication is therefore **not a single coefficient** but a "measure × explanatory variable" results matrix: the territory row should be significantly positive in every column, while the sign or significance of the other variables varies across measures.

---

## 2. Results to replicate

| ID | Target | Description |
|---|---|---|
| R1 | Descriptive statistics | Distribution of each escalation measure (number of MIDs and share escalated) 🔶 |
| R2 | Main regression table(s) | One column per escalation measure, one row per explanatory variable; coefficients, standard errors, significance 🔶 |
| R3 | Selection equation (if any) | First-stage estimates for MID onset, plus ρ (correlation between the two equations' errors) 🔶 |
| R4 | Robustness summary | Which variables are consistently signed and significant across all measures (the article's claim: only territory) ✅ |

**Success criteria** (tiered, since no original replication package exists):
- **Strong replication:** coefficients within ±10% of the published values, with the same or similar N.
- **Substantive replication** (realistic target): every variable has the same **sign and significance** as published under every measure, especially territory being significantly positive in every column.
- **Conclusion-level replication** (minimum): the qualitative conclusion in R4 holds.

---

## 3. Data

### 3.1 Downloaded (in `data/raw/`) ✅

| Dataset | File | Use | Notes |
|---|---|---|---|
| COW MID **v3.10** (1816–2001) | `MID_v3.10.zip` → `mid310/` | Main source of the dependent variables and the territory variable | Current version around 2011; removed from the COW site and recovered from the Wayback Machine. Includes MIDA/MIDB/MIDI/MIDIP |
| Dyadic MID 4.01 | `dyadic_mid_4.01.zip` | Backup for building dyad-level MIDs | Closest match if the article used Maoz's DYMID |
| COW NMC **v4.0** | `NMC_v4_0.csv` | CINC → capability ratio / parity | Matches the 2011 era |
| COW Direct Contiguity 3.2 | `DirectContiguity320.zip` | Contiguity | |
| COW Formal Alliances 4.1 | `Alliance_v4.1_csv.zip` | Alliance ties | v3.03 was current in 2011; small differences possible |
| COW Major Powers / State System 2024 | `MajorPowers2024.zip`, `States2024.zip` | Major-power status; dyad-year frame | Essentially identical to older versions for 1816–2001 |
| Polity IV/5 (p5v2018) | `p5v2018.xls` | Joint democracy | The 2011-era version was Polity IV 2008/2009; early years may have minor revisions |

Quick check of MID v3.10: 2,332 disputes, 1816–2001. `hostlev`: 2=103, 3=569, 4=1,553, 5=107 (war). About 400 disputes have `fatality>0` (272 are missing = −9); 1,147 have `recip=1`.

### 3.2 Possibly needed 🔶

- **Alternative coding of territorial issues:** MID's own `revtype1/revtype2 == 1` (territory) is most likely what the article used. If it used Hensel's ICOW territorial claims instead, those would come from paulhensel.org (coverage is limited to the Americas, Europe, and the Middle East).
- **COW War data v4.0:** needed if "war" is defined by the COW war list rather than `hostlev==5`.
- **EUGene:** MID studies from around 2011 commonly generated dyad-year data (peace years, etc.) with EUGene. EUGene is now hard to run, so equivalent variables would be built in R/Python.
- **ATOP** (if the article used ATOP rather than COW alliance data).

---

## 4. Modeling approach 🔶 (to be corrected against the full text)

### 4.1 Unit of analysis and sample
- **Most likely design:** dyad-years (1816–2001) with selection at the MID-onset stage; alternatively, MID-level data with escalation as the dependent variable.
- The keyword "selection effects" suggests the article **most likely uses a two-stage selection model**:
  1. Selection equation: whether a MID occurs in the dyad-year;
  2. Outcome equation: among MIDs that occur, whether the dispute escalates (under measure k).

### 4.2 Candidate escalation measures (one model per measure)

Measures common in the literature are listed below; which ones the article uses had to be confirmed from the PDF:

| Code | Definition | Source variable |
|---|---|---|
| Y1 War | `hostlev == 5` or ≥1,000 battle deaths | MIDA `hostlev`; or COW War |
| Y2 Fatalities | `fatality ≥ 1` | MIDA `fatality` |
| Y3 Use of force | `hostlev ≥ 4` | MIDA `hostlev` |
| Y4 Reciprocation | `recip == 1` | MIDA `recip` |
| Y5 Ordinal hostility | `hostlev` (2–5, ordered) | MIDA `hostlev` |
| Y6 Highest action | `hiact` (ordered or dichotomized) | MIDA/MIDB `hiact` |
| Y7 Escalation magnitude | Highest minus initial hostility level | Computed from MIDI incident data |

### 4.3 Explanatory variables (standard MID-escalation setup)
- **Territory:** `revtype1/2 == 1` for either side (key variable)
- **Joint democracy:** both Polity scores ≥ 6 (or the weak-link minimum)
- **Capability parity / ratio:** weaker/stronger CINC, or the log ratio
- **Contiguity:** land border, or separated by ≤ 400 miles of water
- **Alliance:** defense pact / any alliance
- **Major-power dyad:** at least one major power
- **Peace years + splines or cubic polynomial** (selection equation; Beck–Katz–Tucker)
- Possibly also: rivalry, MID history, Cold War dummy

### 4.4 Estimation
1. **Main models:** one Heckman probit per Y_k (Stata `heckprob`, R `sampleSelection::selection`); ordered probit or dichotomized probit for ordinal outcomes (Y5, Y6).
2. **Comparison:** ordinary probit/logit on the MID sample only, ignoring selection.
3. Standard errors clustered by dyad.
4. The selection equation needs an **exclusion restriction**, a variable in the selection equation but not in the outcome equation. Peace years or distance are common choices; the article's choice takes precedence.

**Software:** Python (pandas for data) + R (`sampleSelection`, plus `oglmx`/`MASS::polr` if needed). For digit-level comparison, Stata `heckprob` could be used to cross-check.

---

## 5. Computing and comparing results

1. For each Y_k, take β, standard errors, p-values, and N from the outcome equation, and ρ with its test from the selection model.
2. Assemble a replication table in the same layout as the published one (columns = measures, rows = variables), with the published value next to each cell.
3. **Robustness matrix:** mark each variable × measure as `+*` / `−*` / `n.s.`. Core criterion: the territory row is `+*` under every measure.
4. Diagnosing discrepancies: if N does not match, check in order: sample period → how multilateral disputes are split into dyads → missing-value handling (`fatality = −9`) → dataset version differences.
5. Optional extension: rerun with MID 5.0 (the Gibler–Miller–Little revision) in place of MID v3.10 to see whether the territory result still holds.

---

## 6. Workflow

| Step | Task | Output |
|---|---|---|
| 0 | **Obtain the PDF** (UW library SAGE subscription), place it in `paper/`, and revise Sections 3–4 of this plan | Final plan |
| 0b | (Optional) Email the authors for do-files/data (Braithwaite is now at the University of Arizona) | Original replication package |
| 1 | Clean MID 3.10; build MID-level and dyad-level data; compute Y1–Y7 and territory | `data/clean/mid_dyad.csv` |
| 2 | Build the 1816–2001 dyad-year frame; merge NMC, Polity, contiguity, alliances, major-power status, peace years | `data/clean/dyad_year.csv` |
| 3 | Reproduce descriptive statistics (R1) | `output/table1.*` |
| 4 | Estimate Heckman probits and ordinary probits (R2, R3) | `output/table2.*` |
| 5 | Robustness matrix and comparison with the article (R4) | `output/robustness_matrix.*` |
| 6 | Write the replication report: what matched, what did not, and likely reasons | `REPORT.md` |

---

## 7. Feasibility assessment

**Replicable, but most likely at the "substantive" level rather than digit for digit.**

- ✅ **Data: available.** Everything used is public COW/Polity data, including the 2011-era MID v3.10 (recovered from the Wayback Machine and downloaded).
- ✅ **Methods: implementable.** Probit, ordered probit, and Heckman selection models are standard and run directly in R or Stata.
- ⚠️ **Main obstacle: no access to the full text.** Without the PDF, the exact specification (which escalation measures, controls, sample, and exclusion restrictions) could not be verified. **The PDF had to be obtained through the university library.**
- ⚠️ **No official replication package**, and the article may have generated its data with EUGene, so N and coefficients will likely differ slightly from the published values. The core conclusion (territory is robust across all measures) is likely to replicate.
