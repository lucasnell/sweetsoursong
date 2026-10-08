# sweetsoursong — project context

## Project

- **Code:** sweetsoursong
- **Full name:** Regional species coexistence despite local priority effects: the overlooked role of dispersal–community feedback
- **My role:** lead and corresponding author
- **Question:** Can feedback between pollinators and nectar microbes (pollinators avoid bacteria-soured nectar) let competing yeast and bacteria coexist regionally even though priority effects lead to exclusion within flowers and plants?

## People

| Who | Role | What they own |
|---|---|---|
| Lucas Nell | Lead author | Manuscript, R package, figures |
| Christopher Klausmeier | Co-author | New-model formulation and Mathematica code (`definitions_unified.wl`, `exploration_unified_10x.nb`, vacancy-CTMC closure note) |
| Tadashi Fukami | Co-author | Nectar-microbe study system |

## Where things live

- **Manuscript:** `~/GitHub/Stanford/sweetsoursong-ms`, symlinked here as `.ms` (git-ignored). Remote `github.com/lucasnell/sweetsoursong-ms`, linked to an Overleaf project by GitHub sync. Target: *PNAS* Research Report.
- **Manuscript figures:** `~/Box/_sweetandsour/_figures/` (`.ai` sources and `.pdf` exports), copied into `sweetsoursong-ms/figures/` (Box `nectar-model-diagram.pdf` → `fig1-model-diagram.pdf`).
- **New-model code (Chris):** originals in `~/Box/_sweetandsour/ChrisK-nb/`; working copies in `claude-checks/chris-files/`. The 12 MB notebook is git-ignored.
- **This repo:** on `main`, R package `sweetsoursong` 2.0.0, a pure-R port of Chris's model (`R/`, `tests/`), and `_scripts/`, which make the computed panels of Figs 2–6 and SI Table S2 (outputs in `_data/`, `_figures/`, git-ignored). At tag `v1.0.0`: the **old** package (ODE/SDE model of the first submission), which produces none of the current figures.
- **Numerical checks:** `claude-checks/` (scripts with matching `*_output.*` files; summary in `note_for_chris.md`). R-port checks in `claude-checks/r-port/`: Mathematica reference script, comparisons of R with Mathematica and with the manuscript.
- **Drafts and background:** `~/Box/_sweetandsour/_drafts/`, `zzz-background/`, `zzz-outdated/`.

## How this project works

- Model hierarchy in the current manuscript: the K-plant metacommunity in which pollinators fly directly between plants (SI Table S1, `tab:rates-multi`), whose one-plant view with a regional pool is the full per-plant CTMC (Y, B, P; transitions in Fig 1E, parameters in Fig 1D; summation checked in `claude-checks/multi_plant_table_check.R`); the reduced model with no empty flowers plus a vacancy (CTMC) closure for the filling probabilities; the closed metacommunity, where regional pools satisfy P_YR = E[YP] at fixed P_R.
- Manuscript figure panels are the Mathematica exports (decided 2026-10-05); `_scripts/` reproduces their computed content in R. Figure provenance in `exploration_unified_10x.nb`: Fig 2 = section 1.4 (deterministic one-plant model). Figs 3–4 = 2.2, Fig 5 = 2.3, Fig 6 = 2.4 (reduced model, vacancy closure). Section 3 (full model) is used only for the SI check.
- Parameters come from `SetParameters` in the notebook: N = 50, Pmax = 12, c = 500, d = 0.1, m = 0.01, m_B = 0.05, e_Y = 1, e_B = 0.5, c_B∅ = 5, ε = 1e-5, giving P_crit = 1.
- Software: Wolfram 15.0.1 with EcoEvo 1.7.2; R 4.6.1 with `renv`. `wolframscript` needs `WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel`; run one kernel-heavy script at a time (two in parallel crashed).
- R port: tests compare with Mathematica fixtures (`tests/testthat/fixtures/ref_chop.json`, `ref_nochop.json`) made by `claude-checks/r-port/r_port_reference.wl`. `options(sweetsoursong.chop = TRUE)` and `inv_b(method = "chris")` reproduce Chris's code exactly.
- Uncolonized flowers are written `\varnothing` (changed from `\emptyset` on 2026-10-08, `sweetsoursong-ms` `a8ebbfc`).
- The tracked changes of the new-model rewrite were accepted on 2026-10-05 (`sweetsoursong-ms` `5620296`); new edits go in without markup unless the user asks. The `changes` package (`\added`, `\deleted`, `\replaced`, `\comment`, `addedblock`, `defaultcolor=red`) is still loaded in `__ms.tex` if a tracked round is needed again.
- Overleaf word count runs TeXcount, which ignores `\iffalse`. Wrap any hidden code in `%TC:ignore` / `%TC:endignore`. The `%TC:macro` lines at the top of `__ms.tex` skip `\deleted`, the old argument of `\replaced`, and `\comment`. Check locally with `texcount -inc -total __ms.tex`.
- Overleaf sync order: the user pushes from Overleaf → `git fetch` and fast-forward locally → commit → check `origin/main` is an ancestor of `HEAD` → push (with the user's say-so) → the user pulls in Overleaf. The user often edits on Overleaf, so do not edit `sweetsoursong-ms` locally without asking.

## Decisions that are settled

- **2026-10-02 —** Never rewrite history or force-push `sweetsoursong-ms`, and keep figure files tracked. The Overleaf project cannot be unlinked, Overleaf's sync ignores `.gitignore`, and a rewrite on this date broke the sync until it was undone.
- **2026-10-02 —** Report only whether each invasion criterion exceeds 1, not its magnitude. Thresholds are stable across ε = 1e-5 to 1e-9 and Pmax = 12 to 18; magnitudes change by up to 100× with ε (SI Table S2, `tab:numerics`).
- **2026-10-02 —** The abstract says coexistence "across a wide range of regional pollinator abundances" requires the pollinator feedback, not that coexistence requires it outright. With no preference (m_B = m), coexistence still occurs for 1.46 < P_R < 1.57 (`claude-checks/ms_numbers_output.m`).
- **2026-10-05 —** Do not write a release-specific DOI into this repo; the Zenodo badge in `README.md` always points to the latest release. A release's DOI exists only after it is published, so a written one is always one release behind. The manuscript cites the concept DOI, 10.5281/zenodo.15113987, which always resolves to the latest release (2026-10-08).
- **2026-10-05 —** The new model is implemented in pure R (no C++), replacing the old package in this repo; the old one stays at tag `v1.0.0`. Defaults differ from Chris's code in two ways: no `Chop`, and R_B is the invader's own output E[BP]/ε. Neither changes Figs 2–6 or any reported threshold (`claude-checks/r-port/compare_chop_output.txt`, `compare_invb_output.txt`).
- **2026-10-02 —** Main-text results use the reduced model with the vacancy closure. The full model appears only as an SI accuracy check: its thresholds differ by 0.001 for yeast invasion and 0.03 for bacteria invasion.

## Working notes

- Current status: `PROJECT_INDEX.md`
- Live work: `TODO.md`
