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
- **This repo:** R package `sweetsoursong` (`R/`, `src/`) and `_scripts/`, `_mathematica/`, `_data/`, `_figures/`. These implement the **old** model (ODE single/multi-plant and SDE seasonal simulations) from the first submission; they produce none of the current figures.
- **Numerical checks:** `claude-checks/` (scripts with matching `*_output.*` files; summary in `note_for_chris.md`).
- **Drafts and background:** `~/Box/_sweetandsour/_drafts/`, `zzz-background/`, `zzz-outdated/`.

## How this project works

- Model hierarchy in the current manuscript: the full per-plant CTMC (Y, B, P; Table 1); the reduced model with no empty flowers plus a vacancy (CTMC) closure for the filling probabilities; the closed metacommunity, where regional pools satisfy P_YR = E[YP] at fixed P_R.
- Figure provenance in `exploration_unified_10x.nb`: Fig 2 = section 1.4 (deterministic one-plant model). Figs 3–4 = 2.2, Fig 5 = 2.3, Fig 6 = 2.4 (reduced model, vacancy closure). Section 3 (full model) is used only for the SI check.
- Parameters come from `SetParameters` in the notebook: N = 50, Pmax = 12, c = 500, d = 0.1, m = 0.01, m_B = 0.05, e_Y = 1, e_B = 0.5, c_B∅ = 5, ε = 1e-5, giving P_crit = 1.
- Software: Wolfram 15.0.1 with EcoEvo 1.7.2; R 4.4.3 with `renv`.
- Manuscript edits are tracked with the LaTeX `changes` package (`\added`, `\deleted`, `\replaced`, `\comment`). New display equations and tables use `addedblock` / `\addedcolor`. Use `defaultcolor=red`. Accept all edits with `\usepackage[final]{changes}`.
- Old equations and figures are hidden in `\iffalse` behind visible markers, with label shims that print "old N". Leave them until all changes are accepted.
- Overleaf word count runs TeXcount, which ignores `\iffalse`. Wrap hidden code in `%TC:ignore` / `%TC:endignore`. The `%TC:macro` lines at the top of `__ms.tex` skip `\deleted`, the old argument of `\replaced`, and `\comment`. Check locally with `texcount -inc -total __ms.tex`.
- Overleaf sync order: the user pushes from Overleaf → `git fetch` and fast-forward locally → commit → check `origin/main` is an ancestor of `HEAD` → push (with the user's say-so) → the user pulls in Overleaf. The user often edits on Overleaf, so do not edit `sweetsoursong-ms` locally without asking.

## Decisions that are settled

- **2026-10-02 —** Never rewrite history or force-push `sweetsoursong-ms`, and keep figure files tracked. The Overleaf project cannot be unlinked, Overleaf's sync ignores `.gitignore`, and a rewrite on this date broke the sync until it was undone.
- **2026-10-02 —** Report only whether each invasion criterion exceeds 1, not its magnitude. Thresholds are stable across ε = 1e-5 to 1e-9 and Pmax = 12 to 18; magnitudes change by up to 100× with ε (SI Table S1).
- **2026-10-02 —** The abstract says coexistence "across a wide range of regional pollinator abundances" requires the pollinator feedback, not that coexistence requires it outright. With no preference (m_B = m), coexistence still occurs for 1.46 < P_R < 1.57 (`claude-checks/ms_numbers_output.m`).
- **2026-10-02 —** Main-text results use the reduced model with the vacancy closure. The full model appears only as an SI accuracy check: its thresholds differ by 0.001 for yeast invasion and 0.03 for bacteria invasion.

## Working notes

- Current status: `PROJECT_INDEX.md`
- Live work: `TODO.md`
