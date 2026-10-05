# sweetsoursong — index

**Status as of 2026-10-05:** The manuscript has been rewritten for the new pollinator-pool metacommunity model (Figs 1–6), with every edit tracked; the user is editing it on Overleaf. The model is now implemented in R (package 2.0.0 on `main`) and reproduces Chris's Mathematica results and the manuscript numbers. Next milestones: co-author read-through; the user's GitHub and Zenodo release of the R code, which unblocks the data and code statement.

## The question

Priority effects make each flower, and with pollinator feedback each plant, dominated by either yeast or bacteria. Does the same feedback, in which pollinators leave bacteria-dominated flowers more readily, let both microbes coexist across a landscape of plants that share a regional pollinator pool? If so, dispersal–community feedback is a coexistence mechanism that metacommunity theory has largely overlooked.

## Approach

A continuous-time Markov chain for one plant with N flowers (uncolonized, yeast-dominated, or bacteria-dominated) and P pollinators. The plant is coupled to a regional pool through pollinator immigration and emigration. Because pollinators move much faster than flowers turn over, a reduced model assumes no empty flowers and computes filling probabilities with a vacancy Markov chain (closure). A closed metacommunity sets the regional pools equal to the plants' expected pollinator output at fixed regional pollinator abundance P_R. Coexistence means mutual invasibility, judged by invasion criteria at a small invader pool ε.

## Current state

| Workstream | State | Next |
|---|---|---|
| Manuscript | Rewritten with tracked changes; Significance Statement revised; Lerch et al. citations and founder-control sentence added (2026-10-05); word count works on Overleaf | Co-author read-through; accept changes |
| New-model code | R port (package 2.0.0) reviewed, approved by Chris, merged into `main` and pushed: 346 tests against Mathematica pass; `R CMD check` 0 errors, 0 warnings; scripts reproduce Figs 2–6 (computed panels) and Table S1 | User makes the GitHub release `v2.0.0` and Zenodo archive |
| Numerical checks | Done (`claude-checks/`), including founder control (none for m_B ≥ m) | Decide whether to send `note_for_chris.md`; whether to scan founder control at other e_B, c_B∅ |
| Old R package | Replaced on `new-model`; preserved at tag `v1.0.0` (Zenodo 10.5281/zenodo.15113988) | None |

## Key links

- Manuscript repo: `~/GitHub/Stanford/sweetsoursong-ms` (Overleaf-synced), `github.com/lucasnell/sweetsoursong-ms`
- Figures: `~/Box/_sweetandsour/_figures/`
- Chris's model files: `~/Box/_sweetandsour/ChrisK-nb/`
- Checks: `claude-checks/note_for_chris.md`
- Previous code archive (old model): Zenodo DOI 10.5281/zenodo.15113988

## Decision log

| Date | Decision | Why |
|---|---|---|
| 2026-10-02 | Main text uses the reduced model with the vacancy closure; the full model is an SI check only | Its thresholds are within 0.001 (yeast) and 0.03 (bacteria) of the full model's |
| 2026-10-02 | Report only the sign of R − 1, not invasion-criterion magnitudes | Magnitudes are not converged in ε; thresholds are |
| 2026-10-02 | Abstract hedged to "across a wide range of regional pollinator abundances" | Narrow coexistence window (1.46 < P_R < 1.57) with no pollinator preference |
| 2026-10-02 | Track manuscript edits with the `changes` package, in red | User can see every change on Overleaf |
| 2026-10-02 | Word count removed from the title page | Stale after the rewrite |
| 2026-10-02 | Never rewrite `sweetsoursong-ms` history; keep figures tracked | The Overleaf project can't be unlinked, and its sync ignores `.gitignore` |
| 2026-10-05 | Manuscript figure panels stay as the Mathematica exports for now | User's choice; the R scripts reproduce the computed panels |
| 2026-10-05 | Implement the new model in pure R, replacing the old package in this repo (old one at tag `v1.0.0`) | Open-source and reproducible without a Wolfram licence; the numerics are small enough for R |
| 2026-10-05 | R port defaults: no `Chop`; R_B = E[BP]/ε; truncated Poisson weights in log space; with P_BR exactly 0, all mass at y = N | Chop and Chris's InvB lose precision; Figs 2–6 and thresholds are unchanged (`claude-checks/r-port/compare_*_output.txt`) |
| 2026-10-04 | TeXcount skips the hidden Mathematica block and the text in `\deleted`, old `\replaced`, and `\comment` (`%TC:` lines in `__ms.tex`, `02-methods.tex`) | Overleaf's word count errored on `\[EmptySet]` and counted deleted text |
