# Handoff — 2026-10-30

## Session topic

Reproducibility review of this repo and the manuscript, then the fixes. Two Sonnet reviewers (package and scripts; claim provenance) found no wrong numbers; the gaps were provenance, the archive lagging the manuscript, and private files in the public repo. Their reports are in git history only (`23da609`, removed in `d1444de`).

## Key decisions

- Manuscript keeps the Mathematica figure panels; Methods now says so and that the R package reproduces every computed result (random draws in Figs 3A,B and 4 match in distribution only).
- Both `(SI?)` placeholders resolved by a new full-vs-reduced comparison of distributions (option 19), not by narrowing the claim to thresholds. Cited as "SI Appendix, Numerical accuracy"; no new SI table.
- `note_for_chris.md`, `cover-letter/`, `ms-edit/` untracked (still in history and in the v2.0.0 Zenodo record; no history rewrite); these and `reviewers/`, `param-justification/` are in `.gitignore`. Chris agreed to his files (`chris-files/`) being public.
- Release `v2.1.0`; manuscript keeps citing the concept DOI.

## Done

- `sweetsoursong-ms` `c251a76` (pushed): placeholders; SI sentence on reduced vs full distributions; SI Pmax sentence (tail ≤ 4.2e-5, pools ≤ 0.002% reduced, 0.02% full); Table S2 2.89581 → 2.895811; Methods software sentence. Word count 5190.
- `sweetsoursong` `611889e`: `claude-checks/founder_control_tip.R` (tip P_R = 1.5276369729, m_B = 0.0090761076, matches Mathematica to 10 digits), `pmax_check.R` (reduced thresholds change ~1e-6 across Pmax 12/18/24), `full_vs_reduced.R`, expanded `r-port/verify_manuscript_numbers.R` (all match).
- `a2c7641`: README (what scripts reproduce, runtimes, `testthat::test_local()` with `NOT_CRAN=true`, fixture regeneration), `claude-checks/README.md` (statement → script → output), Wolfram scripts use paths relative to themselves, `cached()` keyed on version/params/Chop, one core on Windows, `inst/CITATION`, version 2.1.0, `full_meta_equilibrium()` unnames a named start.
- Verification before tagging: all five `_scripts/` re-ran under 2.1.0; 346 tests pass; `R CMD check --as-cran` 0/0/1 (New submission).
- Tag `v2.1.0` pushed; user published the release. Concept DOI 10.5281/zenodo.15113987 resolves to the v2.1.0 record, 10.5281/zenodo.23817877 (checked 2026-10-30).

## Full vs reduced (claude-checks/full_vs_reduced_output.txt)

Over 0.8 ≤ P_R ≤ 2.8, full minus reduced: 0.27–0.61 uncolonized flowers per plant; yeast-dominated proportion +0.001 to +0.014; mean P per plant ≤ 3e-5; mean P on yeast plants −0.03 to −0.04; mean Y within plant types ≤ 0.3; P_BR −0.47 (moved to the P_0R pool). Wasserstein-1 ≤ 0.75 flowers. Total variation of (Y, P) reaches 0.19 at P_R = 2.8 only because Y = 50 in the reduced model is often Y = 49 plus one empty flower in the full model.

## Open follow-ups

- [ ] User: pull `c251a76` into Overleaf (if not done)
- [ ] Fit the main text to PNAS limits: 4 graphical elements (now 6 figures), 4,000 words (now 5190)
- [ ] Not verified: the no-renv install route at `@v2.1.0` (tested only for v2.0.0); heavy Wolfram scripts not re-run after the path change (mechanism checked by re-running `extract_notebook_inputs.wl` from `/`, identical output)
- [ ] Deferred: Linux system dependencies for the PPM route in README
- [ ] From 2026-10-08, still open: cover-letter questions (signature affiliation, "disease ecology", mention of the *Ecology Letters* submission)

## Context for the next session

- `chop_threshold_check.txt` has no script (interactive Mathematica session); `claude-checks/README.md` says so.
- Local `_data/` still holds v1-era files, including a 2.2 GB `.rds` (git-ignored; clutter only).
- `git commit` after `git rm --cached` commits the whole index, not just the paths added; stage deliberately when splitting commits.
- Earlier handoffs are in git history (`git log -p handoff.md`).
