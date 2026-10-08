# sweetsoursong — to do

**Key dates**

- None recorded yet.

## Manuscript

### In progress

- [ ] Co-author read-through of the rewritten manuscript

### Next

- [ ] Decide whether to repeat the founder-control scan at other e_B, c_B∅ values
- [ ] Ask Chris to confirm that the closed metacommunity assumes effectively infinitely many plants (the Lerch et al. ergodicity clause is in Methods)
- [ ] Check PNAS Research Report limits: 4 medium graphical elements (the main text has 6 figures and no tables), 4,000 words, 50 references
- [ ] Check that the Discussion's argument matches the new model (partly revised on 2026-10-05, `sweetsoursong-ms` `5620296`)
- [ ] Review the Significance Statement (revised on Overleaf 2026-10-05) with co-authors
- [ ] If PNAS needs word counts, recount with `texcount -inc -total __ms.tex` (the title-page note was removed)

### Blocked

### Done

- [x] Fig 1E lists processes, not separate transitions: row 3 is total pollinator immigration, c(m P_∅R + m P_YR + m_B P_BR)/N, and the caption says regional colonization (5, 6) is part of it (user, on Overleaf, 2026-10-08; not yet pushed to GitHub when noted)
- [x] Methods: Figure 1E reference fixed; $a_Y$, $a_B$ described as $N$ times the fill rates (`sweetsoursong-ms`, 2026-10-08)
- [x] Write the full model description in the SI; SI title updated (`sweetsoursong-ms`, 2026-10-08)
- [x] Data and code statement cites the Zenodo concept DOI, 10.5281/zenodo.15113987 (`sweetsoursong-ms` `0d49704`, 2026-10-08)
- [x] Rename `_scripts/05-si-table-s1.R` to `05-si-table-s2.R`, with its outputs, to match the SI numbering (2026-10-08)
- [x] New Fig 1 with parameters (D) and one-plant transitions (E); K-plant table to SI (Table S1); `\varnothing` notation; R package mentioned in Methods (user, `a8ebbfc`, 2026-10-08)
- [x] Add Table 2, transitions in the K-plant metacommunity, after Table 1 (`sweetsoursong-ms` `30657e9`; check `claude-checks/multi_plant_table_check.R`) (2026-10-05)
- [x] Accept all tracked changes; remove the hidden old equations and label shims (user, `5620296`, 2026-10-05)
- [x] Add the Lerch et al. citations and the founder-control sentence to Methods; revise the Significance Statement (user, on Overleaf, 2026-10-05)
- [x] Check for regional founder control (both R < 1): none for m_B ≥ m; region starts at m_B = 0.00908 (`claude-checks/founder_control_*`, 2026-10-05)
- [x] Fix the Overleaf word count (TeXcount errors from the hidden Mathematica block; skip deleted and replaced text), commit `1937795` (2026-10-04)
- [x] Numerical checks of the invasion criteria: ε, Pmax, `Chop`, R_B definition (2026-10-02)
- [x] Rewrite Methods, Results, captions, and SI for the new model, with tracked changes (2026-10-02)
- [x] Fix the figure issues: Fig 2A axis label, Fig 4E P_R value, Fig 1 ∅_j label (2026-10-02)
- [x] Draft the Significance Statement; remove the stale word count (2026-10-02)
- [x] Restore the Overleaf–GitHub sync after the history rewrite (2026-10-02)
- [x] Scaffold project context files; commit and push them with the check scripts and outputs (2026-10-02)

## New-model code

### Next

- [ ] Decide whether to also archive Chris's Mathematica files; if so, get final versions of `definitions_unified.wl` and `exploration_unified_10x.nb` and replace the hard-coded `~/Projects/lucas-codex/*.mx` paths
- [ ] Delete `~/GitHub/Stanford/sweetsoursong-ms-backup.bundle` (79 MB, pre-rewrite backup from 2026-10-02; no longer needed)

### Blocked

- [ ] Send `note_for_chris.md` — **blocked on:** your decision whether and how to send it (since 2026-10-02)

### Done

- [x] GitHub release `v2.0.0` and Zenodo archive; the README badge carries the DOI (user, 2026-10-05)
- [x] Review the R port, merge `new-model` into `main`, and push (2026-10-05)
- [x] Chris approved the R port (2026-10-05)
- [x] Manuscript figure panels stay as the Mathematica exports for now; the R scripts reproduce them (2026-10-05)
- [x] Port the new model to R on branch `new-model` (package 2.0.0): 346 tests against Mathematica pass; scripts reproduce Figs 2-6 and Table S1 (`claude-checks/r-port/verify_manuscript_numbers_output.txt`); old package replaced, kept at tag `v1.0.0` (2026-10-05)
- [x] Remove `Chop` (R port default; `options(sweetsoursong.chop = TRUE)` reproduces Chris's code) (2026-10-05)
- [x] Decide what happens to the old R package: replaced on `new-model`, kept at tag `v1.0.0` (2026-10-05)
