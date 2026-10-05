# sweetsoursong — to do

**Key dates**

- None recorded yet.

## Manuscript

### In progress

- [ ] Read through the tracked-change rewrite (Methods, Results, captions, SI) on Overleaf

### Next

- [ ] Decide whether to repeat the founder-control scan at other e_B, c_B∅ values
- [ ] Ask Chris whether the closed metacommunity assumes effectively infinitely many plants (needed for the optional Lerch et al. ergodicity clause, if it was added)
- [ ] Check PNAS Research Report limits: 4 medium graphical elements (the main text has 6 figures and 1 table), 4,000 words, 50 references
- [ ] Revise the Discussion beyond its figure and wording swaps; its argument still follows the old multi-plant model
- [ ] Review the Significance Statement (revised on Overleaf 2026-10-05) with co-authors
- [ ] Accept all tracked changes (`\usepackage[final]{changes}`) once co-authors have reviewed
- [ ] Update the title-page word-count note in `__ms.tex` (still 3870 / 149). Recount: the Methods and Significance Statement changed on 2026-10-05. Recount with `texcount -inc -total __ms.tex` after accepting the changes

### Blocked

- [ ] Update the data and code statement (Zenodo DOI) — **blocked on:** new-model code being archived (since 2026-10-02)

### Done

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

- [ ] Delete `~/GitHub/Stanford/sweetsoursong-ms-backup.bundle` (79 MB, pre-rewrite backup from 2026-10-02; no longer needed)
- [ ] Get the final versions of `definitions_unified.wl` and `exploration_unified_10x.nb` from Chris and add them to a public repo or archive
- [ ] Remove `Chop` from the stationary-distribution helpers (see `claude-checks/note_for_chris.md`, change 1)
- [ ] Replace the hard-coded `~/Projects/lucas-codex/*.mx` paths with paths relative to `NotebookDirectory[]`
- [ ] Decide what happens to the old R package and `_scripts/`, which reproduce none of the current figures

### Blocked

- [ ] Send `note_for_chris.md` — **blocked on:** your decision whether and how to send it (since 2026-10-02)

### Done

