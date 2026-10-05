# Handoff — 2026-10-05

## Session topic

Where to cite Lerch et al. (2023), `Lerch2023` in `refs.bib`, in Methods, and whether regional founder control (neither species can invade) occurs. Nothing was edited in `sweetsoursong-ms`; the user is editing on Overleaf and pastes suggested text there.

## Lerch et al. citations (suggested, for the "Closed metacommunity" subsection of `02-methods.tex`)

- First sentence: "…coupled only through the regional pollinator pool, following the closed-metacommunity model of Lerch et al.~\cite{Lerch2023}." Their eq. 4 sets regional abundance to the mean of the local stationary distribution, which is our P_YR = E[YP], and they solve it by root-finding.
- Optional, after "expectations are over the stationary distribution of one plant given the pool": "which for many plants is also the distribution of plant states across the metacommunity~\cite{Lerch2023}" (their ergodicity argument). Confirm with Chris that the formulation assumes effectively infinitely many plants.
- Invasion criteria: cite after $R_Y = {\rm E}[Y P] / \varepsilon > 1$. Their λ = log(N̄_j / ε), so λ = log R.

## Founder-control check

- Fig 6 (notebook section 2.4) colours only three outcomes; a both-R < 1 region would be left blank.
- `claude-checks/founder_control_scan.wl` → `founder_control_scan_output.m`: P_R in [0.1, 10] × m_B in [0.001, 1], reduced model, vacancy closure, default parameters. For m_B ≥ m = 0.01 no point has both R < 1. Both R < 1 occurs for m_B ≤ 0.0089 and 1.12 ≤ P_R ≤ 4.82. Each R crosses 1 at most once along P_R in every m_B row.
- `claude-checks/founder_control_tip.wl` → `founder_control_tip_output.m`: the R_Y = 1 and R_B = 1 curves cross at m_B = 0.00908, P_R = 1.528, the tip of that region (about 9% below m). The m_B = m row reproduces the 1.46–1.57 coexistence window.
- Not varied: e_Y, e_B, c_B∅, N.
- Suggested sentence, after "or yeast win ($R_B < 1 < R_Y$)." in the "Closed metacommunity" subsection:
  ```latex
  With $m_B \ge m$, there were no parameter values at which neither species
  could invade, the regional founder control of ref.~\cite{Lerch2023};
  this outcome requires $m_B < 0.0091$
  (scan over $0.1 \le P_R \le 10$).
  ```

---

# Handoff — 2026-10-04

## Session topic

Fixed the Overleaf word count, which failed with 27 TeXcount errors ("Reached end of file while waiting for `\]`"). Commit `1937795` in `sweetsoursong-ms`, pushed and pulled into Overleaf; the user confirmed the count works there.

## What was wrong and what changed

- **Cause:** Overleaf's word count runs TeXcount, which does not honour `\iffalse`. The hidden Mathematica block in `02-methods.tex` (the `\iffalse` before `\[EmptySet]sub = …`) contains 23 `\[EmptySet]` tokens, each read as an unclosed display equation. This also swallowed the rest of Methods, so the old count (3282) was too low.
- **Fix 1:** `%TC:ignore` / `%TC:endignore` around that block. The other `\iffalse` blocks in the main text hold only equations and do not affect the count.
- **Fix 2:** the top of `__ms.tex` now has `%TC:macro` lines so TeXcount skips `\deleted{}`, the old argument of `\replaced{}{}`, and `\comment{}`. The count therefore reflects the text with all changes accepted.

## Current counts (2026-10-04, `texcount`)

- Main text: 3065 words (`texcount -inc -total __ms.tex` in `sweetsoursong-ms`). The title page, abstract, acknowledgments, figures, and SI sit in `%TC:ignore` regions.
- Abstract: 157 words; Significance Statement: 104 words (`texcount -sub=section 00-abstract.tex`). PNAS limits are 250 and 120.
- The title-page note in `__ms.tex` still reads "main text = 3870, abstract = 149" inside `\deleted{}`. I drafted a `\replaced{}{}` update and the user asked me to undo it because they are editing on Overleaf; they may update it there.

## Context for the next session

- The user is actively editing on Overleaf, so Overleaf is ahead of GitHub. Before any local edit to `sweetsoursong-ms`, ask the user to push from Overleaf (Menu → GitHub), then `git fetch` and fast-forward.
- Sync order that worked: user pushes from Overleaf → `git fetch`, stash local edits, `git merge --ff-only origin/main`, pop → commit, check `origin/main` is an ancestor of `HEAD`, push (only with the user's say-so) → user pulls in Overleaf.
- `sweetsoursong-ms` was clean at `1937795` on 2026-10-04, matching `origin/main`.

---

# Handoff — 2026-10-02

## Session topic

Checked Chris Klausmeier's new-model code (`claude-checks/note_for_chris.md`), then rewrote the manuscript in `sweetsoursong-ms` for the new pollinator-pool metacommunity model: Methods, Results, captions for Figs 1–6, the SI, and a Significance Statement. Every edit is tracked with the LaTeX `changes` package. A history rewrite to strip figures from the manuscript repo broke Overleaf's GitHub sync and was undone. Last, the project context files were scaffolded.

## Key decisions

- Report only whether each invasion criterion exceeds 1, since magnitudes are not converged in ε. Recorded in `CLAUDE.md`.
- Abstract hedged to coexistence "across a wide range of regional pollinator abundances". With no pollinator preference there is still a narrow window, 1.46 < P_R < 1.57.
- Keep figures tracked in `sweetsoursong-ms` and never rewrite its history. Overleaf can't unlink and ignores `.gitignore`.
- Tracked changes are red: `\usepackage[commentmarkup=uwave,defaultcolor=red]{changes}` plus `\addedcolor`. The user is applying this in Overleaf.

## Open follow-ups

- [ ] Read through the tracked-change rewrite on Overleaf with co-authors
- [ ] Check PNAS limits on graphical elements (6 figures and 1 table against 4 medium elements)
- [ ] Revise the Discussion's argument for the new model; only figure references and wording were swapped
- [ ] Archive the new-model code and update the Zenodo statement
- [ ] Decide whether to send `note_for_chris.md`

## Context for the next session

- `sweetsoursong-ms` is at `69dfd7a` on GitHub and Overleaf. The user may have newer edits on Overleaf, such as the red colour, so run `git pull` there before editing locally.
- The struck-through old text uses label shims that print "Figure old N"; leave them until all changes are accepted.
- A backup of the manuscript repo from before the rewrite is at `~/GitHub/Stanford/sweetsoursong-ms-backup.bundle`. It's no longer needed and can be deleted.
- My permission settings blocked deleting a remote branch, and moving `main` with deletion of filter-repo's replacement refs. The user ran those.
