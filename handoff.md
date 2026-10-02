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
