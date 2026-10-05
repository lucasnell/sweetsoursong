# Handoff — 2026-10-05

## Session topic

Manuscript upkeep while the user edits on Overleaf (2026-10-04 to 10-05). Fixed the Overleaf word count, which TeXcount broke on the hidden Mathematica block in `02-methods.tex` (commit `1937795` in `sweetsoursong-ms`). Worked out where to cite Lerch et al. (2023), `Lerch2023` in `refs.bib`, in Methods. Checked whether regional founder control (neither species can invade) occurs, since the Methods classify only three outcomes and Fig 6 colours only three. Retrieved the first, narrower Significance Statement draft from git history.

## Key decisions

- TeXcount skips the hidden Mathematica block and the text in `\deleted`, the old argument of `\replaced`, and `\comment`. Recorded in `CLAUDE.md` and the decision log.
- The user edits in Overleaf and pastes suggested text there; nothing in `sweetsoursong-ms` was edited locally after `1937795`. A local update to the title-page word-count note was undone at the user's request.
- No regional founder control for m_B ≥ m at default parameters (`claude-checks/founder_control_scan.wl` → `founder_control_scan_output.m`, `founder_control_tip.wl` → `founder_control_tip_output.m`, 2026-10-05). The R_Y = 1 and R_B = 1 curves cross at m_B = 0.00908, P_R = 1.528; both R < 1 occurs only below that, for m_B ≤ 0.0089 and 1.12 ≤ P_R ≤ 4.82 on the grid. Each R crosses 1 at most once along P_R. Not varied: e_Y, e_B, c_B∅, N.

## Open follow-ups

- [ ] User: paste into the "Closed metacommunity" subsection of `02-methods.tex` on Overleaf (text below): the Lerch et al. citations and the founder-control sentence
- [ ] Confirm with Chris that the closed metacommunity assumes effectively infinitely many plants, before adding the ergodicity clause
- [ ] Decide whether to repeat the founder-control scan at other e_B and c_B∅ values
- [ ] Update the title-page word-count note in `__ms.tex` (still 3870 / 149). On 2026-10-04: main text 3065, abstract 157, Significance Statement 104

## Context for the next session

- Overleaf is likely ahead of GitHub. Before any local edit to `sweetsoursong-ms`, ask the user to push from Overleaf (Menu → GitHub), then `git fetch` and fast-forward. Push only with the user's say-so.
- Suggested text, all inside existing `\added{}` blocks so no new markup is needed:
  - First sentence of "Closed metacommunity": "We model many plants coupled only through the regional pollinator pool, following the closed-metacommunity model of Lerch et al.~\cite{Lerch2023}." Their eq. 4 sets regional abundance to the mean of the local stationary distribution (our P_YR = E[YP]), solved by root-finding.
  - Optional, after "expectations are over the stationary distribution of one plant given the pool": "which for many plants is also the distribution of plant states across the metacommunity~\cite{Lerch2023}" (their ergodicity argument; pending the check with Chris).
  - After "$R_Y = {\rm E}[Y P] / \varepsilon > 1$": `~\cite{Lerch2023}`. Their criterion λ = log(N̄_j / ε) is log R.
  - After "or yeast win ($R_B < 1 < R_Y$)." (search Overleaf for "classified outcomes"):
    ```latex
    With $m_B \ge m$, there were no parameter values at which neither species
    could invade, the regional founder control of ref.~\cite{Lerch2023};
    this outcome requires $m_B < 0.0091$
    (scan over $0.1 \le P_R \le 10$).
    ```
- The founder-control margin is narrow: the crossing at m_B = 0.00908 is about 9% below m = 0.01, consistent with the narrow no-preference coexistence window (1.46 < P_R < 1.57).
- The first Significance Statement draft (118 words, ends "dispersal--community feedback may be an overlooked mechanism of species coexistence") is in `sweetsoursong-ms` commit `dd877ed`. The broader reframe is `b97dfdc`; the current Overleaf version descends from it.
- `wolframscript` needs `WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel` and the binary at `/Applications/Wolfram.app/Contents/MacOS/wolframscript`. `NotebookImport` with `"InputText"` fails without a front end; read input cells with `Get` on the `.nb` and `ToExpression[boxes, StandardForm, HoldComplete]`.
- Earlier handoffs (2026-10-02, 2026-10-04) are in git history (`git log -p handoff.md`).
