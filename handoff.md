# Handoff — 2026-10-05

## Session topic

Manuscript upkeep while the user edits on Overleaf (2026-10-04 to 10-05), then a port of Chris Klausmeier's Mathematica model to R. Earlier in the session: fixed the Overleaf word count (`sweetsoursong-ms` commit `1937795`), placed the Lerch et al. (2023) citations, checked for regional founder control, drafted an intro paragraph replacing commented-out old text, and reorganized the project notes to the agentic-starter templates. The user has since pasted the Lerch et al. text, the founder-control sentence, and a revised Significance Statement into Overleaf.

The R port was developed on branch `new-model` and merged into `main` (pushed 2026-10-05). Package `sweetsoursong` 2.0.0 replaces the old package; the old one stays at tag `v1.0.0`. Plan: `~/.claude/plans/what-ways-could-we-polymorphic-pretzel.md`. User choices: pure R, R package with testthat, replace in this repo, reproduce Chris's code exactly first and then fix known issues.

## Key decisions

- R port defaults differ from Chris's code in two ways, each in its own commit: no `Chop` (`23dbe04`) and R_B = E[BP]/ε (`fa5e098`). Neither changes Figs 2–6 or any threshold: `claude-checks/r-port/compare_chop_output.txt`, `compare_invb_output.txt`. `options(sweetsoursong.chop = TRUE)` and `inv_b(method = "chris")` reproduce his code.
- Two numerical fixes not in Chris's code: truncated Poisson weights in log space (`2d7eaa8`; Mathematica switches to extended precision on underflow, R did not, which put R_Y off by ~3× at m_B ≥ 0.79, P_R ≥ 7.9), and with P_BR exactly 0 all mass at y = N (`503417e`; his product formula zeroes that absorbing state, but his `FindRoot` never lands exactly on the bound).
- Closed-metacommunity roots: Newton from the start value (as `FindRoot`), accepted if interior; otherwise the nearest candidate among sign changes and the bounds. This reproduces all of Chris's roots, including starts on the bound (P_R = 2.6, start 130).
- Watershed: Mathematica's `WatershedComponents[..., Method -> "Basins"]` here equals steepest descent over 8 neighbours with raw differences; R reproduces its labels cell for cell.
- Founder control (both R < 1) does not occur for m_B ≥ m at default parameters; it starts at m_B = 0.00908, P_R = 1.528 (`claude-checks/founder_control_*`).

## Open follow-ups

- [x] User reviewed the R port; merged into `main` and pushed (2026-10-05)
- [x] Chris approved the port
- [x] Manuscript figure panels stay as the Mathematica exports for now
- [ ] User: GitHub release `v2.0.0` and Zenodo archive (user handles these)
- [ ] Update the data and code statement once the release has a DOI
- [ ] Update the title-page word-count note in `__ms.tex` (user said they will do this)

## Context for the next session

- Verification, all saved in `claude-checks/r-port/`:
  - 346 testthat checks against Mathematica fixtures (`tests/testthat/fixtures/ref_chop.json`, `ref_nochop.json`, from `r_port_reference.wl`) pass; about 2 minutes with `devtools::test()`. Full-model tests compare at 1e-6 because `FindRoot` stops at about 8 digits on the resident equilibrium.
  - `check_founder_grid_output.txt`: R_Y matches Mathematica at all 7381 founder-control grid points to 5e-14; outcomes identical.
  - `verify_manuscript_numbers_output.txt`: SI Table S1, Fig 3 caption pools, Fig 5 coexistence range, and the m_B = m window all match the manuscript to printed digits.
  - `R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 note ("New submission").
- Figure scripts make the computed panels only; Fig 2A and Fig 5A are schematics, and Illustrator annotations are not reproduced.
- Without Chop, the full model's 3-D watershed found 363 spurious modes in Mathematica (round-off noise ~1e-19 in tails); in R, clipping negative round-off to 0 removes them. `decompose_distribution(floor =)` is a further safeguard.
- `wolframscript`: set `WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel`; run one heavy script at a time (two in parallel crashed once); `Print` output appears only when the script exits; `log` is protected in that context.
- Sparse stationary solves: replacing a row with ones makes `Matrix::lu` ~100× slower (dense row); the code fixes one state instead, chosen by inflow/outflow ratio, then re-solves at the most probable state.
- Notebook input cells, readable without a front end: `claude-checks/r-port/notebook_inputs.txt` (`extract_notebook_inputs.wl`).
- Overleaf is likely ahead of GitHub for `sweetsoursong-ms`; ask the user to push from Overleaf before any local edit there.
- Earlier handoffs are in git history (`git log -p handoff.md`).
