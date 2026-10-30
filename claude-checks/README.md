# Numerical checks

Scripts behind statements in the manuscript that are not shown in a figure,
and comparisons of the R package with the Mathematica code it was ported
from. Each script writes the output file named beside it. R scripts run
from the project root (`Rscript claude-checks/<script>.R`) with the
package installed. Wolfram scripts run from any directory with
`wolframscript -file claude-checks/<script>.wl`; they need Wolfram 15.0.1
and the EcoEvo package 1.7.2 (<https://github.com/cklausme/EcoEvo>), and
they load Christopher Klausmeier's functions from
`chris-files/definitions_unified.wl`. Default parameters throughout
(Methods; `ss_params()` in R, `SetParameters` in Mathematica) unless a
script says otherwise. "Reduced model" is the model without uncolonized
flowers, with the vacancy closure, used for all main-text results.


## Statements in the manuscript

| Statement | Where | Script | Output |
|---|---|---|---|
| Thresholds P_R = 0.73 and 2.90; Fig 3 pools; Fig 4 outcomes; Fig 5 coexistence range and pollinator means; Fig 2 equilibria; Fig 6 outcome counts and the alternative-stable-state curve; SI Table S2 values | Results, figure legends, SI Table S2 | `r-port/verify_manuscript_numbers.R` (reads the CSVs made by `_scripts/`) | `r-port/verify_manuscript_numbers_output.txt` |
| No preference (m_B = m): mutual invasion only for 1.46 < P_R < 1.57; P_crit = 1 | Results, abstract | `ms_numbers.wl` | `ms_numbers_output.m` |
| Founder control (neither species invades) requires m_B < 0.0091, none for m_B ≥ m | Methods | `founder_control_scan.wl`, `founder_control_tip.wl`; in R, `founder_control_tip.R` and `r-port/check_founder_grid.R` | `founder_control_scan_output.m`, `founder_control_tip_output.m`, `founder_control_tip_output.txt`, `r-port/check_founder_grid_output.txt` |
| Thresholds do not depend on ε; magnitudes of R_Y and R_B do | Methods, SI Numerical accuracy | `eps_ctmc.wl` (reduced model); `_scripts/05-si-table-s2.R` | `eps_ctmc_output.txt`; `_data/si-table-s2.csv` |
| Thresholds do not depend on P_max; probability of P = P_max and change in pool outputs at P_max = 18 | Methods, SI Numerical accuracy | `pmax_check.R` (reduced and full models); `eps_pmax_full.wl` (full model) | `pmax_check_output.txt`; `eps_pmax_full_output.m` |
| Full-model rows of SI Table S2; the reduced model shifts thresholds by 0.001 and 0.03 | SI Table S2 | `_scripts/05-si-table-s2.R`; `eps_pmax_full.wl` | `_data/si-table-s2.csv`; `eps_pmax_full_output.m` |
| Reduced and full models give similar distributions (pools, mode weights, within-mode means) | SI Numerical accuracy | `full_vs_reduced.R` | `full_vs_reduced_output.txt`, `full_vs_reduced_output.csv` |
| R_B as the invader's own output and as n P_R − E[YP], related by R_B − 1 = (m / m_B)(R_B′ − 1) | SI Invasion criteria | `invB_check.wl`; `r-port/compare_invb.R` | `invB_check_output.txt`; `r-port/compare_invb_output.txt` |
| K-plant transitions (SI Table S1) sum to the one-plant transitions (Fig 1E) | SI Full model description | `multi_plant_table_check.R` | `multi_plant_table_check_output.txt` |


## Other checks

- `fig1e_row3_alternative.R` → `fig1e_row3_alternative_output.txt`: if
  regional colonization (Fig 1E rows 5 and 6) did not also add a
  pollinator, thresholds would shift by −0.0047 (R_Y) and −0.030 (R_B).
  Supports the wording of the Fig 1E legend.
- `chop_threshold_check.txt`: Mathematica's `Chop` in Chris Klausmeier's
  full-model code biases its thresholds (0.729677 vs 0.731249 for R_Y).
  The values were computed interactively after `eps_pmax_full.wl`, so no
  script reproduces this file; `r-port/compare_chop.R` covers the same
  question for the reduced model.


## `chris-files/`

Copies of Christopher Klausmeier's code, unchanged:
`definitions_unified.wl` (model functions used by every Wolfram script
here) and `vacancy_ctmc_closure_note.pdf` (derivation of the vacancy
closure). His notebook `exploration_unified_10x.nb`, which made the
manuscript figures, is not included.


## `r-port/`: the R package against Mathematica

- `r_port_reference.wl` → `tests/testthat/fixtures/ref_chop.json`,
  `ref_nochop.json` and `r_port_reference_*_output.m`: reference values
  for the package tests (see the main README).
- `compare_chop.R` → `compare_chop_output.txt`: Figs 3–6 with and without
  `Chop`.
- `compare_invb.R` → `compare_invb_output.txt`: Fig 6 with either
  definition of R_B.
- `export_founder_grid.wl` → `founder_control_grid.csv`, then
  `check_founder_grid.R` → `check_founder_grid_output.txt`: R against
  Mathematica at all 7381 points of the founder-control scan.
- `verify_manuscript_numbers.R` → `verify_manuscript_numbers_output.txt`:
  numbers printed in the manuscript against the R outputs.
- `extract_notebook_inputs.wl` → `notebook_inputs.txt`: the input cells of
  Chris Klausmeier's notebook as text, from which the figure settings in
  `_scripts/` were taken. Needs the notebook, which is not included.
- `ref_chop.log`, `ref_nochop.log`: logs of the fixture runs.
