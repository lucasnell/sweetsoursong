# Claim provenance and reproducibility review

Reviewed 2026-10-30. Repo `sweetsoursong` at `551363a` (v2.0.0 tag plus later commits). Manuscript `sweetsoursong-ms` at `b4dafab` (2026-10-29), working tree clean; Overleaf may be ahead. "Run" means I ran it this session; "read" means inferred from reading files.

## Verdict

The numbers in the manuscript are, with one rounding exception, correct and each traces to a script and saved output. The weak points are the links between them:

- Every published figure panel comes from a notebook that is not in the public repo, and the manuscript does not say so.
- The R scripts' numeric outputs are in git-ignored `_data/`, so a reader of the repo sees a verification report but not the numbers it checked.
- The Zenodo archive that the Data statement cites predates several files the SI and README refer to.
- Several claims in the SI and Methods rest only on Mathematica output, or on a hand-written summary in `note_for_chris.md`.
- Two `(SI?)` placeholders remain in the text.

No numerical claim I tested was wrong. One printed value is a rounding artefact (finding 6).

## Claim table

Manuscript locations are `file:line` in `sweetsoursong-ms`. "Public?" means tracked in the public repo. "R" is the R package, "WL" is Mathematica.

| # | Claim | Location | Supporting script | Saved output | Public? | Status |
|---|---|---|---|---|---|---|
| 1 | Parameters N=50, c=500, d=0.1, m=0.01, m_B=0.05, e_Y=1, e_B=0.5, c_B∅=5, ε=1e-5, Pmax=12, P_crit=1 | 05-methods:50; Fig 1D | `R/` `ss_params()`; `ms_numbers.wl` (P_crit) | `ms_numbers_output.m`; run: `ss_params()` returns these | yes | Verified by running (defaults match). P_crit = 5·50/(500·0.5) = 1 by arithmetic |
| 2 | Bacteria exclude yeast at P_R < 0.73; yeast exclude bacteria at P_R > 2.90; coexist between | 03-results:18; 98-figures:146 | `_scripts/05-si-table-s2.R`; `verify_manuscript_numbers.R` | `r-port/verify_manuscript_numbers_output.txt` (0.7308551, 2.895810); `eps_ctmc_output.txt` (WL) | yes | Verified (R and WL agree to 7 digits) |
| 3 | "No preference (m_B = m = 0.01): mutual invasion only for 1.46 < P_R < 1.57" | 03-results:20 | `ms_numbers.wl`; `verify_manuscript_numbers.R` | `ms_numbers_output.m` (1.462006, 1.566834); verify output (R matches) | yes | Verified. The R check compares to constants typed into the script, not to the .m file |
| 4 | Proportion of yeast plants goes 0 to 1 over the range; mean P about 3 on yeast plants, about 0.7 on bacteria plants | 03-results:18 | `_scripts/03-fig5.R` | `_data/fig5-continuation.csv` (git-ignored) | no | Run: p(Y) 0.004 to 0.997; P on Y 2.90 to 3.29; P on B 0.68 to 0.73. Not in any tracked output |
| 5 | Coexistence region lies largely inside the alternative-stable-state region (dashed curves, m_B/m > max(P_R/P_crit, P_crit/P_R)) | 03-results:20; Figs 6, S1 captions | `_scripts/04-fig6.R` (draws the curve) | `_data/fig6-grid-*.rds` (git-ignored) | no | Curve is analytic (SI eq S7). "Largely inside" is a visual judgement; no overlap fraction is computed anywhere |
| 6 | Alternative stable states in Fig 2 (P_R = 2, P_YR:P_BR = 9:1, m_B = 0.01, 0.015, 0.05) | 98-figures:72; 03-results:10-12 | `_scripts/01-fig2.R` | `_data/fig2-equilibria.csv` (git-ignored) | no | Run (log in `_data/`): B: 2 stable states; C: 1; D: 2; E: 2. Hand check of conditions gives the same. Not in tracked output |
| 7 | Fig 3: P_R = 1.8, P_YR = 71.87, P_BR = 18.13; simulations stay near their start for 1000 days | 98-figures:90-95; 03-results:14 | `_scripts/02-fig3-4.R`; `verify_manuscript_numbers.R` | verify output (P_YR, P_BR match); simulation in `_data/fig3-gillespie.rds` | partly | P_YR and P_BR verified. Run: R simulations stay in y 37-50 and 0-14. The published panels used Mathematica random draws, so trajectories differ by design |
| 8 | Fig 4: outcomes and mode weights at P_R = 0.6, 1, 1.8, 2.6, 3 | Fig 4 caption | `_scripts/02-fig3-4.R` | verify output (Fig 4 section); `r-port/compare_chop_output.txt` | yes | Verified against WL fixtures (test suite, other reviewer) |
| 9 | Fig 5 continuation, steps of 0.01 from P_R = 1, up to 2.89 and down to 0.74 | 05-methods:44; SI line 318-320 | `_scripts/03-fig5.R` | verify output: 216 points, all R_Y > 1 and R_B > 1 | yes | Verified |
| 10 | Newton start values 25, 25, 85, 130, 130 for P_R = 0.6, 1, 1.8, 2.6, 3 | SI:316-317 | `examples` table in `_scripts/` (not opened); `tests/testthat/test-examples.R` | fixtures | yes | Inferred by reading; consistent with the notebook inputs file |
| 11 | Fig 6 grids: outcomes for P_R 0.1-10 and m_B up to 0.1 and 1 | Fig 6 and S1 captions | `_scripts/04-fig6.R`; `compare_chop.R` | `compare_chop_output.txt`: counts 2460/1582/1524/0 (A), 4183/6826/1091/0 (B) | yes | Verified. The grid starts at m_B = m, so it cannot show m_B < m |
| 12 | No parameter values with m_B >= m where neither species can invade | 05-methods:47 | `check_founder_grid.R`; `founder_control_scan.wl` | `check_founder_grid_output.txt`; `founder_control_grid.csv` | yes | Verified. Run: the 377 "neither" points have m_B <= 0.00891 |
| 13 | Founder control "requires m_B < 0.0091 (scan over 0.1 <= P_R <= 10)" | 05-methods:47 | `founder_control_tip.wl` | `founder_control_tip_output.m` (m_B = 0.009076) | yes | WL only. R reproduces the grid, not the 0.009076 tip. The grid's m_B steps (61 log-spaced values) are too coarse to give 0.0091 alone |
| 14 | Invasion thresholds independent of ε (1e-5 to 1e-9) | 05-methods:47; SI:386-387 | `eps_ctmc.wl`; `05-si-table-s2.R` | `eps_ctmc_output.txt`; `verify_manuscript_numbers_output.txt` | yes | Verified (R and WL). See finding 6 for one printed digit |
| 15 | Thresholds independent of Pmax (12 to 18), as a claim about the main-text (reduced) model | 05-methods:47; SI:387-388 | none for the reduced model | none | no | Not supported by any saved output. Run: reduced thresholds at Pmax 12, 18, 24 are 0.7308551, 0.7308542, 0.7308542 (Y) and 2.89581 (B) in all three. Table S2 shows Pmax only for the full model |
| 16 | Table S2 entries (5 rows) | SI:411-415 | `05-si-table-s2.R`; `eps_ctmc.wl`; `eps_pmax_full.wl` | `si-table-s2.csv` (git-ignored); `verify_manuscript_numbers_output.txt`; `eps_pmax_full_output.m` | partly | Verified except row 3 (finding 6) |
| 17 | Vacancy closure shifts thresholds by 0.001 (Y) and 0.03 (B) relative to the full model | SI:389-390 | `05-si-table-s2.R` | Table S2 values: 0.730855 - 0.729677 = 0.0012; 2.895810 - 2.867489 = 0.0283 | partly | Arithmetic from the table agrees with rounding (0.0012 reads as "0.001", 0.0283 as "0.03"). The same SI sentence says 0.001 while `note_for_chris.md` says 0.0012; one significant figure for the Y shift hides that it is 20% larger than printed |
| 18 | At coexistence equilibria, probability of P = Pmax at most 4e-5; pool outputs change at most 0.03% for Pmax 18 | SI:396-399 | `eps_pmax_full.wl` (full model) | `eps_pmax_full_output.m`; summary typed by hand into `note_for_chris.md` §2 | yes | The SI does not say these are full-model numbers. Run for the reduced model: mass at P = 12 is 1.7e-5, 4.0e-5, 4.2e-5 (P_R = 1, 1.8, 2.6); at Pmax 18 about 2e-9; P_YR changes by under 1e-5 relative. So "at most 4e-5" is slightly low for the reduced model (4.18e-5 at P_R = 2.6) |
| 19 | Magnitudes of R_Y, R_B not converged in ε (2.7e7 at 1e-5, 7.4e7 at 1e-9 for R_Y at P_R = 1.8) | SI:391-395 | `eps_ctmc.wl` | `eps_ctmc_output.txt` | yes | Verified in WL output (2.692e7, 7.398e7). No R version |
| 20 | Full model has 17,238 states for N = 50, Pmax = 12 | SI:358 | none | none | n/a | Run by arithmetic: 1326 (Y+B <= 50) x 13 = 17,238 |
| 21 | R_B - 1 = (m/m_B)(R_B' - 1); Fig 6 identical under either | SI:343-349 | `invB_check.wl`; `compare_invb.R` | `invB_check_output.txt`; `compare_invb_output.txt` | yes | Verified (read; outputs not re-run) |
| 22 | K-plant Table S1 sums to Fig 1E rates | SI:113-115; 02-modeling-framework:11 | `multi_plant_table_check.R` | `multi_plant_table_check_output.txt` (4.5e-16) | yes | Verified. The script and output are not referenced from the manuscript, only from a LaTeX comment (SI:31) |
| 23 | Alternative reading of Fig 1E (separate events) would shift thresholds by -0.0047 / -0.030 | PROJECT_INDEX decision log only | `fig1e_row3_alternative.R` | `fig1e_row3_alternative_output.txt` | yes | Not in the manuscript; fine |
| 24 | Reduced model "does not meaningfully affect our results (SI?)" | 02-modeling-framework:19 | none cited | none | n/a | Placeholder. The only evidence is the threshold shifts in row 17, which do not address Figs 3-6 distributions |
| 25 | Reduced-model results "closely match the full model (SI?)" | 05-methods:32 | as row 17 | as row 17 | n/a | Placeholder; SI Numerical accuracy section exists and should be cited |
| 26 | "An R package that replicates all results" | 05-methods:61 | `_scripts/` | none | partly | Overclaim; see finding 2 |
| 27 | Data and code archived at Zenodo 10.5281/zenodo.15113987 | __ms.tex:205-206 | n/a | n/a | n/a | Run (curl): the concept DOI redirects to record 23170901, v2.0.0, 2026-10-05, one zip. See finding 3 |

## Findings

### Blocker

None for correctness. Findings 1 to 3 should be fixed before submission because a reviewer who follows the Data statement will hit them.

### Major

**1. The published panels come from a notebook that is not archived, and the manuscript does not say so.**
Evidence: figure files in `sweetsoursong-ms/figures/` are Illustrator-assembled Mathematica exports (decided 2026-10-05). The source, `exploration_unified_10x.nb`, is git-ignored (`.gitignore:53`), so Zenodo v2.0.0 does not contain it. `definitions_unified.wl` is tracked, but the notebook's `DumpGet` caches point to `~/Projects/lucas-codex/*.mx` (`note_for_chris.md`, suggested change 4). 05-methods:61 says the R package "replicates all results", and the Data statement says "data and code supporting the results" are archived. No text tells the reader the published panels are not the R output. No image-level comparison of R and published panels exists: the saved comparisons (`compare_chop`, `compare_invb`, `check_founder_grid`, `verify_manuscript_numbers`, the 346 tests) compare numbers at specific points (equilibria, weights, thresholds, grid outcomes), not panel contents. Panels with random draws (Fig 4 left, Fig 5A illustration, Fig 3 time series) cannot match Mathematica draws.
Fix: (a) change 05-methods:61 to say what is true, for example "The R package reproduces the model, the invasion thresholds, and the computed panels of Figs 2-6 and Table S2; the published panels were drawn from Wolfram Language output, and simulated draws differ"; (b) archive the notebook (after Chris replaces the `.mx` paths) or state it is available on request; (c) add one saved side-by-side comparison, such as overlay of R and notebook Fig 6 outcome grids, or at least a CSV of the notebook's Fig 6 grid against the R grid. A grid-level comparison already exists for the founder-control scan (`check_founder_grid_output.txt`, 7381 points) and could be extended to the Fig 6 grids by exporting the notebook's grids.

**2. The R scripts' own outputs are git-ignored, so the public repo holds the checks but not what was checked.**
Evidence: `_data/` and `_figures/` are ignored. `verify_manuscript_numbers.R` reads `_data/si-table-s2.csv`, `fig3-4-equilibria.csv`, `fig5-continuation.csv`, and `fig6-grid-*.rds`. Claims 4, 5, 6, and 16 rest on these. The tracked `verify_manuscript_numbers_output.txt` repeats the R values, which covers claims 2, 3, 7, 9, 16 but not 4, 5, or 6 (the pollinator means, the P on each plant type, the Fig 2 equilibria). `cached()` in `_scripts/00-preamble.R` reloads an `.rds` whenever the file exists, regardless of package version or parameters, so stale caches can pass silently. The verify script's header says it compared against the manuscript at `origin/main` on 2026-10-05; the manuscript has since changed and the values are hard-coded in the script, so it will not notice future edits.
Fix: track the small CSVs (`fig2-equilibria.csv`, `fig3-4-equilibria.csv`, `fig5-continuation.csv`, `si-table-s2.csv`) by un-ignoring them, or write them into `claude-checks/r-port/`; have `cached()` store `packageVersion()` and the `pars` list and refuse a mismatch; add the Results-text numbers in row 4 to `verify_manuscript_numbers.R`.

**3. The Zenodo archive the Data statement cites is behind the repo and the SI.**
Evidence (run): `git ls-tree v2.0.0` has `_scripts/05-si-table-s1.R` (renamed to `05-si-table-s2.R` on 2026-10-08) and lacks `claude-checks/multi_plant_table_check.R` (cited in the SI comment, SI:31) and `fig1e_row3_alternative.R`. The archived `verify_manuscript_numbers_output.txt` predates the Table S2 renumbering. The README at v2.0.0 also differs. The concept DOI resolves to "latest", so a reviewer today sees v2.0.0 and a later release would silently change what the DOI shows. This is the tradeoff behind the 2026-10-05 decision, which was about writing a DOI into the repo; the manuscript could still cite a version DOI.
Fix: cut a v2.0.1 (or similar) release before submission, after un-ignoring or tracking the outputs from finding 2. Then either keep the concept DOI and accept that it drifts, or cite the version DOI in the Data statement. This is a choice for you; the concept DOI is convenient for authors and the version DOI is what reviewers can check against the text.

**4. Claims supported by Mathematica output only, or by a hand summary.**
- Claim 13 (m_B < 0.0091): WL only (`founder_control_tip.wl`). No R reproduction of the tip.
- Claim 15 (Pmax independence for the reduced model): nothing saved. My run supports it (above), but no script in the repo does.
- Claim 18: the summary figures (1.6 to 3.8e-5, 0.03%) are typed into `note_for_chris.md` from `eps_pmax_full_output.m`; no script computes them, and the SI does not label them as full-model.
- Claim 19: WL only.
- Claim 22: not cited in the manuscript.
Fix: add one R script, for example `claude-checks/r-port/pmax_reduced_check.R`, that saves the reduced-model Pmax thresholds and the P = Pmax mass (my scratch versions are in the scratchpad directory and take about a minute). Add "full model" to SI:396. Add R reproduction of the 0.00908 tip with `uniroot` over (P_R, m_B).

**5. Two unresolved `(SI?)` placeholders in the main text.**
Evidence: `02-modeling-framework.tex:19` ("does not meaningfully affect our results (SI?)") and `05-methods.tex:32` ("(SI?)"). Both should point to the "Numerical accuracy" subsection and Table S2. The first claim is wider than the evidence: the SI compares thresholds only (0.001 and 0.03 shifts), not the stationary distributions in Figs 3 to 5. A 0.03 shift in the bacteria threshold (1% of 2.9) is small but it is a shift in a headline number (2.90).
Fix: 02:19 "...and shifts the invasion thresholds by 0.001 (yeast) and 0.03 (bacteria) (SI Appendix, Table S2)"; 05-methods:32 "(SI Appendix, Table S2)". Delete "meaningfully" unless distributions are compared.

### Minor

**6. Table S2, reduced model, ε = 1e-9, R_B = 1 prints 2.89581; the computed value is 2.895811003.**
Evidence: `eps_ctmc_output.txt` lines 17-18: the 6-digit print "2.89581" came from a rerun after a failed FindRoot (the file records the failed "2.9" entry just above it). R gives 2.895811003 (`verify_manuscript_numbers_output.txt`, `05-si-table-s2.log`). At 1e-5 and 1e-7 the table prints 6 decimals (2.895810, 2.895811), so 2.89581 reads as 2.895810. In `verify_manuscript_numbers.R` the check for this row uses `dig_b = 5` where all other rows use 6, which makes it print "match"; at 6 digits it would be DIFFERENT. This is a rounding display issue, not a numerical one.
Fix: print 2.895811 in the table; set `dig_b` to 6 in the script and rerun it.

**7. SI:386-388 states ranges that hold for one model.**
"Changed by less than 1e-5 across ε = 1e-5 to 1e-9" holds for the reduced model; the full model was run for 1e-5 to 1e-7 only (noise floor limits it, `note_for_chris.md` §1). "Less than 2e-5 for Pmax 12 to 18" holds for the full model (7e-6 and 1.6e-5) and, by my run, for the reduced model (9e-7 and under 1e-6). Fix: name the model in each clause.

**8. The "largely inside" and "neither sufficient nor necessary" statements (Fig 6, S1 captions; 03-results:20) have no computed support.**
They read off the figure. A one-line count (fraction of coexist cells inside the dashed region; number of exclusion cells inside it; number of coexist cells outside) from the existing grid would make them checkable.

**9. `claude-checks/` is hard to use as a reproducibility record for an outsider.**
- No README or index in `claude-checks/` (the only map is the Sources comments in `99-supplement.tex:372-379` and `note_for_chris.md`). Script-to-claim mapping exists only in my table above.
- Eight `.wl` scripts hard-code `/Users/lucasnell/GitHub/Stanford/sweetsoursong/...` (`eps_ctmc.wl:6`, `eps_pmax_full.wl:10,14`, `founder_control_scan.wl:11`, `founder_control_tip.wl:6`, `invB_check.wl:3`, `ms_numbers.wl:7`, `export_founder_grid.wl:4`, `extract_notebook_inputs.wl:4`, `r_port_reference.wl:13`). They cannot run elsewhere without edits; they need Wolfram 15.0.1 and EcoEvo 1.7.2.
- Software versions are recorded for the R scripts (package version, sometimes R version) but not in `.wl` outputs, except `invB_check_output.txt` ("Wolfram 15.0.1"). `eps_ctmc_output.txt` has no header with date or version. Several outputs have the date only in the filename of a note.
- `eps_ctmc_output.txt` retains a failed FindRoot value ("2.9") with a comment; correct but confusing.
- Stale text: `R/rates.R:3` and `R/full.R:17` refer to "Table 1", which no longer exists (now Fig 1E and SI eq S2); `multi_plant_table_check.R` header and comments call the K-plant table "Table 2" and the one-plant table "Table 1". `.Rbuildignore` lists `_mathematica`, which does not exist.
- `compare_chop.R`, `compare_invb.R`, and the verify script depend on `_data/` being present.
Fix: add `claude-checks/README.md` with a table (script, output, software, claim); replace the absolute paths with `NotebookDirectory[]` or a `root` variable at the top of each script; add a date and version header to each `.wl` output; update the Table 1/2 wording.

**10. Things tracked in a public repo that probably should not be, or need consent.**
Tracked (`git ls-files`): `claude-checks/cover-letter/make_tracked.py` (script for the PNAS cover letter, with the docx names; no letter text), `claude-checks/ms-edit/trackdel.py` (tooling for the manuscript's tracked-changes round), `claude-checks/note_for_chris.md` and `chris-files/` (Chris's `definitions_unified.wl` and closure-note PDF), `handoff.md`, `CLAUDE.md`, `PROJECT_INDEX.md`, `TODO.md`.
- `note_for_chris.md` is a critique of Chris's code (it reports a `Chop` bias of about 0.2% in P_R in his full-model thresholds) that TODO.md lists as unsent ("blocked on your decision whether to send"). It is public and inside the Zenodo v2.0.0 archive, so he can read it there before you decide. This is the most sensitive item.
- `chris-files/` republishes a co-author's code and a PDF note under the repo's GPL-3 licence. Chris is listed as author in `DESCRIPTION`, but I found no statement that he agreed to publishing the Mathematica source or the note.
- `cover-letter/`, `ms-edit/`, `handoff.md`, `TODO.md`, `PROJECT_INDEX.md`, `CLAUDE.md` are working scaffolding, with Box paths, the PNAS cover letter workflow, and user decisions. They do no harm but they are noise for a reviewer, and `.Rbuildignore` correctly keeps them out of the package.
- Currently untracked (so not public): `claude-checks/reviewers/suggested_reviewers.md` (reviewer suggestions with conflict checks and gender notes) and `claude-checks/param-justification/` (Zotero full-text keyword output). Keep both out of git; add them to `.gitignore` or move them out of the repo, because a later `git add -A` would publish them and they would reach the next Zenodo archive.
Fix: ask Chris about `chris-files/` and `note_for_chris.md` before the next release; move `cover-letter/`, `ms-edit/`, `reviewers/`, `param-justification/` out of the repo (or ignore them) before the release in finding 3.

**11. Model consistency: no mismatches found, with two documentation gaps.**
See "Checked and fine". The gaps: Fig 1D lists neither ε nor Pmax, which the Methods and Table S2 use (Methods 05:50 states Pmax; ε appears only in 05:47 text); and the Fig 1E caption says regional colonization rows 5 and 6 are part of row 3, while their "Change" column shows only Y→Y+1 and B→B+1 (no P change), so a reader taking the table literally still gets the separate-event chain whose thresholds shift by -0.0047 and -0.030. The caption warns about this; a "P → P+1 also" note in the row would remove the ambiguity.

### Nit

- `05-methods:50` cites "T. Fukami, unpublished results" for parameter values; the `param-justification/` folder (untracked) suggests an attempt to find published support. Commented-out justification text remains at 05-methods:57-64.
- `verify_manuscript_numbers.R` heads "sweetsoursong-ms, origin/main on 2026-10-05"; update when rerun.
- Fig 3 caption shows 71.87 and 18.13 while R gives 71.8733 and 18.1267; both round correctly.
- `_data/` contains `.mx` and old figure files (`Bplot-u*.mx`, `invasion-*.pdf`) from the previous model; they are ignored but confusing if someone zips the folder.
- The R DESCRIPTION maintainer email differs from the one on your account; not a reproducibility issue.

## Checked and fine

- Parameter defaults in `ss_params()` equal the manuscript values (run).
- SI eq S2 (per-plant transitions) matches `flower_rates()` term by term (read): `p_plus` is the non-colonizing arrival rate, `py_plus` and `pb_plus` are the colonizing arrivals, with the same factors as SI eq S2.
- Fig 1E rows 3 to 6 are consistent with SI eq S2: immigration row 3 equals the S2 `P immigration` rate plus rows 5 and 6 (read; the algebra `m <PY>(1 - e_Y ∅/N) + m e_Y <PY> ∅/N = m <PY>` holds).
- SI eq S6 and the vacancy rates S7 (with ∅ = 1) match the reduced-model fill rates a_Y, a_B in Methods eq. 1 (a_Y/N and a_B/N) and Fig 1E (read).
- P_crit = c_B∅ N / [c (e_Y - e_B)] = 1 for the defaults (arithmetic).
- Fig 2 conditions P_R m_R/m_B < P_crit < P_R m_R/m for m_B = 0.01, 0.015, 0.05 at P_R = 2 with a 9:1 pool, giving 0.56, 1.4, and 2.8 times P_crit boundaries; the equilibria in `_data/01-fig2.log` agree (e.g. D: y = 0, 26.45, 50 with P = 0.56, 0.97, 2.8).
- Table S2 matches R outputs in four of five rows at the printed digits, and the Y-threshold column in all five.
- Newton continuation in Fig 5 (216 points, 0.74 to 2.89) all have two modes and both invasion criteria above 1 (run, `_data/fig5-continuation.csv`).
- Data statement DOI resolves (run, curl) to a v2.0.0 record containing `sweetsoursong-v2.0.0.zip`.
- R and Mathematica agree on the founder-control grid at all 7381 points (output file, read; not rerun).
- The SI and manuscript state the same model and parameter list; no stale old-model equations remain in the files I read.

## Could not check

- Any Wolfram or Mathematica result: not run, per instructions. Mathematica-only claims (13, 19, the WL half of 2, 14, 18) were checked against the saved `.m` and `.txt` outputs and against R outputs where they exist.
- The 11 MB notebook (git-ignored); I read neither it nor the published figure panels beyond Fig 1, which I read as an image.
- Whether the published panels 2 to 6 match the R figures visually. I did not generate the R figures or compare PDFs.
- The 346 tests, the `renv` install path, and `_scripts/` as a clean run: the other reviewer's scope. I used the existing `_data/` caches (dated 2026-10-05 and 2026-10-08), which I did not regenerate, so they confirm what the scripts produced then, not that the scripts reproduce from scratch.
- Whether Overleaf is ahead of `b4dafab`.
- Whether Chris has agreed to the public copies of his files.
- `eps_pmax_full_output.m` internals beyond the tail I read (the 0.03% and 4e-5 figures).
- Intro (`01-intro.tex`) and Discussion claims beyond what I grepped; they are verbal, not quantitative.
