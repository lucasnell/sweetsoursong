# Checks on definitions_unified.wl and exploration_unified_10x.nb

Parameters throughout are `SetParameters` from the notebook (n = 50, c = 500, d = 0.1, m = 0.01, mB = 0.05, eY = 1, eB = 0.5, cB∅ = 5). Nothing in `definitions_unified.wl` was changed; all checks are in `claude-checks/`; the copies of Chris's files used are in `claude-checks/chris-files/`.

## Summary

| Issue | Thresholds affected? | Magnitudes affected? |
|---|---|---|
| ε = 10⁻⁵ not in the linear regime | No | Yes, by up to >100× |
| Pmax = 12 cutoff | No | Yes, for Y invading at high PR |
| `Chop` in `FullGetPDF` | Yes, ~0.2% in PR | Yes, near thresholds |
| Reduced-model `InvB` definition | No | Differs from E[PB]/ε by a fixed factor |

## 1. ε dependence

Thresholds are independent of ε; invasion-criterion magnitudes are not.

| Model | ε | PR at InvY = 1 | PR at InvB = 1 |
|---|---|---|---|
| CTMC closure | 10⁻⁵, 10⁻⁷, 10⁻⁹ | 0.730855 (all three) | 2.895810–2.895811 |
| Full | 10⁻⁵, 10⁻⁶, 10⁻⁷ | 0.7296768 (all three) | 2.867489 (all three) |

CTMC-closure magnitudes:

| PR | InvY, ε = 10⁻⁵ | 10⁻⁷ | 10⁻⁹ | InvB, ε = 10⁻⁵ | 10⁻⁷ | 10⁻⁹ |
|---|---|---|---|---|---|---|
| 1.0 | 1699.5 | 1699.7 | 1699.7 | 4.97×10⁶ | 4.97×10⁸ | 4.97×10¹⁰ |
| 1.8 | 2.69×10⁷ | 7.27×10⁷ | 7.40×10⁷ | 3.83×10⁶ | 6.67×10⁶ | 6.72×10⁶ |
| 2.6 | 6.09×10⁷ | 3.38×10⁹ | 7.42×10⁹ | 69.22 | 69.22 | 69.22 |

- Where the value grows as ε falls, the invader is not rare in the stationary distribution at ε = 10⁻⁵. At PR = 1.8, InvY × ε = 269, larger than the total pool n·PR = 90.
- InvB at PR = 1 scales as 1/ε: the invader takes over the stationary distribution at every ε tested.
- Only the comparison with 1 is robust. The `{InvY, InvB}` values printed at the example PRs should not be reported as magnitudes.
- **Trap:** with `Chop` left on, ε ≤ 10⁻⁷ gives CTMC InvY = 0 at every PR tested. The vacancy-filling probability at y = 0 is about 0.02 ε, which `Chop` sets to zero. The table above was computed inside `Block[{Chop = Identity}, …]`.

Script: `eps_ctmc.wl` → `eps_ctmc_output.txt`; full model `eps_pmax_full.wl` (+ `eps_pmax_full_defs.wl`) → `eps_pmax_full_output.m`. The full-model runs use a direct sparse solve of TM·π = 0 instead of Arnoldi + `Chop`. It reproduces `FullInvY`/`FullInvB` at PR = 1 to 5×10⁻⁵ relative. Its noise floor (min π ≈ −4×10⁻⁹) limits ε to ≥ 10⁻⁷.

## 2. Pmax cutoff (full model)

Stationary mass at P = Pmax, as recommended in section 10 of the vacancy note:

| Distribution | π(P = 12) | π(P = 18) |
|---|---|---|
| Coexistence equilibria, PR = 1, 1.8, 2.6 | 1.6–3.8×10⁻⁵ | ~2×10⁻⁹ |
| B invading Y monoculture, PR = 0.6–3 | ≤ 5.5×10⁻⁵ | ≤ 3×10⁻⁹ |
| Y invading B monoculture, PR = 1.8 | 0.042 | 0.0013 |
| Y invading B monoculture, PR = 2.6 | 0.21 | 0.034 |
| Y invading B monoculture, PR = 3.0 | 0.28 | 0.072 |

- Equilibrium pool outputs change by ≤ 0.03% between Pmax = 12 and 18.
- FullInvY at high PR depends on the cutoff (PR = 2.6: 4.74×10⁷ → 5.86×10⁷; PR = 3: 5.01×10⁷ → 6.54×10⁷) and is not converged at Pmax = 18.
- Thresholds: 0.729677 → 0.729670 (Y), 2.867489 → 2.867473 (B).
- Not yet checked: the Poisson tail at Pmax in the CTMC closure for the same Y-invasion inputs. E[P | y] reaches 13.9 at y = 49, PR = 3.

## 3. `Chop` biases the full-model thresholds

`FullGetPDF` applies `Chop` (zeroes |x| < 10⁻¹⁰) to the stationary vector, which removes part of the invader's mass.

| | PR at FullInvY = 1 | PR at FullInvB = 1 |
|---|---|---|
| Notebook (Arnoldi + `Chop`, ε = 10⁻⁵) | 0.731249 | 2.86269 |
| Without `Chop` | 0.729677 | 2.867489 |

At PR = 0.729677, `FullInvY` = 0.9489 as written and 1.0000072 inside `Block[{Chop = Identity}, …]`. The direct solve gives 1.000007, so the discrepancy comes from `Chop`, not Arnoldi. The CTMC-closure thresholds at ε = 10⁻⁵ are unaffected (0.730855 reproduced without `Chop`). Without the bias, the CTMC-closure vs full-model gap is 0.0012 for the Y threshold and 0.028 for the B threshold. The notebook values give 0.0004 and 0.033.

Consequence: the full-model invasion contours in notebook section 3.5 carry this bias. The CTMC-closure figures (Fig. 6) do not.

Output: `chop_threshold_check.txt`.

## 4. Reduced-model `InvB` (no action needed for thresholds)

`InvB = (n PR − InOutPY[n PR − ε])/ε` is not E[PB]/ε, the invader's own pool output used in `FullInvB` and in Lerch et al. (2023). The pollinator balance m·E[PY] + mB·E[PB] = m·pyr + mB·pbr, which holds exactly in the reduced model, gives

E[PB]/ε − 1 = (m/mB)·(InvB − 1).

Both therefore cross 1 at the same PR. Their magnitudes differ (m/mB = 0.2): at PR = 2.6, InvB = 69.2 against E[PB]/ε = 14.6 (full model FullInvB: 10.05). Use E[PB]/ε if reduced and full magnitudes are compared.

Script: `invB_check.wl` → `invB_check_output.txt`.

## Suggested changes

1. Remove `Chop` from `FullGetPDF`, `FullGetYPDF`, `FullGetBPDF`, `GetYDistribution`, and `CTMCClosureReplacementProbabilityVectors`, clipping negative round-off to 0 before building `EmpiricalDistribution`. Re-run section 3.5.
2. If invasion magnitudes are reported anywhere, compute the ε → 0 limit by first-order perturbation of the stationary distribution rather than at finite ε.
3. Set Pmax from the tail check in section 10 of the vacancy note, including the invasion inputs, not only the equilibria.
4. The `DumpSave`/`DumpGet` paths (`~/Projects/lucas-codex/*.mx`) don't exist on my machine; a path relative to `NotebookDirectory[]` would make the bifurcation results portable.
