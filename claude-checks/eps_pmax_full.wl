(* Full three-species model: (1) mass at the pollinator cutoff P = Pmax, and
   (2) dependence of the invasion criteria on epsilon and Pmax.
   Parameters from SetParameters in exploration_unified_10x.nb.
   Invader distributions use a direct sparse solve of TM.pi = 0 (no Arnoldi, no Chop),
   so probability mass of order epsilon is retained. Residents use Chris's
   FullInOutB / FullInOutY (monoculture chains; their mass is not small).
   Direct-solve noise floor: min(pi) about -4e-9, so eps is kept >= 1e-7.
   definitions_unified.wl is unchanged. *)

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks/chris-files"];
Get["definitions_unified.wl"];
n = 50; c = 500; d = 0.1; {m, mB} = {0.01, 0.05}; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5;

Get["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks/eps_pmax_full_defs.wl"];

(* coexistence equilibria reported in the notebook (Pmax = 12), section 3.3 *)
coexPools = {{1.0, {18.2848, 31.304}}, {1.8, {71.8926, 17.6553}}, {2.6, {125.109, 4.39406}}};

out = <||>;
Do[
  Pmax = pm; Block[{PR}, SetFullTM];

  (* sanity: direct solver vs Chris's FullInvY at eps = 1e-5, PR = 1 *)
  If[pm == 12,
   out["sanity"] = Block[{PR = 1., \[Epsilon] = 10^-5}, {FullInvY, invYd[10^-5][[1]], FullInvB, invBd[10^-5][[1]]}]];

  (* (1) P = Pmax mass at coexistence pools and at invasion inputs *)
  out[{"coex", pm}] = Table[
    Block[{PR = cp[[1]]}, Module[{s = statFull[cp[[2]]]}, {cp[[1]], pmaxMass[s], inOutDirect[s]}]],
    {cp, coexPools}];
  out[{"inv", pm}] = Table[
    Block[{PR = pr}, {pr, invYd[10^-5], invBd[10^-5]}],
    {pr, {0.6, 1.0, 1.8, 2.6, 3.0}}];

  (* (2) critical PRs for each epsilon *)
  out[{"crit", pm}] = Table[
    {N@eps,
     pr /. FindRoot[gY[pr, eps] == 1, {pr, 0.72, 0.74}],
     pr /. FindRoot[gB[pr, eps] == 1, {pr, 2.84, 2.88}]},
    {eps, {10^-5, 10^-6, 10^-7}}];
  , {pm, {12, 18}}];

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks"];
Put[out, "eps_pmax_full_output.m"];
out
