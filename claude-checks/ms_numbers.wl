(* Numbers cited in the manuscript text that are not printed on the figures.
   Reduced model, CTMC closure, SetParameters from exploration_unified_10x.nb
   (n = 50, Pmax = 12, c = 500, d = 0.1, m = 0.01, mB = 0.05, eY = 1, eB = 0.5,
   cB0 = 5, eps = 1e-5). Chop disabled (see note_for_chris.md, section 1).
   definitions_unified.wl is unchanged. *)

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks/chris-files"];
Get["definitions_unified.wl"];
n = 50; Pmax = 12; c = 500; d = 0.1; {m, mB} = {0.01, 0.05}; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5; \[Epsilon] = 10^-5;
SetCTMCClosureTM;
Clear[PR];

out = <||>;

(* (1) invasion criteria at the PR values labelled in Fig. 4 panel E (2.8) and used in the notebook (3) *)
out["fig4E"] = Table[Block[{PR = pr, Chop = Identity}, {pr, InvY, InvB}], {pr, {2.8, 3.0}}];

(* (2) no pollinator preference, mB = m: critical PRs and invasion criteria at PR = 1.5 *)
out["mBeqm"] = Block[{mB = m, Chop = Identity},
  {PR /. FindRoot[InvY == 1, {PR, 1.45}],
   PR /. FindRoot[InvB == 1, {PR, 1.55}],
   Block[{PR = 1.5}, {InvY, InvB}]}];

(* (1b) B invader's own pool output E[(n - y) P]/eps, which avoids the affine rescaling in InvB
   (InvB as coded can be negative; see note_for_chris.md, section 4) *)
out["fig4E_InvB_direct"] = Table[Block[{PR = pr},
   Module[{pyr = n PR - \[Epsilon], pbr = \[Epsilon], dY},
    dY = GetYDistribution[{pyr, pbr}];
    {pr, InvB, Expectation[(n - y) PMeanGivenY[y, {pyr, pbr}], y \[Distributed] dY]/\[Epsilon]}]],
  {pr, {2.8, 2.9, 3.0}}];

(* (3) Pcrit for the closed one-plant model *)
out["Pcrit"] = cB\[EmptySet] n/(c (eY - eB));

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks"];
Put[out, "ms_numbers_output.m"];
out
