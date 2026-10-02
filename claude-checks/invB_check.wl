(* Compare reduced-model InvB as coded (n PR - E[y P]) with the direct B-carrier output E[(n-y) P].
   CTMC closure, parameters from SetParameters in exploration_unified_10x.nb. Run 2026-10-02. *)
SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks/chris-files"];
Get["definitions_unified.wl"];
n = 50; Pmax = 12; c = 500; d = 0.1; {m, mB} = {0.01, 0.05}; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5; \[Epsilon] = 10^-5;
SetCTMCClosureTM;
res = Table[
  Module[{dY, pyr = n PR - \[Epsilon], pbr = \[Epsilon]},
   dY = GetYDistribution[{pyr, pbr}];
   {PR, InvB, Expectation[(n - y) PMeanGivenY[y, {pyr, pbr}], y \[Distributed] dY]/\[Epsilon]}],
  {PR, {0.6, 1.0, 2.6, 2.8, 2.9}}];
Print[TableForm[res, TableHeadings -> {None, {"PR", "InvB (code)", "E[(n-y) P]/eps"}}]];
