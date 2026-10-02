(* Epsilon dependence of the reduced CTMC-closure invasion criteria.
   Parameters from SetParameters in exploration_unified_10x.nb.
   Chop (zeroes |x| < 1e-10) is disabled with Block so that small-epsilon
   vacancy-filling probabilities are not truncated. definitions_unified.wl is unchanged. *)

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks/chris-files"];
Get["definitions_unified.wl"];
n = 50; Pmax = 12; c = 500; d = 0.1; {m, mB} = {0.01, 0.05}; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5;
SetCTMCClosureTM;
Clear[PR];

epsList = {10^-5, 10^-7, 10^-9};

(* invasion criteria at fixed PR; Chop on (as in Chris's code) and off *)
ctmcInvTable = Flatten[Table[
   Block[{\[Epsilon] = eps, PR = pr},
    {N@eps, pr, InvY, InvB, Block[{Chop = Identity}, InvY], Block[{Chop = Identity}, InvB]}],
   {eps, epsList}, {pr, {1.0, 1.8, 2.6}}], 1];

(* critical PRs, Chop off *)
ctmcCrit = Table[
   Block[{\[Epsilon] = eps, Chop = Identity},
    {N@eps,
     PR /. FindRoot[InvY == 1, {PR, 0.73}],
     PR /. FindRoot[InvB == 1, {PR, 2.9}]}],
   {eps, epsList}];

SetDirectory["/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks"];
Export["eps_ctmc_output.txt",
  StringJoin[
   "Reduced CTMC closure, Pmax = 12. Columns: eps, PR, InvY, InvB (Chop on), InvY, InvB (Chop off)\n",
   ExportString[ctmcInvTable, "Table"],
   "\n\nCritical PR (Chop off). Columns: eps, PR where InvY = 1, PR where InvB = 1\n",
   ExportString[ctmcCrit, "Table"], "\n"], "Text"];
{ctmcInvTable, ctmcCrit}
