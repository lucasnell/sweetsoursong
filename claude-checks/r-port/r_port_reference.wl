(* Reference values for testing the R port of Chris's model code.
   Usage: wolframscript -file r_port_reference.wl chop     (Chris's code as is)
          wolframscript -file r_port_reference.wl nochop   (Chop -> Identity)
   Output: tests/testthat/fixtures/ref_<mode>.json and
           claude-checks/r-port/r_port_reference_<mode>_output.m
   Parameters: SetParameters from exploration_unified_10x.nb.
   definitions_unified.wl is unchanged. Notebook settings (FindRoot starts,
   PR values, continuation steps) are copied from notebook_inputs.txt. *)

mode = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "chop"];
If[! MemberQ[{"chop", "nochop"}, mode], Print["mode must be chop or nochop"]; Exit[1]];

root = "/Users/lucasnell/GitHub/Stanford/sweetsoursong";
SetDirectory[root <> "/claude-checks/chris-files"];
Get["definitions_unified.wl"];

SetParameters := (n = 50; Pmax = 12; c = 500; d = 0.1; {m, mB} = {0.01, 0.05};
  {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5; \[Epsilon] = 10^-5;);

(* run body with Chop disabled in nochop mode *)
SetAttributes[withMode, HoldAll];
withMode[body_] := If[mode === "nochop", Block[{Chop = Identity}, body], body];

out = <||>;
t0 = AbsoluteTime[];
log[s_] := Print[Round[AbsoluteTime[] - t0, 0.1], "s  ", s];

withMode[
  ClearParameters; SetReducedEcoEvoModel; SetParameters; SetCTMCClosureTM;

  (* ---- 1. full-model rates fF at selected states ---- *)
  rateStates = {{10, 20, 3}, {0, 0, 0}, {49, 1, 12}, {25, 25, 1}, {0, 50, 5}};
  ratePools = {20., 15., 15.}; (* P0R, PYR, PBR *)
  out["rates"] = Table[
    Block[{Y = s[[1]], B = s[[2]], P = s[[3]], P\[EmptySet]R = ratePools[[1]],
        PYR = ratePools[[2]], PBR = ratePools[[3]]},
      <|"Y" -> s[[1]], "B" -> s[[2]], "P" -> s[[3]],
        "P0R" -> ratePools[[1]], "PYR" -> ratePools[[2]], "PBR" -> ratePools[[3]],
        "rates" -> N[Values[fF]], "names" -> Keys[fF]|>],
    {s, rateStates}];
  log["rates"];

  (* ---- 2. vacancy fill probabilities ---- *)
  vacCases = {{0, 49, {25., 25.}}, {10, 39, {25., 25.}}, {49, 0, {25., 25.}},
    {25, 24, {25., 25.}}, {0, 49, {10.^-5, 50. - 10.^-5}}, {49, 0, {150. - 10.^-5, 10.^-5}}};
  out["vacancy"] = Table[
    Module[{v = CTMCClosureReplacementProbabilityVectors[vc[[1]], vc[[2]], vc[[3]]],
        sys = BuildCTMCClosureVacancySystem[vc[[1]], vc[[2]], vc[[3]]]},
      <|"y" -> vc[[1]], "b" -> vc[[2]], "pyr" -> vc[[3, 1]], "pbr" -> vc[[3, 2]],
        "hY" -> v["Y"], "hB" -> v["B"],
        "yRepl" -> sys["YReplacementRates"], "bRepl" -> sys["BReplacementRates"],
        "pDown" -> sys["PDownRates"], "pUp" -> sys["PUpRates"]|>],
    {vc, vacCases}];
  log["vacancy"];

  (* ---- 3. hazards and stationary pi(y) at several pools ---- *)
  hazCases = {{1., 25.}, {0.6, 10.^-5}, {2.8, 2.8*50 - 10.^-5}, {1.5, 30.}};
  out["hazards"] = Table[
    Block[{PR = hc[[1]]},
      Module[{rows = CTMCClosureDiagnostics[{hc[[2]], n PR - hc[[2]]}], dist},
        dist = GetYDistribution[{hc[[2]], n PR - hc[[2]]}];
        <|"PR" -> PR, "pyr" -> hc[[2]],
          "up" -> Lookup[rows, "Up"], "down" -> Lookup[rows, "Down"],
          "tail" -> Lookup[rows, "P0TailMass"],
          "piY" -> Table[PDF[dist, y], {y, 0, n}],
          "eYP" -> InOutPY[hc[[2]]]|>]],
    {hc, hazCases}];
  log["hazards"];

  (* ---- 4. invasion criteria ---- *)
  invPR = {0.3, 0.6, 1., 1.5, 1.8, 2.6, 2.8, 3., 5.};
  out["inv"] = Table[Block[{PR = pr}, <|"PR" -> pr, "mB" -> mB, "invY" -> InvY, "invB" -> InvB,
      "eYPdirectB" -> Module[{dY = GetYDistribution[{n PR - \[Epsilon], \[Epsilon]}]},
         Expectation[(n - y) PMeanGivenY[y, {n PR - \[Epsilon], \[Epsilon]}], y \[Distributed] dY]]|>],
    {pr, invPR}];
  out["inv_mBeqm"] = Block[{mB = m}, Table[Block[{PR = pr},
      <|"PR" -> pr, "mB" -> mB, "invY" -> InvY, "invB" -> InvB|>], {pr, {1.45, 1.5, 1.55}}]];
  log["inv"];

  (* thresholds, as in notebook 2.1 and 2.5 *)
  Clear[PR];
  out["thresholds"] = <|
    "invY_mB05" -> (PR /. FindRoot[InvY == 1, {PR, 0.7}]),
    "invB_mB05" -> (PR /. FindRoot[InvB == 1, {PR, 2.9}]),
    "invY_mB01" -> Block[{mB = m}, PR /. FindRoot[InvY == 1, {PR, 1.45}]],
    "invB_mB01" -> Block[{mB = m}, PR /. FindRoot[InvB == 1, {PR, 1.55}]]|>;
  log["thresholds"];

  (* ---- 5. closed-metacommunity examples (notebook 2.2) ---- *)
  exCases = {{0.6, 25}, {1., 25}, {1.8, 85}, {2.6, 130}, {3., 130}};
  out["examples"] = Table[
    Block[{PR = ec[[1]]},
      Module[{eq, pyrStar, dist, decomp, dat, idat, markers, wsc, domain, xmin, ymin},
        eq = FindRoot[InOutPY[pyr] == pyr, {pyr, ec[[2]]}];
        pyrStar = pyr /. eq;
        dist = GetPDF[pyrStar];
        decomp = DecomposeDistribution[dist];
        (* watershed labels, same steps as DecomposeDistribution *)
        domain = dist["Domain"];
        {xmin, ymin} = domain[[{1, 2}, 1]];
        dat = Normal@SparseArray[Map[Plus[#, {1 - xmin, 1 - ymin}] &, dist[[2, 2]]] -> dist[[2, 1]]];
        idat = Image[-dat];
        markers = MinDetect[idat];
        wsc = WatershedComponents[idat, markers, Method -> "Basins"];
        <|"PR" -> PR, "start" -> ec[[2]], "pyr" -> pyrStar,
          "invY" -> InvY, "invB" -> InvB,
          "mean" -> Mean[dist],
          "joint" -> Table[PDF[dist, {y, p}], {y, 0, n}, {p, 0, Pmax}],
          "domainY" -> domain[[1]], "domainP" -> domain[[2]],
          "datDims" -> Dimensions[dat],
          "markers" -> ImageData[markers],
          "labels" -> wsc,
          "weights" -> decomp[[1]],
          "compMeans" -> (Mean /@ decomp[[2]])|>]],
    {ec, exCases}];
  log["examples"];

  (* ---- 6. Fig 5 continuation (notebook 2.3) ---- *)
  Block[{PR = 1}, eq0 = FindRoot[InOutPY[pyr] == pyr, {pyr, 25}]];
  contRow[] := Module[{eq, dist, md},
    eq = FindRoot[InOutPY[pyr] == pyr, {pyr, pyr0}];
    pyr0 = pyr /. eq;
    dist = GetPDF[pyr0];
    md = DecomposeDistribution[dist];
    <|"PR" -> PR, "pyr" -> pyr0, "mean" -> Mean[dist], "weights" -> md[[1]],
      "compMeans" -> (Mean /@ md[[2]]), "invY" -> InvY, "invB" -> InvB|>];
  pyr0 = pyr /. eq0;
  resUp = Table[contRow[], {PR, 1., 2.89, 0.01}];
  pyr0 = pyr /. eq0;
  resDown = Table[contRow[], {PR, 0.99, 0.74, -0.01}];
  Clear[PR];
  out["fig5"] = SortBy[Join[resUp, resDown], #["PR"] &];
  log["fig5"];

  (* ---- 7. deterministic one-plant model (notebook 1.4, Fig 2 B-E) ---- *)
  detCases = {{"noimm", 0.01}, {"noimm", 0.015}, {"noimm", 0.05}, {"reduced", 0.05}};
  out["fig2"] = Table[
    Module[{eq, res},
      ClearParameters;
      If[dc[[1]] === "noimm", SetNoImmigrationEcoEvoModel, SetReducedEcoEvoModel];
      SetParameters; mB = dc[[2]];
      {PYR, PBR} = 2 n {0.9, 0.1};
      eq = SelectValid[SolveEcoEq[]];
      res = <|"model" -> dc[[1]], "mB" -> dc[[2]],
        "eqY" -> (Y /. eq), "eqP" -> (P /. eq),
        "stable" -> (EcoStableQ /@ eq),
        "eigen" -> (ReIm /@ Flatten[EcoEigenvalues /@ eq])|>;
      Clear[PYR, PBR];
      res],
    {dc, detCases}];
  ClearParameters; SetReducedEcoEvoModel; SetParameters; SetCTMCClosureTM;
  log["fig2"];

  (* ---- 8. full model (notebook 3.2, 3.3) ---- *)
  SetFullTM;
  out["full"] = Table[
    Block[{PR = pr}, <|"PR" -> pr, "invY" -> FullInvY, "invB" -> FullInvB|>],
    {pr, {0.6, 1., 2.6, 3.}}];
  log["full inv"];
  Block[{PR = 1.},
    Module[{eq, dist, decomp},
      eq = FindRoot[FullInOutYB[{pyr, pbr}] == {pyr, pbr}, {pyr, 18}, {pbr, 31}];
      dist = FullGetPDF[{pyr, pbr} /. eq];
      decomp = DecomposeDistribution[dist];
      out["full_example"] = <|"PR" -> PR, "pyr" -> (pyr /. eq), "pbr" -> (pbr /. eq),
        "mean" -> Mean[dist], "weights" -> decomp[[1]],
        "compMeans" -> (Mean /@ decomp[[2]])|>]];
  log["full example"];
  Clear[PR];
  out["full_thresholds"] = <|
    "invY" -> (PR /. FindRoot[FullInvY == 1, {PR, 0.74}, Evaluated -> False]),
    "invB" -> (PR /. FindRoot[FullInvB == 1, {PR, 2.85}, Evaluated -> False])|>;
  log["full thresholds"];
];

out["meta"] = <|"mode" -> mode, "version" -> $Version, "date" -> DateString["ISODate"],
  "params" -> <|"n" -> 50, "Pmax" -> 12, "c" -> 500, "d" -> 0.1, "m" -> 0.01, "mB" -> 0.05,
    "eY" -> 1, "eB" -> 0.5, "cB0" -> 5, "eps" -> 10.^-5|>|>;

(* EcoStableQ returns True/False; JSON handles both. Replace any non-numeric leftovers. *)
clean = out /. {x_Real :> x, x_Rational :> N[x], Indeterminate -> "NaN", ComplexInfinity -> "Inf"};
Export[root <> "/tests/testthat/fixtures/ref_" <> mode <> ".json", clean, "JSON", "Compact" -> True];
Put[out, root <> "/claude-checks/r-port/r_port_reference_" <> mode <> "_output.m"];
log["done"];
