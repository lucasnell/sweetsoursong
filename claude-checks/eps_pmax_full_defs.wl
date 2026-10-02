ClearAll[statFull, inOutDirect, pmaxMass, residentB, residentY, invYd, invBd];

(* stationary distribution of the full chain, direct solve; returns {pi, states, residual, min pi} *)
statFull[{pyr_?NumericQ, pbr_?NumericQ}] := Module[{A, A2, len, v, states},
  A = SparseArrayReplace[TM, {PYR -> pyr, PBR -> pbr, P\[EmptySet]R -> n PR - pyr - pbr}];
  len = Length[A];
  A2 = A; A2[[1]] = ConstantArray[1., len];
  v = LinearSolve[A2, N@UnitVector[len, 1]];
  states = Values[as];
  {v, states, Norm[A . v, Infinity], Min[v]}];

inOutDirect[s_] := {s[[1]] . (s[[2]][[All, 1]] s[[2]][[All, 3]]),
                    s[[1]] . (s[[2]][[All, 2]] s[[2]][[All, 3]])};
pmaxMass[s_] := Total@Pick[s[[1]], s[[2]][[All, 3]], Pmax];

residentB[] := pbr /. FindRoot[FullInOutB[pbr] == pbr, {pbr, 0.99 PR n}];
residentY[] := pyr /. FindRoot[FullInOutY[pyr] == pyr, {pyr, 0.99 PR n}];

(* invasion criteria with explicit epsilon; also return P = Pmax mass of the invasion distribution *)
invYd[eps_] := Module[{s = statFull[{eps, residentB[]}]}, {inOutDirect[s][[1]]/eps, pmaxMass[s], s[[3]], s[[4]]}];
invBd[eps_] := Module[{s = statFull[{residentY[], eps}]}, {inOutDirect[s][[2]]/eps, pmaxMass[s], s[[3]], s[[4]]}];

ClearAll[gY, gB];
gY[pr_?NumericQ, eps_] := Block[{PR = pr}, invYd[eps][[1]]];
gB[pr_?NumericQ, eps_] := Block[{PR = pr}, invBd[eps][[1]]];

