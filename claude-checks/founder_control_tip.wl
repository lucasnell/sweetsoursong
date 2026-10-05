(* Follow-up to founder_control_scan.wl: (1) the mB = m row of the scan, against the
   thresholds in ms_numbers_output.m (1.462, 1.567); (2) the point where R_Y = 1 and
   R_B = 1 cross, which is the tip of the neither-invades region; (3) the smallest
   |log10 R| in the scan, as a check that underflow warnings could not flip a sign.
   Same model and parameters as founder_control_scan.wl. *)
dir = "/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks";
SetDirectory[dir <> "/chris-files"];
Get["definitions_unified.wl"];
n = 50; Pmax = 12; c = 500; d = 0.1; m = 0.01; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5; \[Epsilon] = 10^-5;
SetCTMCClosureTM; Clear[PR, mB];
iy[pr_?NumericQ, mb_?NumericQ] := Block[{PR = pr, mB = mb, Chop = Identity}, InvY];
ib[pr_?NumericQ, mb_?NumericQ] := Block[{PR = pr, mB = mb, Chop = Identity}, InvB];
scan = Get[dir <> "/founder_control_scan_output.m"];
row = SelectFirst[scan["grid"], Abs[#[[1, 2]] - 0.01] < 10^-9 &];
cls = Function[r, Which[r[[3]] < 1 && r[[4]] > 1, "B", r[[3]] > 1 && r[[4]] > 1, "C", r[[3]] > 1 && r[[4]] < 1, "Y", True, "N"]];
out = <||>;
out["mBeqmRowCoexistPR"] = MinMax[Select[row, cls[#] === "C" &][[All, 1]]];
tip = FindRoot[{iy[pr, mb] == 1, ib[pr, mb] == 1}, {pr, 1.5}, {mb, 0.0095}];
out["tip"] = {pr, mb} /. tip;
out["tipCheck"] = {iy @@ ({pr, mb} /. tip), ib @@ ({pr, mb} /. tip)};
flat = Flatten[scan["grid"], 1];
out["minAbsLog10R"] = {Min[Abs[Log10[flat[[All, 3]]]]], Min[Abs[Log10[Abs[flat[[All, 4]]]]]]};
out["nonPositiveRB"] = Count[flat[[All, 4]], x_ /; x <= 0];
SetDirectory[dir];
Put[out, "founder_control_tip_output.m"];
Print[out];
