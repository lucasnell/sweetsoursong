(* Is there any (PR, mB) where neither species can invade (R_Y < 1 and R_B < 1,
   regional founder control)? Fig 6 (notebook section 2.4) colours only the other three
   outcomes, over PR in [0.1, 10] and mB in [m, 1].
   Reduced model, CTMC closure, SetParameters from exploration_unified_10x.nb
   (n = 50, Pmax = 12, c = 500, d = 0.1, m = 0.01, eY = 1, eB = 0.5, cB0 = 5, eps = 1e-5).
   Chop disabled (see note_for_chris.md, section 1). InvB as coded has the same sign of
   R - 1 as the direct B-carrier output (invB_check_output.txt), so InvB is used here.
   Grid: 121 log-spaced PR in [0.1, 10] x 61 log-spaced mB in [0.001, 1].
   definitions_unified.wl is unchanged. *)

dir = "/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks";
LaunchKernels[];
ParallelEvaluate[
  SetDirectory[dir <> "/chris-files"];
  Get["definitions_unified.wl"];
  n = 50; Pmax = 12; c = 500; d = 0.1; m = 0.01; {eY, eB} = {1, 0.5}; cB\[EmptySet] = 5; \[Epsilon] = 10^-5;
  SetCTMCClosureTM;
  Clear[PR, mB];
  inv[pr_, mb_] := Block[{PR = pr, mB = mb, Chop = Identity}, {InvY, InvB}];
];

prs = 10^Subdivide[-1., 1., 120];
mbs = 10^Subdivide[-3., 0., 60];

grid = ParallelTable[{pr, mb, Sequence @@ inv[pr, mb]}, {mb, mbs}, {pr, prs},
   Method -> "FinestGrained"];

class[{_, _, ry_, rb_}] := Which[
   ry < 1 && rb > 1, "bacteria",
   ry > 1 && rb > 1, "coexist",
   ry > 1 && rb < 1, "yeast",
   True, "neither"];

flat = Flatten[grid, 1];
neither = Select[flat, class[#] === "neither" &];

(* per mB row: number of sign changes of R_Y - 1 and R_B - 1 along PR *)
crossings = Table[{row[[1, 2]],
    Count[Differences[Sign[row[[All, 3]] - 1]], Except[0]],
    Count[Differences[Sign[row[[All, 4]] - 1]], Except[0]]}, {row, grid}];

out = <|
  "counts" -> Counts[class /@ flat],
  "countsInFig6Range" -> Counts[class /@ Select[flat, #[[2]] >= 0.01 - 10^-12 &]],
  "neitherMBRange" -> If[neither === {}, None, MinMax[neither[[All, 2]]]],
  "neitherPRRange" -> If[neither === {}, None, MinMax[neither[[All, 1]]]],
  "neither" -> neither,
  "crossingsPerRow" -> crossings,
  "nonFinite" -> Select[flat, ! AllTrue[#[[3 ;; 4]], NumericQ] &],
  "grid" -> grid|>;

Put[out, dir <> "/founder_control_scan_output.m"];
Print[KeyDrop[out, {"grid", "neither", "crossingsPerRow"}]];
Print["rows with >1 crossing: ", Select[crossings, #[[2]] > 1 || #[[3]] > 1 &]];
