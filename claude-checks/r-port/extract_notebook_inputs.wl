(* claude-checks directory, from this script's location. *)
checks = FileNameDrop[ExpandFileName[$InputFileName], -2];
(* Extract all section headings and input cells (as InputForm text) from
   exploration_unified_10x.nb, in order, without a front end.
   Output: notebook_inputs.txt. Used to port figure settings to R. *)
nb = Get[checks <> "/chris-files/exploration_unified_10x.nb"];
cells = Cases[nb, Cell[_, "Title" | "Section" | "Subsection" | "Subsubsection" | "Input" | "Text", ___], Infinity];
txt[c_] := Module[{b = First[c], s = c[[2]]},
  If[s === "Input",
    Quiet@Check[ToString[ToExpression[b, StandardForm, HoldComplete], InputForm], "<<parse fail>>"],
    If[StringQ[b], b, Quiet@Check[ToString[ToExpression[b, StandardForm, HoldForm], InputForm], ToString[b, InputForm]]]]];
out = StringRiffle[Table["[" <> c[[2]] <> "] " <> txt[c], {c, cells}], "\n\n"];
Export[checks <> "/r-port/notebook_inputs.txt", out, "Text"];
Print[Length[cells]];
