(* Export the grid of founder_control_scan_output.m (Mathematica, Chop off,
   Chris's InvB) to CSV for comparison with the R port.
   Output: founder_control_grid.csv with columns pr, mB, invY, invB. *)
dir = "/Users/lucasnell/GitHub/Stanford/sweetsoursong/claude-checks";
scan = Get[dir <> "/founder_control_scan_output.m"];
Export[dir <> "/r-port/founder_control_grid.csv",
  Prepend[Flatten[scan["grid"], 1], {"pr", "mB", "invY", "invB"}], "CSV"];
Print[Length[Flatten[scan["grid"], 1]]];
