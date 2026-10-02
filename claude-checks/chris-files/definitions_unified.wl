(* ::Package:: *)

(* Definitions for the sweetsoursong metacommunity models.

Full stochastic model:
  SetFullTM, FullGetPDF[], FullInOutYB[], FullInOutY[], FullInOutB[]

Reduced stochastic closures:
  SetMeanFieldClosureTM, SetPoissonHazardClosureTM, SetCTMCClosureTM
  GetYDistribution[], GetPDF[], InOutPY[]

Deterministic EcoEvo models:
  SetFullEcoEvoModel, SetReducedEcoEvoModel, SetNoImmigrationEcoEvoModel
  
*)


(* ::Section:: *)
(*Setup*)


Needs["EcoEvo`"];


Unprotect[Color];
Color[Y]=RGBColor[0.880722, 0.611041, 0.142051];Color[B]=RGBColor[0.368417, 0.506779, 0.709798];Color[P]=RGBColor[0.5, 0, 0.5];
Protect[Color];


(* ::Section::Closed:: *)
(*Utilities*)


(* inspired by Carl Wohl https://mathematica.stackexchange.com/a/154287/6358 *)
ClearAll[SparseArrayReplace];
SparseArrayReplace[s_SparseArray, rule_] := With[
  {
    elem = Replace[s["NonzeroValues"], rule, -1],
    def = Replace[s["Background"], rule]
  },
  SparseArray[
    Automatic,
    s["Dimensions"],
    def,
    {1, {s["RowPointers"], s["ColumnIndices"]}, elem}
  ]
];


(* based on <https://mathematica.stackexchange.com/a/83404/6358> by MichaelE2 *)
ClearAll[ConstantCoefficient];
ConstantCoefficient[list_List, form_] := (ConstantCoefficient[#, form]) & /@ list;
ConstantCoefficient[expr_, form_] := If[expr =!= 0, Flatten[CoefficientList[expr, form]][[1]], 0];


(* simulate a continuous-time Markov process with the Gillespie algorithm *)
ClearAll[GillespieSSA];
GillespieSSA[
  proc_List,
  ics : ({__Rule} | _Association),
  {mint_?NumberQ, maxt_?NumberQ},
  opts___?OptionQ
] := Module[
  {
    maxsteps, dstep, debug, in, vars, initialValues, rates, rules,
    null, a, b, symRates, rep, stepList, compiled, times, rest
  },

  maxsteps = Evaluate[MaxSteps /. Flatten[{opts, Options[GillespieSSA]}]];
  dstep = Evaluate[DStep /. Flatten[{opts, Options[GillespieSSA]}]];
  debug = Evaluate[Debug /. Flatten[{opts, Options[GillespieSSA]}]];

  stepList = N@RandomVariate[ExponentialDistribution@1, maxsteps + 1];

  in = Normal[ics];
  vars = in[[All, 1]];
  If[debug, Print["vars=", vars]];
  initialValues = in[[All, 2]];

  rules = proc[[All, 1]];
  If[debug, Print["rules=", rules]];

  rates = proc[[All, 2]];
  If[debug, Print["rates=", rates]];

  null = Table[var -> var, {var, vars}];
  a = Table[
    Table[
      Coefficient[Table[var\[Prime] /. Join[rule, null], {var\[Prime], vars}], var],
      {var, vars}
    ],
    {rule, rules}
  ];
  If[debug, Print["a=", a]];
  b = Table[ConstantCoefficient[vars /. Join[rule, null], vars], {rule, rules}];
  If[debug, Print["b=", b]];

  Block[{count},
    rep = Thread[vars -> Table[Indexed[count, i], {i, Length@vars}]];
    symRates = rates /. rep;
    compiled = ReleaseHold[
      Hold@Compile[
        {
          {init, _Integer, 1},
          {aa, _Integer, 3},
          {bb, _Integer, 2},
          {dtList, _Real, 1},
          {min, _Real},
          {max, _Real},
          {iter, _Integer},
          {resol, _Integer}
        },
        Module[
          {
            count = init, rates, i = 1, c = 1, t = min, dt,
            range, r, rateSum, data = Internal`Bag[]
          },
          rates = "SymbolicRates";
          rateSum = N@Total@rates;
          range = Range@Length@rates;
          Internal`StuffBag[data, Internal`Bag[Join[{t}, N@count]]];
          While[rateSum > 0. && t < max && i <= iter,
            i++;
            dt = dtList[[i]]/rateSum;
            t = Min[max, t + dt];
            If[Mod[i, resol] == 0,
              c++;
              Internal`StuffBag[data, Internal`Bag[Join[{t}, N@count]]];
            ];
            r = RandomChoice[rates -> range];
            count = Max[0, #] & /@ (aa[[r]] . count + bb[[r]]);
            rates = "SymbolicRates";
            rateSum = N@Total@rates
          ];

          If[Mod[i, resol] != 0 || (t < max && i < iter),
            c++;
            Internal`StuffBag[data, Internal`Bag[Join[{max}, N@count]]]
          ];

          Table[Internal`BagPart[Internal`BagPart[data, j], All], {j, c}]
        ],
        Parallelization -> True,
        RuntimeAttributes -> Listable,
        RuntimeOptions -> "Speed",
        CompilationOptions -> {
          "InlineExternalDefinitions" -> True,
          "InlineCompiledFunctions" -> True
        }
      ] /. "SymbolicRates" -> symRates
    ];

    {times, rest} = {First@#, Rest@#} & @ Transpose @ compiled[
      Round@initialValues,
      Round@a,
      Round@b,
      N@stepList,
      N@mint,
      N@maxt,
      Round@maxsteps,
      dstep
    ];

    Thread[vars -> (Interpolation[Transpose@{times, #}, InterpolationOrder -> 0] & /@ rest)]
  ]
];

Options[GillespieSSA] = {MaxSteps -> 10^7, DStep -> 1, Debug -> False};
GillespieSSA[proc_List, ics : ({__Rule} | _Association), maxt_?NumberQ, opts___?OptionQ] :=
  GillespieSSA[proc, ics, {0, maxt}, opts];


(* ::Section:: *)
(*Rate Definitions*)


ClearAll[f, fF, fR, proc, Psub, \[EmptySet]sub, no\[EmptySet]sub, $ActiveTM];

$ActiveTM = None;

(* Keep regional-pool symbols symbolic if this file is reloaded mid-analysis. *)
Block[{PYR, PBR},

\[EmptySet]sub = {\[EmptySet] -> n - Y - B}; (* empty flowers as a function of Y and B *)
no\[EmptySet]sub = {\[EmptySet] -> 0, P\[EmptySet]R -> 0, B -> n - Y}; (* assume no empty flowers *)

(* \[EmptySet]-Y-B-P model transition rates *)
f = <|
  "P-" -> c P (m \[EmptySet]/n + m Y/n + mB B/n), (* P emigration *)
  "Y+" -> c P (1 - m) eY Y/n \[EmptySet]/n, (* Y birth *)
  "B+" -> c P (1 - mB) eB B/n \[EmptySet]/n + cB\[EmptySet] B \[EmptySet]/n, (* B birth *)
  "P+" -> c (m P\[EmptySet]R/n + m (1 - eY \[EmptySet]/n) PYR/n + mB (1 - eB \[EmptySet]/n) PBR/n), (* P immigration only *)
  "PY+" -> c m eY PYR/n \[EmptySet]/n, (* P+Y immigration *)
  "PB+" -> c mB eB PBR/n \[EmptySet]/n, (* P+B immigration *)
  "Y-" -> d Y, (* Y death *)
  "B-" -> d B (* B death *)
|>;

(* full three-species model rates *)
fF = Simplify /@ (f /. \[EmptySet]sub);

(* simplified no-empty-flowers reduced-model transition rates *)
With[
  {fillingRate = Total@Lookup[f, {"Y+", "B+", "PY+", "PB+"}]},
  fR = <|
    "P-" -> f["P-"] /. no\[EmptySet]sub, (* P emigration *)
    "Y+" -> Simplify[f["B-"] f["Y+"]/fillingRate] /. no\[EmptySet]sub, (* B death, Y birth (local) *)
    "B+" -> Simplify[f["Y-"] f["B+"]/fillingRate] /. no\[EmptySet]sub, (* Y death, B birth (local) *)
    "PY+" -> Simplify[f["B-"] f["PY+"]/fillingRate] /. no\[EmptySet]sub, (* B death, Y immigration *)
    "PB+" -> Simplify[f["Y-"] f["PB+"]/fillingRate] /. no\[EmptySet]sub, (* Y death, B immigration *)
    "P+" -> Simplify[f["P+"] /. no\[EmptySet]sub] (* P immigration *)
  |>
];

(* Conditional mean of the fast pollinator cloud. *)
Psub = {P -> (mB PBR + m PYR)/(mB (n - Y) + m Y)};

];


(* ::Section:: *)
(*Deterministic EcoEvo Models*)


ClearAll[SetFullEcoEvoModel, SetReducedEcoEvoModel, SetNoImmigrationEcoEvoModel];

SetFullEcoEvoModel := SetModel[{
  Pop[Y] -> {Equation :> fF["Y+"] + fF["PY+"] - fF["Y-"], Color -> Color[Y]},
  Pop[B] -> {Equation :> fF["B+"] + fF["PB+"] - fF["B-"], Color -> Color[B]},
  Pop[P] -> {Equation :> fF["P+"] + fF["PY+"] + fF["PB+"] - fF["P-"], Color -> Color[P]},
  Parameters :> {
    c > 0, d > 0, 0 <= m <= 1, 0 <= mB <= 1, cB\[EmptySet] >= 0,
    0 <= eY <= 1, 0 <= eB <= 1, PR >= 0, P\[EmptySet]R >= 0,
    PYR >= 0, PBR >= 0, n > 0
  }
}];

SetReducedEcoEvoModel := SetModel[{
  Pop[Y] -> {Equation :> Simplify[fR["Y+"] + fR["PY+"] - fR["B+"] - fR["PB+"]], Color -> RGBColor[0.880722, 0.611041, 0.142051]},
  Pop[P] -> {Equation :> Simplify[fR["P+"] - fR["P-"]], Color -> RGBColor[0.5, 0, 0.5]},
  Parameters :> {
    c > 0, d > 0, 0 <= m <= 1, 0 <= mB <= 1, cB\[EmptySet] >= 0,
    0 <= eY <= 1, 0 <= eB <= 1, PR >= 0, PYR >= 0, PBR >= 0,
    n > 0
  }
}];

SetNoImmigrationEcoEvoModel := SetModel[{
  Pop[Y] -> {Equation :> Simplify[fR["Y+"] - fR["B+"]], Color -> RGBColor[0.880722, 0.611041, 0.142051]},
  Pop[P] -> {Equation :> Simplify[fR["P+"] - fR["P-"]], Color -> RGBColor[0.5, 0, 0.5]},
  Parameters :> {
    c > 0, d > 0, 0 <= m <= 1, 0 <= mB <= 1, cB\[EmptySet] >= 0,
    0 <= eY <= 1, 0 <= eB <= 1, PR >= 0, PYR >= 0, PBR >= 0,
    n > 0
  }
}];


(* ::Section:: *)
(*Reduced Closure Helpers*)


ClearAll[
  PMeanGivenY, PoissonAverageP,
  BuildReducedClosureTM,
  SetMeanFieldClosureTM, SetPoissonHazardClosureTM,
  GetYDistribution, GetPDF, PoissonizeP, InOutPY
];

(* Remove definitions left by older versions when this file is reloaded. *)
ClearAll[SetReducedClosureTMFromRates, yUpRate, yDownRate];

PMeanGivenY[y_, {pyr_, pbr_}] := (mB pbr + m pyr)/(mB (n - y) + m y);
PMeanGivenY[y_, pyr_?NumericQ] := PMeanGivenY[y, {pyr, n PR - pyr}];


(* Average a transition hazard over P | Y ~ Poisson[E[P | Y]].

   The no-empty-flower Y hazards in this model are linear fractions in P,
   (a P + b)/(s P + t).  The expectation is computed analytically via

     E[1/(s P + t)] =
       Exp[-\[Lambda]] Hypergeometric1F1[t/s, 1 + t/s, \[Lambda]]/t.

   Keep the raw Hypergeometric1F1 form here: aggressive Simplify can rewrite
   this to generalized Gamma[a, z0, z1], which is fragile inside Compile. *)
PoissonAverageP::nonlin =
  "Expected a transition hazard with numerator and denominator at most linear in P; got `1`.";

PoissonAverageP[expr_Plus] := Module[
  {terms = PoissonAverageP /@ List @@ expr},
  If[AnyTrue[terms, # === $Failed &], $Failed, Total[terms]]
];

PoissonAverageP[expr_] := Module[
  {\[Lambda], num, den, a, b, s, t},
  \[Lambda] = PMeanGivenY[Y, {PYR, PBR}];
  {num, den} = NumeratorDenominator[Together[expr]];

  If[
    ! PolynomialQ[num, P] || ! PolynomialQ[den, P] ||
      Exponent[num, P] > 1 || Exponent[den, P] > 1,
    Message[PoissonAverageP::nonlin, HoldForm[expr]];
    $Failed,

    a = Coefficient[num, P, 1];
    b = num /. P -> 0;
    s = Coefficient[den, P, 1];
    t = den /. P -> 0;

    Which[
      TrueQ[num == 0],
        0,
      TrueQ[s == 0],
        (a \[Lambda] + b)/t,
      True,
        a/s + (b - a t/s) Exp[-\[Lambda]] Hypergeometric1F1[t/s, 1 + t/s, \[Lambda]]/t
    ]
  ]
];


BuildReducedClosureTM::length =
  "Expected up- and down-rate vectors of length `1`; got lengths `2` and `3`.";

BuildReducedClosureTM[upRates_List, downRates_List] := Module[
  {stateCount = n + 1, diagonalRates},

  If[Length[upRates] =!= stateCount || Length[downRates] =!= stateCount,
    Message[BuildReducedClosureTM::length, stateCount, Length[upRates], Length[downRates]];
    Return[$Failed]
  ];

  as = AssociationThread[Range[stateCount], Range[0, n]];
  diagonalRates = -(
    Append[Most[upRates], 0] +
      Prepend[Rest[downRates], 0]
  );

  TM = SparseArray[
    {
      Band[{2, 1}] -> Most[upRates],
      Band[{1, 2}] -> Rest[downRates],
      Band[{1, 1}] -> diagonalRates
    },
    {stateCount, stateCount},
    0.
  ]
];


SetMeanFieldClosureTM := Block[{PYR, PBR},
  Module[{upRate, downRate},
    upRate = Simplify[(fR["Y+"] + fR["PY+"]) /. Psub];
    downRate = Simplify[(fR["B+"] + fR["PB+"]) /. Psub];

    If[
      BuildReducedClosureTM[
        Table[upRate /. Y -> y, {y, 0, n}],
        Table[downRate /. Y -> y, {y, 0, n}]
      ] === $Failed,
      Return[$Failed]
    ];

    proc = {
      {{P -> P - 1}, fR["P-"]},
      {{Y -> Y + 1}, upRate},
      {{Y -> Y - 1}, downRate},
      {{P -> P + 1}, fR["P+"]}
    };

    $ActiveTM = "MeanFieldClosure";
  ]
];


SetPoissonHazardClosureTM := Block[{PYR, PBR},
  Module[{upRate, downRate},
    upRate = Simplify[PoissonAverageP[fR["Y+"] + fR["PY+"]]];
    downRate = Simplify[PoissonAverageP[fR["B+"] + fR["PB+"]]];

    If[
      BuildReducedClosureTM[
        Table[upRate /. Y -> y, {y, 0, n}],
        Table[downRate /. Y -> y, {y, 0, n}]
      ] === $Failed,
      Return[$Failed]
    ];

    proc = {
      {{P -> P - 1}, fR["P-"]},
      {{Y -> Y + 1}, upRate},
      {{Y -> Y - 1}, downRate},
      {{P -> P + 1}, fR["P+"]}
    };

    $ActiveTM = "PoissonHazardClosure";
  ]
];


(* CTMC vacancy closure.

   After a Y or B death opens a single empty site, this closure solves the
   little pollinator-vacancy CTMC exactly over P = 0..Pmax before asking which
   plant type eventually fills the vacancy.  Pollinator turnover uses the
   original model rates directly. *)

ClearAll[
  BuildCTMCClosureVacancySystem,
  CTMCClosureReplacementProbabilityVectors,
  CTMCClosureBirthPDistribution,
  CTMCClosurePoissonWeights,
  CTMCClosureHazards,
  CTMCClosureOpeningBecomesYProbability,
  CTMCClosureDiagnostics,
  CTMCClosureStationaryWeights,
  CTMCClosureRateExpression,
  BuildCTMCClosureTM,
  SetCTMCClosureTM,
  ctmcClosureRows,
  ctmcClosureTailMasses,
  ctmcClosurePool,
  ctmcClosureOptions
];

BuildCTMCClosureVacancySystem::badstate =
  "The post-death state {Y,B} = `1` is outside the y-only domain.";

Options[BuildCTMCClosureVacancySystem] = {
  "Pmax" -> Automatic
};

BuildCTMCClosureVacancySystem[
  y_Integer?NonNegative,
  b_Integer?NonNegative,
  pool : {pyr_?NumericQ, pbr_?NumericQ},
  OptionsPattern[]
] := Module[
  {
    pmax, ps, empty, rules,
    yReplacementRates, bReplacementRates,
    pDownRates, pUpRates, totalRates,
    matrix, row, totalRate
  },

  If[y + b > n,
    Message[BuildCTMCClosureVacancySystem::badstate, {y, b}];
    Return[$Failed]
  ];

  pmax = Round[OptionValue["Pmax"] /. Automatic -> Pmax];
  ps = Range[0, pmax];
  empty = n - y - b;

  rules[p_] := {
    Y -> y, B -> b, P -> p, \[EmptySet] -> empty,
    PYR -> pyr, PBR -> pbr, P\[EmptySet]R -> 0
  };

  yReplacementRates = Table[
    N[f["Y+"] /. rules[p]] + If[p < pmax, N[f["PY+"] /. rules[p]], 0.],
    {p, ps}
  ];
  bReplacementRates = Table[
    N[f["B+"] /. rules[p]] + If[p < pmax, N[f["PB+"] /. rules[p]], 0.],
    {p, ps}
  ];
  pDownRates = Table[
    If[p > 0, N[f["P-"] /. rules[p]], 0.],
    {p, ps}
  ];
  pUpRates = Table[
    If[p < pmax, N[f["P+"] /. rules[p]], 0.],
    {p, ps}
  ];
  totalRates = yReplacementRates + bReplacementRates + pDownRates + pUpRates;

  matrix = ConstantArray[0., {pmax + 1, pmax + 1}];
  Do[
    row = p + 1;
    totalRate = totalRates[[row]];

    If[
      Chop[totalRate] == 0,
      matrix[[row, row]] = 1.;,
      matrix[[row, row]] = totalRate;
      If[p > 0, matrix[[row, row - 1]] = -pDownRates[[row]]];
      If[p < pmax, matrix[[row, row + 1]] = -pUpRates[[row]]];
    ],
    {p, ps}
  ];

  <|
    "PValues" -> ps,
    "Matrix" -> matrix,
    "YReplacementRates" -> yReplacementRates,
    "BReplacementRates" -> bReplacementRates,
    "PDownRates" -> pDownRates,
    "PUpRates" -> pUpRates,
    "TotalRates" -> totalRates,
    "Pmax" -> pmax,
    "Pool" -> pool
  |>
];

Options[CTMCClosureReplacementProbabilityVectors] =
  Options[BuildCTMCClosureVacancySystem];

CTMCClosureReplacementProbabilityVectors[
  y_Integer?NonNegative,
  b_Integer?NonNegative,
  pool : {_?NumericQ, _?NumericQ},
  OptionsPattern[]
] := Module[{system, rhs, sol},
  system = BuildCTMCClosureVacancySystem[
    y, b, pool,
    "Pmax" -> OptionValue["Pmax"]
  ];
  If[system === $Failed, Return[$Failed]];

  rhs = Transpose@{
    system["YReplacementRates"],
    system["BReplacementRates"]
  };
  sol = LinearSolve[system["Matrix"], rhs];

  <|
    "Y" -> Clip[Chop[sol[[All, 1]]], {0., 1.}],
    "B" -> Clip[Chop[sol[[All, 2]]], {0., 1.}]
  |>
];

Options[CTMCClosureBirthPDistribution] =
  Options[BuildCTMCClosureVacancySystem];

CTMCClosureBirthPDistribution[
  y_Integer?NonNegative,
  b_Integer?NonNegative,
  pool : {_?NumericQ, _?NumericQ},
  OptionsPattern[]
] := Module[{system, absorptionRates},
  system = BuildCTMCClosureVacancySystem[
    y, b, pool,
    "Pmax" -> OptionValue["Pmax"]
  ];
  If[system === $Failed, Return[$Failed]];

  absorptionRates =
    system["YReplacementRates"] + system["BReplacementRates"];

  Clip[
    Chop@LinearSolve[
      system["Matrix"],
      DiagonalMatrix[absorptionRates]
    ],
    {0., 1.}
  ]
];

Options[CTMCClosurePoissonWeights] = {
  "Pmax" -> Automatic,
  "NormalizeTruncatedPoisson" -> True
};

CTMCClosurePoissonWeights[
  y_Integer?NonNegative,
  pool : {_?NumericQ, _?NumericQ},
  OptionsPattern[]
] := Module[
  {pmax, \[Lambda], rawWeights, representedMass, weights},

  pmax = Round[OptionValue["Pmax"] /. Automatic -> Pmax];
  \[Lambda] = N@PMeanGivenY[y, pool];
  rawWeights = PDF[PoissonDistribution[\[Lambda]], #] & /@ Range[0, pmax];
  representedMass = Total[rawWeights];
  weights = If[
    TrueQ[OptionValue["NormalizeTruncatedPoisson"]] && representedMass > 0,
    rawWeights/representedMass,
    rawWeights
  ];

  <|
    "Weights" -> weights,
    "Lambda" -> \[Lambda],
    "RepresentedMass" -> representedMass,
    "TailMass" -> Max[0., 1. - representedMass],
    "Pmax" -> pmax
  |>
];

Options[CTMCClosureHazards] = {
  "Pmax" -> Automatic,
  "NormalizeTruncatedPoisson" -> True
};

CTMCClosureHazards[
  y_Integer?NonNegative,
  pool : {_?NumericQ, _?NumericQ},
  OptionsPattern[]
] := Catch[
  Module[
    {
      pmax, normalize, b, p0Info, p0Weights,
      yVectors, bVectors, yFillAfterBDeath, bFillAfterYDeath
    },

    pmax = Round[OptionValue["Pmax"] /. Automatic -> Pmax];
    normalize = OptionValue["NormalizeTruncatedPoisson"];
    b = n - y;

    p0Info = CTMCClosurePoissonWeights[
      y,
      pool,
      "Pmax" -> pmax,
      "NormalizeTruncatedPoisson" -> normalize
    ];
    p0Weights = p0Info["Weights"];

    yVectors = If[
      b > 0,
      CTMCClosureReplacementProbabilityVectors[
        y,
        b - 1,
        pool,
        "Pmax" -> pmax
      ],
      <|"Y" -> ConstantArray[0., pmax + 1], "B" -> ConstantArray[0., pmax + 1]|>
    ];

    bVectors = If[
      y > 0,
      CTMCClosureReplacementProbabilityVectors[
        y - 1,
        b,
        pool,
        "Pmax" -> pmax
      ],
      <|"Y" -> ConstantArray[0., pmax + 1], "B" -> ConstantArray[0., pmax + 1]|>
    ];

    If[yVectors === $Failed || bVectors === $Failed, Throw[$Failed]];

    yFillAfterBDeath = p0Weights . yVectors["Y"];
    bFillAfterYDeath = p0Weights . bVectors["B"];

    <|
      "Y" -> y,
      "B" -> b,
      "Up" -> If[b > 0, N[d b yFillAfterBDeath], 0.],
      "Down" -> If[y > 0, N[d y bFillAfterYDeath], 0.],
      "Net" -> If[b > 0, N[d b yFillAfterBDeath], 0.] - If[y > 0, N[d y bFillAfterYDeath], 0.],
      "YFillAfterBDeath" -> yFillAfterBDeath,
      "BFillAfterYDeath" -> bFillAfterYDeath,
      "P0Mean" -> p0Info["Lambda"],
      "P0TailMass" -> p0Info["TailMass"],
      "Pmax" -> pmax
    |>
  ]
];

CTMCClosureHazards[
  y_Integer?NonNegative,
  pyr_?NumericQ,
  opts : OptionsPattern[]
] := CTMCClosureHazards[y, {pyr, n PR - pyr}, opts];

Options[CTMCClosureOpeningBecomesYProbability] =
  Options[CTMCClosureHazards];

CTMCClosureOpeningBecomesYProbability[
  y_Integer?NonNegative,
  pool : {_?NumericQ, _?NumericQ},
  OptionsPattern[]
] := Module[{row, b = n - y},
  row = CTMCClosureHazards[
    y,
    pool,
    "Pmax" -> OptionValue["Pmax"],
    "NormalizeTruncatedPoisson" -> OptionValue["NormalizeTruncatedPoisson"]
  ];
  If[row === $Failed, Return[$Failed]];

  N[
    (b row["YFillAfterBDeath"] +
      y (1 - row["BFillAfterYDeath"]))/n
  ]
];

CTMCClosureDiagnostics[
  pool : {_?NumericQ, _?NumericQ},
  opts : OptionsPattern[CTMCClosureHazards]
] := Table[CTMCClosureHazards[y, pool, opts], {y, 0, n}];

CTMCClosureDiagnostics[
  pyr_?NumericQ,
  opts : OptionsPattern[CTMCClosureHazards]
] := CTMCClosureDiagnostics[{pyr, n PR - pyr}, opts];

CTMCClosureStationaryWeights[upRates_List, downRates_List] := Module[
  {weights, scale},

  weights = ConstantArray[0., Length[upRates]];
  weights[[1]] = 1.;

  Do[
    weights[[y + 1]] = If[
      downRates[[y + 1]] > 0,
      weights[[y]] upRates[[y]]/downRates[[y + 1]],
      0.
    ];

    scale = Max[weights];
    If[scale > 10.^100, weights /= scale],
    {y, 1, Length[upRates] - 1}
  ];

  If[Total[weights] > 0, weights/Total[weights], weights]
];

CTMCClosureRateExpression[values_List] := Piecewise[
  Table[{values[[y + 1]], Y == y}, {y, 0, Length[values] - 1}],
  0.
];

BuildCTMCClosureTM[
  pool : {pyr_?NumericQ, pbr_?NumericQ},
  opts : OptionsPattern[CTMCClosureHazards]
] := Module[
  {rows, upRates, downRates, pRateRules},

  rows = CTMCClosureDiagnostics[pool, opts];
  upRates = Lookup[rows, "Up"];
  downRates = Lookup[rows, "Down"];

  If[BuildReducedClosureTM[upRates, downRates] === $Failed,
    Return[$Failed]
  ];

  ctmcClosureRows = rows;
  ctmcClosureTailMasses = Lookup[rows, "P0TailMass"];
  ctmcClosurePool = pool;
  pRateRules = {PYR -> pyr, PBR -> pbr};

  (* Gillespie simulations keep the CTMC-derived Y hazards, but let the fast
     pollinator cloud fluctuate as in the other reduced closures. *)
  proc = {
    {{P -> P - 1}, fR["P-"] /. pRateRules},
    {{Y -> Y + 1}, CTMCClosureRateExpression[upRates]},
    {{Y -> Y - 1}, CTMCClosureRateExpression[downRates]},
    {{P -> P + 1}, fR["P+"] /. pRateRules}
  };

  TM
];

SetCTMCClosureTM := (
  ctmcClosureOptions = {};
  Clear[TM, proc, ctmcClosureRows, ctmcClosureTailMasses, ctmcClosurePool];
  $ActiveTM = "CTMCClosure";
);


GetYDistribution::activetm =
  "The active transition matrix is `1`. Run SetMeanFieldClosureTM, SetPoissonHazardClosureTM, or SetCTMCClosureTM first.";

GetYDistribution[{pyr_?NumericQ, pbr_?NumericQ}] := Block[{PYR, PBR},
  Module[{TM\[Prime], vec},
    If[! MemberQ[{"MeanFieldClosure", "PoissonHazardClosure", "CTMCClosure"}, $ActiveTM],
      Message[GetYDistribution::activetm, $ActiveTM];
      Return[$Failed]
    ];

    If[$ActiveTM === "CTMCClosure",
      If[! ValueQ[ctmcClosureOptions], ctmcClosureOptions = {}];
      BuildCTMCClosureTM[{pyr, pbr}, Sequence @@ ctmcClosureOptions];
      vec = CTMCClosureStationaryWeights[
        Lookup[ctmcClosureRows, "Up"],
        Lookup[ctmcClosureRows, "Down"]
      ];
      Return[distY = EmpiricalDistribution[vec -> Values[as]]]
    ];

    TM\[Prime] = Chop@N@SparseArrayReplace[TM, {PYR -> pyr, PBR -> pbr}];
    vec = Chop@Eigenvectors[TM\[Prime], -1, Method -> "Arnoldi"][[1]];
    vec /= Total[vec];
    distY = EmpiricalDistribution[vec -> Values[as]]
  ]
];

GetYDistribution[pyr_?NumericQ] := GetYDistribution[{pyr, n PR - pyr}];


(* Reconstruct the fast pollinator cloud around the reduced Y process.
   The plotted distribution is truncated at Pmax, so Pmax should be large
   enough that the omitted Poisson tail is negligible. *)
PoissonizeP[distY_DataDistribution, {pyr_?NumericQ, pbr_?NumericQ}] := Module[
  {weights, states},
  weights = Flatten@Table[
    PDF[distY, y] PDF[PoissonDistribution[N@PMeanGivenY[y, {pyr, pbr}]], p],
    {y, 0, n},
    {p, 0, Pmax}
  ];
  states = Flatten[Table[{y, p}, {y, 0, n}, {p, 0, Pmax}], 1];
  EmpiricalDistribution[weights -> states]
];

PoissonizeP[distY_DataDistribution, pyr_?NumericQ] :=
  PoissonizeP[distY, {pyr, n PR - pyr}];


GetPDF[{pyr_?NumericQ, pbr_?NumericQ}] := PoissonizeP[GetYDistribution[{pyr, pbr}], {pyr, pbr}];

GetPDF[pyr_?NumericQ] := GetPDF[{pyr, n PR - pyr}];


InOutPY[pyr_?NumericQ] := Module[{localDistY = GetYDistribution[pyr]},
  If[localDistY === $Failed, Return[$Failed]];
  Expectation[y PMeanGivenY[y, pyr], y \[Distributed] localDistY]
];


(* invasion criteria *)
InvY:=InOutPY[\[Epsilon]]/\[Epsilon];
InvB:=(n PR-InOutPY[n PR-\[Epsilon]])/\[Epsilon]


(* ::Section:: *)
(*Full Stochastic Model*)


ClearAll[
  SetFullTM, SetFullMonocultureTMs,
  FullGetPDF, FullGetYPDF, FullGetBPDF,
  FullInOutYB, FullInOutY, FullInOutB,
  fullMonocultureStateCount
];

FullGetPDF::activetm =
  "The active transition matrix is `1`. Run SetFullTM before using full stochastic-model functions.";

FullGetYPDF::activetm =
  "Full monoculture transition matrices are unavailable. Run SetFullTM first.";

FullGetBPDF::activetm =
  "Full monoculture transition matrices are unavailable. Run SetFullTM first.";

fullMonocultureStateCount[] := (n + 1) (Pmax + 1);

SetFullTM := Block[{P\[EmptySet]R, PYR, PBR},
  Clear[TM, as, proc, TMY, TMB, asY, asB];

  as = Association@Flatten@Table[
    1 + Y + B (n + 3/2) - B^2/2 + P (n + 2) ((n + 1)/2) -> {Y, B, P},
    {P, 0, Pmax},
    {B, 0, n},
    {Y, 0, n - B}
  ];

  ind = Block[{y, b, p, y\[Prime], b\[Prime], p\[Prime]},
    With[{code = Round@{
        1 + y\[Prime] + b\[Prime] (n + 3/2) - b\[Prime]^2/2 + p\[Prime] (n + 2) ((n + 1)/2),
        1 + y + b (n + 3/2) - b^2/2 + p (n + 2) ((n + 1)/2)}
      },
      Compile[{{y, _Integer}, {b, _Integer}, {p, _Integer}, {y\[Prime], _Integer}, {b\[Prime], _Integer}, {p\[Prime], _Integer}},
        code,
        CompilationTarget -> "WVM",
        RuntimeAttributes -> Listable,
        Parallelization -> True
      ]
    ]
  ];

  TM = SparseArray[
    Flatten@Join[
      Table[ind[Y, B, P, Y, B, P - 1] -> fF["P-"], {P, 1, Pmax}, {B, 0, n}, {Y, 0, n - B}],
      Table[ind[Y, B, P, Y + 1, B, P] -> fF["Y+"], {P, 0, Pmax}, {B, 0, n}, {Y, 0, n - B - 1}],
      Table[ind[Y, B, P, Y, B + 1, P] -> fF["B+"], {P, 0, Pmax}, {Y, 0, n}, {B, 0, n - Y - 1}],
      Table[ind[Y, B, P, Y, B, P + 1] -> fF["P+"], {P, 0, Pmax - 1}, {B, 0, n}, {Y, 0, n - B}],
      Table[ind[Y, B, P, Y + 1, B, P + 1] -> fF["PY+"], {P, 0, Pmax - 1}, {B, 0, n}, {Y, 0, n - B - 1}],
      Table[ind[Y, B, P, Y, B + 1, P + 1] -> fF["PB+"], {P, 0, Pmax - 1}, {Y, 0, n}, {B, 0, n - Y - 1}],
      Table[ind[Y, B, P, Y - 1, B, P] -> fF["Y-"], {P, 0, Pmax}, {B, 0, n}, {Y, 1, n - B}],
      Table[ind[Y, B, P, Y, B - 1, P] -> fF["B-"], {P, 0, Pmax}, {B, 1, n}, {Y, 0, n - B}],
      Table[ind[Y, B, P, Y, B, P] -> -Total[Values[fF]], {P, 0, Pmax - 1}, {B, 0, n}, {Y, 0, n - B}],
      Table[ind[Y, B, P, Y, B, P] -> -Total[Lookup[fF, {"P-", "Y+", "B+", "Y-", "B-"}]],
        {P, Pmax, Pmax}, {B, 0, n}, {Y, 0, n - B}]
    ],
    {(n + 2) (n + 1) ((Pmax + 1)/2), (n + 2) (n + 1) ((Pmax + 1)/2)},
    0.
  ];

  proc = {
    {{P -> P - 1}, fF["P-"]},
    {{Y -> Y + 1}, fF["Y+"]},
    {{B -> B + 1}, fF["B+"]},
    {{P -> P + 1}, fF["P+"]},
    {{P -> P + 1, Y -> Y + 1}, fF["PY+"]},
    {{P -> P + 1, B -> B + 1}, fF["PB+"]},
    {{Y -> Y - 1}, fF["Y-"]},
    {{B -> B - 1}, fF["B-"]}
  };

  SetFullMonocultureTMs[];
  $ActiveTM = "Full";
];


SetFullMonocultureTMs[] := Module[
  {
    indYP, indBP,
    fFY, fFB
  },

  asY = Association@Flatten@Table[
    1 + Y + P (n + 1) -> {Y, P},
    {P, 0, Pmax},
    {Y, 0, n}
  ];

  asB = Association@Flatten@Table[
    1 + B + P (n + 1) -> {B, P},
    {P, 0, Pmax},
    {B, 0, n}
  ];

  indYP = Block[
    {y, p, y\[Prime], p\[Prime]},
    With[
      {code = Round@{1 + y\[Prime] + p\[Prime] (n + 1), 1 + y + p (n + 1)}},
      Compile[
        {{y, _Integer}, {p, _Integer}, {y\[Prime], _Integer}, {p\[Prime], _Integer}},
        code,
        CompilationTarget -> "WVM",
        RuntimeAttributes -> Listable,
        Parallelization -> True
      ]
    ]
  ];

  indBP = Block[
    {b, p, b\[Prime], p\[Prime]},
    With[
      {code = Round@{1 + b\[Prime] + p\[Prime] (n + 1), 1 + b + p (n + 1)}},
      Compile[
        {{b, _Integer}, {p, _Integer}, {b\[Prime], _Integer}, {p\[Prime], _Integer}},
        code,
        CompilationTarget -> "WVM",
        RuntimeAttributes -> Listable,
        Parallelization -> True
      ]
    ]
  ];

  (* Y monoculture: B is absent, and there is no regional B pool. *)
  fFY = fF /. {B -> 0, PBR -> 0};

  TMY = SparseArray[
    Flatten@Join[
      Table[indYP[Y, P, Y, P - 1] -> fFY["P-"], {P, 1, Pmax}, {Y, 0, n}],
      Table[indYP[Y, P, Y + 1, P] -> fFY["Y+"], {P, 0, Pmax}, {Y, 0, n - 1}],
      Table[indYP[Y, P, Y, P + 1] -> fFY["P+"], {P, 0, Pmax - 1}, {Y, 0, n}],
      Table[indYP[Y, P, Y + 1, P + 1] -> fFY["PY+"], {P, 0, Pmax - 1}, {Y, 0, n - 1}],
      Table[indYP[Y, P, Y - 1, P] -> fFY["Y-"], {P, 0, Pmax}, {Y, 1, n}],
      Table[indYP[Y, P, Y, P] -> -Total[Lookup[fFY, {"P-", "Y+", "P+", "PY+", "Y-"}]],
        {P, 0, Pmax - 1}, {Y, 0, n}],
      Table[indYP[Y, P, Y, P] -> -Total[Lookup[fFY, {"P-", "Y+", "Y-"}]],
        {P, Pmax, Pmax}, {Y, 0, n}]
    ],
    {fullMonocultureStateCount[], fullMonocultureStateCount[]},
    0.
  ];

  (* B monoculture: Y is absent, and there is no regional Y pool. *)
  fFB = fF /. {Y -> 0, PYR -> 0};

  TMB = SparseArray[
    Flatten@Join[
      Table[indBP[B, P, B, P - 1] -> fFB["P-"], {P, 1, Pmax}, {B, 0, n}],
      Table[indBP[B, P, B + 1, P] -> fFB["B+"], {P, 0, Pmax}, {B, 0, n - 1}],
      Table[indBP[B, P, B, P + 1] -> fFB["P+"], {P, 0, Pmax - 1}, {B, 0, n}],
      Table[indBP[B, P, B + 1, P + 1] -> fFB["PB+"], {P, 0, Pmax - 1}, {B, 0, n - 1}],
      Table[indBP[B, P, B - 1, P] -> fFB["B-"], {P, 0, Pmax}, {B, 1, n}],
      Table[indBP[B, P, B, P] -> -Total[Lookup[fFB, {"P-", "B+", "P+", "PB+", "B-"}]],
        {P, 0, Pmax - 1}, {B, 0, n}],
      Table[indBP[B, P, B, P] -> -Total[Lookup[fFB, {"P-", "B+", "B-"}]],
        {P, Pmax, Pmax}, {B, 0, n}]
    ],
    {fullMonocultureStateCount[], fullMonocultureStateCount[]},
    0.
  ];
];


FullGetPDF[{p\[EmptySet]r_?NumericQ, pyr_?NumericQ, pbr_?NumericQ}] :=
  Block[{P\[EmptySet]R, PYR, PBR},
    Module[{TM\[Prime], vec},
      If[$ActiveTM =!= "Full",
        Message[FullGetPDF::activetm, $ActiveTM];
        Return[$Failed]
      ];

      TM\[Prime] = SparseArrayReplace[TM, {PYR -> pyr, PBR -> pbr, P\[EmptySet]R -> p\[EmptySet]r}];
      vec = Chop@Eigenvectors[TM\[Prime], -1, Method -> "Arnoldi"][[1]];
      vec /= Total[vec];
      fullDist = EmpiricalDistribution[vec -> Values[as]]
    ]
  ];

FullGetPDF[{pyr_?NumericQ, pbr_?NumericQ}] := FullGetPDF[{n PR - pyr - pbr, pyr, pbr}];


FullGetYPDF[pyr_?NumericQ] :=
  Block[{P\[EmptySet]R, PYR, PBR},
    Module[{TMY\[Prime], vec},
      If[! ValueQ[TMY] || Dimensions[TMY] =!= {fullMonocultureStateCount[], fullMonocultureStateCount[]},
        Message[FullGetYPDF::activetm];
        Return[$Failed]
      ];

      TMY\[Prime] = SparseArrayReplace[TMY, {PYR -> pyr, PBR -> 0, P\[EmptySet]R -> n PR - pyr}];
      vec = Chop@Eigenvectors[TMY\[Prime], -1, Method -> "Arnoldi"][[1]];
      vec /= Total[vec];
      fullDistY = EmpiricalDistribution[vec -> Values[asY]]
    ]
  ];

FullGetBPDF[pbr_?NumericQ] :=
  Block[{P\[EmptySet]R, PYR, PBR},
    Module[{TMB\[Prime], vec},
      If[! ValueQ[TMB] || Dimensions[TMB] =!= {fullMonocultureStateCount[], fullMonocultureStateCount[]},
        Message[FullGetBPDF::activetm];
        Return[$Failed]
      ];

      TMB\[Prime] = SparseArrayReplace[TMB, {PYR -> 0, PBR -> pbr, P\[EmptySet]R -> n PR - pbr}];
      vec = Chop@Eigenvectors[TMB\[Prime], -1, Method -> "Arnoldi"][[1]];
      vec /= Total[vec];
      fullDistB = EmpiricalDistribution[vec -> Values[asB]]
    ]
  ];

FullInOutYB[{pyr_?NumericQ, pbr_?NumericQ}] := Expectation[{p y, p b}, {y, b, p} \[Distributed] FullGetPDF[{pyr, pbr}]];

FullInOutY[pyr_?NumericQ] := Expectation[p y, {y, p} \[Distributed] FullGetYPDF[pyr]];

FullInOutB[pbr_?NumericQ] := Expectation[p b, {b, p} \[Distributed] FullGetBPDF[pbr]];


(* invasion criteria *)
FullInvY:=Module[{eqB},
	eqB=FindRoot[FullInOutB[pbr]==pbr,{pbr,0.99 PR n}];
	FullInOutYB[{\[Epsilon],pbr/.eqB}][[1]]/\[Epsilon]
];

FullInvB:=Module[{eqY},
	eqY=FindRoot[FullInOutY[pyr]==pyr,{pyr,0.99 PR n}];
	FullInOutYB[{pyr/.eqY,\[Epsilon]}][[2]]/\[Epsilon]
];


(* ::Section:: *)
(*Distribution Plotting and Decomposition*)


ClearAll[DecomposeDistribution];
DecomposeDistribution[dist_DataDistribution] := Module[
  {
    domain, xmin, ymin, nx, ny, dat, idat, markers, wsc,
    ncomps, comps, weights, indices
  },
  Which[
    dist["Dimension"] === 2,
    domain = dist["Domain"];
    {xmin, ymin} = domain[[{1, 2}, 1]];
    nx = Length[domain[[1]]];
    ny = Length[domain[[2]]];
    dat = SparseArray[Map[Plus[#, {1 - xmin, 1 - ymin}] &, dist[[2, 2]]] -> dist[[2, 1]]];
    idat = Image[-dat];
    markers = MinDetect[idat];
    wsc = WatershedComponents[idat, markers, Method -> "Basins"];
    ncomps = Length[Union[Flatten[wsc]]];
    comps = Table[
      Table[If[wsc[[i, j]] == comp, dat[[i, j]], 0.], {i, nx}, {j, ny}],
      {comp, ncomps, 1, -1}
    ];
    weights = Total[comps, {2, 4}];
    indices = Flatten[Table[{x, y}, {x, domain[[1]]}, {y, domain[[2]]}], 1];
    Return[{weights, Table[EmpiricalDistribution[Flatten[comp] -> indices], {comp, comps}]}],

    dist["Dimension"] === 3,
    dat = SparseArray[Map[Plus[#, {1, 1, 1}] &, dist[[2, 2]]] -> dist[[2, 1]]];
    idat = Image3D[-dat];
    markers = MinDetect[idat];
    wsc = WatershedComponents[idat, markers, Method -> "Basins"];
    ncomps = Length[Union[Flatten[wsc]]];
    comps = Table[
      Table[
        If[wsc[[y + 1, b + 1, p + 1]] == comp, dat[[y + 1, b + 1, p + 1]], 0.],
        {y, 0, n},
        {b, 0, n},
        {p, 0, Pmax}
      ],
      {comp, ncomps, 1, -1}
    ];
    weights = Total[comps, {2, 4}];
    indices = Flatten[Table[{y, b, p}, {y, 0, n}, {b, 0, n}, {p, 0, Pmax}], 2];
    Return[{weights, Table[EmpiricalDistribution[Flatten[comp] -> indices], {comp, comps}]}],

    True,
    Print["only 2D and 3D distributions currently supported"]
  ]
];


ClearAll[PlotDistribution, PlotDistribution1D, PlotDistribution2D, PlotDistribution3D];

PlotDistribution[dist_, opts___?OptionQ] :=
  Which[
    dist["Dimension"] === 1,
    PlotDistribution1D[dist, opts],
    dist["Dimension"] === 2,
    PlotDistribution2D[dist, opts],
    dist["Dimension"] === 3,
    PlotDistribution3D[dist, opts],
    True,
    Print["only 1D, 2D, and 3D distributions currently supported"]
  ];

PlotDistribution1D[dist_, opts___?OptionQ] :=
  ListPlot[
    Transpose[{dist["Domain"], dist["Weights"]}],
    Evaluate[Sequence @@ {opts}],
    PlotRange -> All
  ];

PlotDistribution2D[dist_, opts___?OptionQ] :=
  ArrayPlot[
    Reverse@Chop@Normal[SparseArray[Map[Plus[#, {1, 1}] &, dist[[2, 2]]] -> dist[[2, 1]]]],
    DataRange -> {{0, Pmax}, {0, n}},
    Evaluate[Sequence @@ {opts}],
    AspectRatio -> 1,
    FrameLabel -> {{"P", None}, {"Y", "B"}},
    FrameTicks -> {
      {Table[{y, y}, {y, 0, n, 10}], Table[{y, n - y}, {y, 0, n, 10}]},
      {Table[p, {p, 0, Pmax, 1}], (*Table[p, {p, 0, Pmax, 1}]*)None}
    },
    Mesh -> False,
    PlotRange -> {{-0.5, Pmax + 0.5}, {-0.5, n + 0.5}, All},
    PlotRangePadding -> Scaled[0.02]
  ];

PlotDistribution3D[dist_, opts___?OptionQ] :=
  ListDensityPlot3D[
    SparseArray[Map[Plus[#, {1, 1, 1}] &, dist[[2, 2]]] -> dist[[2, 1]]],
    BoxRatios -> {0.4, 1, 1},
    ViewPoint -> {1, -2, -2},
    ViewVertical -> {1, 0, 0},
    AxesLabel -> {"P", "B", "Y"},
    DataRange -> {{0, Pmax}, {0, n}, {0, n}},
    Evaluate[Sequence @@ {opts}]
  ];
