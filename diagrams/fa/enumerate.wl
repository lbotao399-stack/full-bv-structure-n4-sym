(* enumerate.wl -- FeynArts enumeration + Paint (TeX) for the N4_Costello_BV diagrams.
   Drawing-only model N4Draw (model/N4Draw.gen, model/N4Draw.mod).
   Run:  WolframKernel -noinit -noprompt -run 'CheckAbort[Check[Get[".../enumerate.wl"],Quit[2]],Quit[2]];Quit[0]'
*)
ClearAll["Global`*"];
$HistoryLength = 0;
$FeynCalcStartupMessages = False;
$LoadAddOns = {"FeynArts"};
Needs["FeynCalc`"];

root = "/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/diagrams";
faDir = "/Users/libotao/Library/Wolfram/Applications/FeynCalc/FeynArts/";
If[FeynArts`$FeynArtsDir =!= faDir, Print["FAIL: wrong FeynArts copy"]; Quit[2]];

sha[p_] := IntegerString[FileHash[p, "SHA256"], 16, 64];
gate = <|
  "kernel" -> $Version,
  "system" -> $SystemID,
  "feyncalc_version" -> FeynCalc`$FeynCalcVersion,
  "feyncalc_dir" -> FeynCalc`$FeynCalcDirectory,
  "feyncalc_sha256" -> sha[FileNameJoin[{FeynCalc`$FeynCalcDirectory, "FeynCalc.m"}]],
  "feynarts_version" -> FeynArts`$FeynArtsVersion,
  "feynarts_dir" -> FeynArts`$FeynArtsDir,
  "feynarts_sha256" -> sha[FileNameJoin[{FeynArts`$FeynArtsDir, "FeynArts.m"}]],
  "model_mod_sha256" -> sha[root <> "/model/N4Draw.mod"],
  "model_gen_sha256" -> sha[root <> "/model/N4Draw.gen"]
|>;
checks = <|
  "kernel" -> ($Version == "15.0.0 for Mac OS X ARM (64-bit) (May 19, 2026)"),
  "feyncalc-version" -> (FeynCalc`$FeynCalcVersion == "10.2.1"),
  "feyncalc-hash" -> (gate["feyncalc_sha256"] == "4238350af2ff91a6687ea937446b9e6077318914593afe963d4d10026dfc5165"),
  "feynarts-version" -> (FeynArts`$FeynArtsVersion == "FeynArts 3.12 (27 Mar 2025)"),
  "feynarts-hash" -> (gate["feynarts_sha256"] == "7b8a693596aff5279b0f221cec264594241adb34418c848983dd18cbfc11437a")
|>;
Print["GATE ", InputForm[checks]];
If[!And @@ Values[checks], Print["FAIL: runtime gate"]; Quit[2]];

FeynArts`$ModelPath = {root <> "/model", faDir <> "Models"};
FeynArts`$ShapeDataDir = root <> "/fa/ShapeData";
FeynArts`$FAVerbose = 0;
init = FeynArts`InitializeModel["N4Draw", FeynArts`GenericModel -> "N4Draw", FeynArts`Reinitialize -> True];
Print["model initialized: ", init, "  classes: ", Length[FeynArts`M$ClassesDescription],
  "  couplings: ", Length[FeynArts`M$CouplingMatrices]];
If[init =!= True, Quit[2]];
gate["model_classes"] = Length[FeynArts`M$ClassesDescription];
gate["model_couplings"] = Length[FeynArts`M$CouplingMatrices];

pList = {S[1], S[2], S[3], S[4]};
pinList = {S[11], S[12], S[13], S[14]};
kList = {S[21], S[22], S[23], S[24]};
srcList = {S[31], S[32]};  (* sources never run on internal lines *)

base[f_] := f /. -x_ :> x;
isInt[pr_] := MatchQ[pr, Propagator[Internal | _FALoop][__]];
countSp[props_, lst_] := Count[props, pr_ /; isInt[pr] && MemberQ[lst, base[pr[[3]]]]];

toS[x_] := ToString[x, InputForm];

(* flatten a classes-level TopologyList into {head, top, genGraph, clsGraph} *)
flat[ins_] := Module[{out = {}},
  Scan[Function[tg,
      Scan[Function[gg,
          Scan[Function[cg, AppendTo[out, {Head[ins], tg[[1]], gg[[1]], cg}]],
            List @@ gg[[2]]]],
        List @@ tg[[2]]]],
    List @@ ins];
  out];

diagProps[{h_, top_, gg_, cg_}] := List @@ (top /. (List @@ cg) /. (List @@ gg));

vstr[Vertex[e_, c_][n_]] := "V" <> ToString[e] <> "c" <> ToString[c] <> "_" <> ToString[n];
vstr[Vertex[e_][n_]] := "V" <> ToString[e] <> "_" <> ToString[n];

record[d_, job_, idx_] := Module[{props = diagProps[d], cg = d[[4]], fps},
  fps = FieldPoints[d[[4]], d[[2]]];
  <|"job" -> job, "index" -> idx,
    "symmetry_factor" -> toS[cg[[0, 1]]],
    "propagators" -> Map[<|"type" -> toS[#[[0, 1]]], "from" -> vstr[#[[1]]], "to" -> vstr[#[[2]]],
        "field" -> toS[#[[3]]]|> &, props],
    "fieldpoints" -> Map[toS, fps]|>];

single[{h_, top_, gg_, cg_}] := h[top -> Insertions[Generic][gg -> Insertions[Classes][cg]]];

paintOne[d_, file_] := Paint[single[d], PaintLevel -> {Classes}, ColumnsXRows -> 1,
  SheetHeader -> None, Numbering -> None, AutoEdit -> False, FieldNumbers -> False,
  ImageSize -> {200, 200}, DisplayFunction -> (Export[file, #] &)];

allRecords = {};
summary = {};

runJob[job_String, tops_, proc_, excl_, sel_, paint_:True] := Module[{ins, ds, keep, k = 0},
  ins = InsertFields[tops, proc, InsertionLevel -> {Classes}, Model -> "N4Draw",
    GenericModel -> "N4Draw", ExcludeParticles -> Join[excl, srcList]];
  ds = If[ins === {} || Length[ins] === 0, {}, flat[ins]];
  keep = Select[ds, sel[diagProps[#]] &];
  AppendTo[summary, <|"job" -> job, "process" -> toS[proc], "inserted" -> Length[ds],
    "kept" -> Length[keep]|>];
  Print[job, "  ", toS[proc], "  inserted=", Length[ds], "  kept=", Length[keep]];
  Do[
    k++;
    AppendTo[allRecords, record[d, job, k]];
    If[paint, paintOne[d, root <> "/tex/" <> job <> "_" <> IntegerString[k, 10, 2] <> ".tex"]],
    {d, keep}];
];

noneSel = (countSp[#, pinList] == 0 && countSp[#, kList] == 0) &;
onePinSel = (countSp[#, pinList] == 1 && countSp[#, kList] == 0) &;
oneKSel = (countSp[#, kList] == 1 && countSp[#, pinList] == 0) &;

(* ---------- (A) one-loop 1PI two-point graphs, every degree-0 even sector ---------- *)
t2 = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Reducible];
sectors2 = {
  {"tp_AVAV", {S[1]} -> {S[1]}},
  {"tp_AVAW", {S[1]} -> {S[2]}},
  {"tp_AWAW", {S[2]} -> {S[2]}},
  {"tp_lamlam", {S[3]} -> {S[3]}},
  {"tp_cAsV", {-S[4]} -> {-S[4]}},
  {"tp_cAsW", {-S[4]} -> {-S[22]}}
};
Do[runJob[s[[1]], t2, s[[2]], Join[pinList, kList], noneSel], {s, sectors2}];

(* ---------- (B) route (i): eps d/deps boundary versions (one L_eps edge) ---------- *)
Do[runJob["bi_" <> StringDrop[s[[1]], 3], t2, s[[2]], kList, onePinSel], {s, sectors2}];

(* ---------- (C) one-loop 1PI graphs with one O insertion ---------- *)
tO0 = CreateTopologies[1, 1 -> 0, ExcludeTopologies -> Reducible];
tO1 = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Reducible];
tO2 = CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Reducible];
sectorsO = {
  {"ins_vac", tO0, {S[31]} -> {}},
  {"ins_AV", tO1, {S[31]} -> {S[1]}},
  {"ins_AW", tO1, {S[31]} -> {S[2]}},
  {"ins_AVAV", tO2, {S[31]} -> {S[1], S[1]}},
  {"ins_AVAW", tO2, {S[31]} -> {S[1], S[2]}},
  {"ins_AWAW", tO2, {S[31]} -> {S[2], S[2]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[pinList, kList], noneSel], {s, sectorsO}];
(* omega < 0 sectors: enumerated, not painted *)
sectorsOneg = {
  {"insneg_lamlam", tO2, {S[31]} -> {S[3], S[3]}},
  {"insneg_cAsV", tO2, {S[31]} -> {-S[4], -S[4]}},
  {"insneg_cAsW", tO2, {S[31]} -> {-S[4], -S[22]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[pinList, kList], noneSel, False], {s, sectorsOneg}];

(* ---------- (C') one-particle-reducible one-loop graphs with one O insertion ---------- *)
rO0 = CreateTopologies[1, 1 -> 0, ExcludeTopologies -> Irreducible];
rO1 = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Irreducible];
rO2 = CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Irreducible];
sectorsR = {
  {"red_vac", rO0, {S[31]} -> {}},
  {"red_AV", rO1, {S[31]} -> {S[1]}},
  {"red_AW", rO1, {S[31]} -> {S[2]}},
  {"red_AVAV", rO2, {S[31]} -> {S[1], S[1]}},
  {"red_AVAW", rO2, {S[31]} -> {S[1], S[2]}},
  {"red_AWAW", rO2, {S[31]} -> {S[2], S[2]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[pinList, kList], noneSel], {s, sectorsR}];

(* ---------- (D) route (i) boundary versions of the O-insertion graphs ---------- *)
Do[runJob["biO_" <> StringDrop[s[[1]], 4], s[[2]], s[[3]], kList, onePinSel], {s, sectorsO}];

(* ---------- (E) route (ii): 1PI one-loop graphs with O, S^(1)(Delta) and exactly one K_eps edge ---------- *)
tD0 = CreateTopologies[1, 2 -> 0, ExcludeTopologies -> Reducible];
tD1 = CreateTopologies[1, 2 -> 1, ExcludeTopologies -> Reducible];
tD2 = CreateTopologies[1, 2 -> 2, ExcludeTopologies -> Reducible];
sectorsD = {
  {"bii_vac", tD0, {S[31], S[32]} -> {}},
  {"bii_AV", tD1, {S[31], S[32]} -> {S[1]}},
  {"bii_AW", tD1, {S[31], S[32]} -> {S[2]}},
  {"bii_AVAV", tD2, {S[31], S[32]} -> {S[1], S[1]}},
  {"bii_AVAW", tD2, {S[31], S[32]} -> {S[1], S[2]}},
  {"bii_AWAW", tD2, {S[31], S[32]} -> {S[2], S[2]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], pinList, oneKSel], {s, sectorsD}];

(* ---------- (F) route (ii): tree graphs with the O counterterm vertex and one K_eps edge ---------- *)
ctD2 = CreateCTTopologies[1, 2 -> 2];
ctSel = (countSp[#, kList] == 1 && countSp[#, pinList] == 0 && countSp[#, pList] == 0) &;
sectorsCT = {
  {"biict_AVAV", ctD2, {S[31], S[32]} -> {S[1], S[1]}},
  {"biict_AVAW", ctD2, {S[31], S[32]} -> {S[1], S[2]}},
  {"biict_AWAW", ctD2, {S[31], S[32]} -> {S[2], S[2]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[pList, pinList], ctSel], {s, sectorsCT}];

Export[root <> "/data/fa_records.json", allRecords, "RawJSON"];
Export[root <> "/data/fa_summary.json", <|"gate" -> gate, "jobs" -> summary|>, "RawJSON"];
Print["DONE records=", Length[allRecords]];
