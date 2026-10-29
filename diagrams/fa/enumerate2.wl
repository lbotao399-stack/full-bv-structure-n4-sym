(* enumerate2.wl -- round 2: FeynArts enumeration + Paint (TeX) of
     (T1) one-loop 1PI graphs with P(eps,infty) edges only, cubic and quartic external sets;
     (T2) QME-defect graphs with external c and A|_V:
          - one special edge K_eps - K_0 joining two distinct vertices (from (1/2)({I,I}_eps - {I,I}_0)),
            all other edges P(eps,infty), all one-loop connected graphs (1PI and reducible);
          - one K_eps self-loop (from hbar Delta_eps I), all other edges P(eps,infty).
   Drawing-only model N4Draw2 (model/N4Draw.gen, model/N4Draw2.mod).
   Run:  WolframKernel -noinit -noprompt -run 'CheckAbort[Check[Get[".../enumerate2.wl"],Quit[2]],Quit[2]];Quit[0]'
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
  "model_mod_sha256" -> sha[root <> "/model/N4Draw2.mod"],
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
init = FeynArts`InitializeModel["N4Draw2", FeynArts`GenericModel -> "N4Draw", FeynArts`Reinitialize -> True];
Print["model initialized: ", init, "  classes: ", Length[FeynArts`M$ClassesDescription],
  "  couplings: ", Length[FeynArts`M$CouplingMatrices]];
If[init =!= True, Quit[2]];
gate["model_classes"] = Length[FeynArts`M$ClassesDescription];
gate["model_couplings"] = Length[FeynArts`M$CouplingMatrices];

pList = {S[1], S[2], S[3], S[4]};
kList = {S[21], S[22], S[23], S[24]};
dList = {S[51], S[52], S[53], S[54]};
extOnly = {S[41], S[42], S[43]};   (* external antifield legs; never internal *)

base[f_] := f /. -x_ :> x;
isInt[pr_] := MatchQ[pr, Propagator[Internal | _FALoop][__]];
countSp[props_, lst_] := Count[props, pr_ /; isInt[pr] && MemberQ[lst, base[pr[[3]]]]];
selfLoopSp[props_, lst_] := Count[props, pr_ /; isInt[pr] && MemberQ[lst, base[pr[[3]]]] && pr[[1]] === pr[[2]]];

toS[x_] := ToString[x, InputForm];

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
  ins = InsertFields[tops, proc, InsertionLevel -> {Classes}, Model -> "N4Draw2",
    GenericModel -> "N4Draw", ExcludeParticles -> Join[excl, extOnly]];
  ds = If[ins === {} || Length[ins] === 0, {}, flat[ins]];
  keep = Select[ds, sel[diagProps[#]] &];
  AppendTo[summary, <|"job" -> job, "process" -> toS[proc], "topologies" -> Length[tops],
    "inserted" -> Length[ds], "kept" -> Length[keep]|>];
  Print[job, "  ", toS[proc], "  topologies=", Length[tops], "  inserted=", Length[ds], "  kept=", Length[keep]];
  Do[
    k++;
    AppendTo[allRecords, record[d, job, k]];
    If[paint, paintOne[d, root <> "/tex/" <> job <> "_" <> IntegerString[k, 10, 2] <> ".tex"]],
    {d, keep}];
];

pOnlySel = (countSp[#, kList] == 0 && countSp[#, dList] == 0) &;

(* ---------- T1 cubic: one-loop 1PI, P edges only ---------- *)
t12 = CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Reducible];
t21 = CreateTopologies[1, 2 -> 1, ExcludeTopologies -> Reducible];
cubic = {
  {"c3_AVAVAV",   t12, {S[1]} -> {S[1], S[1]}},
  {"c3_AVAWAW",   t12, {S[1]} -> {S[2], S[2]}},
  {"c3_lamlamAV", t12, {S[3]} -> {S[3], S[1]}},
  {"c3_lamlamAW", t12, {S[3]} -> {S[3], S[2]}},
  {"c3_cAVAsV",   t12, {-S[4]} -> {S[1], -S[4]}},
  {"c3_cAWAsW",   t12, {-S[4]} -> {S[2], -S[41]}},
  {"c3_clamlams", t12, {-S[4]} -> {S[3], -S[42]}},
  {"c3_cccs",     t21, {-S[4], -S[4]} -> {-S[43]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[kList, dList], pOnlySel], {s, cubic}];

(* ---------- T1 quartic: one-loop 1PI, P edges only ---------- *)
t22 = CreateTopologies[1, 2 -> 2, ExcludeTopologies -> Reducible];
quartic = {
  {"q4_AVAVAVAV", t22, {S[1], S[1]} -> {S[1], S[1]}},
  {"q4_AVAVAWAW", t22, {S[1], S[1]} -> {S[2], S[2]}},
  {"q4_AWAWAWAW", t22, {S[2], S[2]} -> {S[2], S[2]}}
};
Do[runJob[s[[1]], s[[2]], s[[3]], Join[kList, dList], pOnlySel], {s, quartic}];

(* ---------- T2: QME defect, external c and A|_V ---------- *)
t11all = CreateTopologies[1, 1 -> 1];
Print["T2 topologies (1->1, all): ", Length[t11all]];
(* (1/2)({I,I}_eps - {I,I}_0): exactly one K_eps - K_0 edge between two distinct vertices *)
bracketSel = (countSp[#, dList] == 1 && selfLoopSp[#, dList] == 0 && countSp[#, kList] == 0) &;
runJob["qd_cAV", t11all, {-S[4]} -> {S[1]}, kList, bracketSel];
(* hbar Delta_eps I: exactly one K_eps edge, a self-loop at one vertex *)
deltaSel = (countSp[#, kList] == 1 && selfLoopSp[#, kList] == 1 && countSp[#, dList] == 0) &;
runJob["qdD_cAV", t11all, {-S[4]} -> {S[1]}, dList, deltaSel];

Export[root <> "/data/fa_records_r2.json", allRecords, "RawJSON"];
Export[root <> "/data/fa_summary_r2.json", <|"gate" -> gate, "jobs" -> summary|>, "RawJSON"];
Print["DONE records=", Length[allRecords]];
