"""Merge parts/*/*.json into dr_survey.json and dr_survey.txt.

Records describing the same arXiv work (written by two section workers) are merged into one
record: the first in section order is the primary record, the other is kept verbatim under
"alias_records" (no content dropped).
"""
import json, glob, os, re, sys, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rec_io import rec_txt, KEYS

ROOT = "/Users/libotao/N4_BV_DR/literature"
SECTION_TITLES = {
    "1a": "1a  N=4 (and deformed) spectrum from Feynman graphs / supergraphs in DR-type schemes",
    "1b": "1b  Integrability-based high-loop spectrum and the identity of the '9 loops' claims",
    "1c": "1c  Correlator integrands (Lagrangian insertion, hidden symmetry, graphical bootstrap)",
    "1d": "1d  Amplitudes and UV divergences (FDH/DRED-type continuation, critical dimensions); dual conformal Ward identity",
    "1e": "1e  N=1 beta functions in DRED at 3-4 loops (and manifestly supersymmetric contrasts)",
    "2":  "2   DRED consistency, FDH/HV, epsilon-scalars and evanescent couplings",
    "3":  "3   BMHV ('t Hooft-Veltman/Breitenlohner-Maison) and symmetry-restoring counterterms; SUSY applications",
    "4":  "4   Algebraic renormalization of SUSY gauge theories with constant SUSY ghosts",
    "5":  "5   Standard Model at >= 3 loops and gamma5",
    "6":  "6   Supporting literature for the synthesis: BV/QAP in DR, conformal symmetry in d != 4, all-order protection",
}
ORDER = ["1a", "1b", "1c", "1d", "1e", "2", "3", "4", "5", "6"]
ID_RE = re.compile(r"^S(\d)([a-e]?)\.(\d+)_")
ARX_RE = re.compile(r"(\d{4}\.\d{4,5}|[a-z\-]+(?:\.[A-Z]{2})?/\d{7})")


def skey(rid):
    m = ID_RE.match(rid)
    if not m:
        raise ValueError(rid)
    sec = m.group(1) + m.group(2)
    return (ORDER.index(sec), int(m.group(3)), rid), sec


def arx(r):
    m = ARX_RE.search(r.get("arxiv", ""))
    return m.group(1) if m else None


def main():
    recs = []
    for f in glob.glob(os.path.join(ROOT, "parts", "*", "*.json")):
        with open(f) as fh:
            r = json.load(fh)
        missing = [k for k in KEYS if k not in r]
        if missing:
            raise ValueError(f"{f}: missing {missing}")
        recs.append(r)
    recs.sort(key=lambda r: skey(r["id"])[0])

    by_arx, merged, aliases = {}, [], {}
    for r in recs:
        a = arx(r)
        if a and a in by_arx:
            prim = by_arx[a]
            prim.setdefault("alias_records", []).append(r)
            aliases[r["id"]] = prim["id"]
            continue
        r = dict(r)
        r["arxiv_id_normalized"] = a or ""
        if a:
            by_arx[a] = r
        merged.append(r)

    for r in merged:
        r.setdefault("alias_records", [])
        r["also_recorded_as"] = [x["id"] for x in r["alias_records"]]

    with open(os.path.join(ROOT, "dr_survey.json"), "w") as fh:
        json.dump(merged, fh, ensure_ascii=False, indent=1)

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out = []
    out.append("DR-TYPE REGULARIZATION AND SYMMETRY AT HIGH LOOP ORDER -- LITERATURE SURVEY")
    out.append(f"Generated {now} from parts/*/*.json by tools/merge.py")
    out.append(f"Works: {len(merged)}   (records written: {len(recs)}; {len(aliases)} duplicate records merged as aliases)")
    out.append("")
    out.append("Field legend")
    out.append("  Claim           what the work claims (paraphrase)")
    out.append("  Scheme          regularization, subtraction, gamma5 / epsilon-scalar treatment")
    out.append("  Proved          established inside the work (within its stated assumptions)")
    out.append("  Assumed         used without proof")
    out.append("  Checked part.   verified only for some quantities / orders / numerically")
    out.append("  Evidence        section, equation, table numbers as seen in the source")
    out.append("  Access          what was actually read: abstract only / full text (ar5iv, arxiv html, PDF) / secondary source")
    out.append("  'eq. no. not verified' marks numbers that were not seen in the text.")
    out.append("")
    out.append("CONTENTS")
    cur = None
    for r in merged:
        sec = skey(r["id"])[1]
        if sec != cur:
            out.append("  " + SECTION_TITLES[sec])
            cur = sec
        al = f"  [also: {', '.join(r['also_recorded_as'])}]" if r["also_recorded_as"] else ""
        out.append(f"      {r['id']:<44s} {r['arxiv'][:40]}{al}")
    out.append("")
    cur = None
    for r in merged:
        sec = skey(r["id"])[1]
        if sec != cur:
            out.append("=" * 100)
            out.append(SECTION_TITLES[sec])
            out.append("=" * 100)
            out.append("")
            cur = sec
        out.append(rec_txt(r).rstrip("\n"))
        for a in r["alias_records"]:
            out.append(f"  ALSO RECORDED AS {a['id']} (section {a['section']}); that record, verbatim:")
            out.append("\n".join("    " + line for line in rec_txt(a).rstrip("\n").split("\n")))
        out.append("")
    out.append("=" * 100)
    out.append("INDEX BY arXiv ID")
    for r in sorted(merged, key=lambda r: (r["arxiv_id_normalized"] == "", r["arxiv_id_normalized"])):
        out.append(f"  {r['arxiv_id_normalized'] or '(no arXiv id)':<22s} {r['id']}")
    if aliases:
        out.append("")
        out.append("ALIASES (duplicate records merged into the primary record)")
        for k, v in sorted(aliases.items()):
            out.append(f"  {k:<44s} -> {v}")
    with open(os.path.join(ROOT, "dr_survey.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(len(recs), "records,", len(merged), "works,", len(aliases), "aliases")
    for k, v in sorted(aliases.items()):
        print("  alias", k, "->", v)


if __name__ == "__main__":
    main()
