"""Write a survey record as <ID>.json and <ID>.txt (layout of parts/RECORD_FORMAT.txt)."""
import json, os, datetime

KEYS = ["id", "section", "authors", "title", "year", "arxiv", "journal", "claim", "loop_order",
        "theory", "scheme", "symmetry_argument", "symmetry_status", "evidence", "red_flags",
        "cross_checks", "source_access", "confidence", "notes"]


def _lst(xs, ind="      - "):
    if not xs:
        return " -"
    return "\n" + "\n".join(ind + str(x) for x in xs)


def rec_txt(r):
    s = r["symmetry_status"]
    out = []
    out.append(f"[{r['id']}] {r['authors']} -- {r['title']} ({r['year']})")
    out.append(f"  arXiv: {r['arxiv']}   Journal: {r['journal']}")
    out.append(f"  Section: {r['section']}")
    out.append(f"  Claim: {r['claim']}")
    out.append(f"  Loop order / theory: {r['loop_order']} / {r['theory']}")
    out.append(f"  Scheme: {r['scheme']}")
    out.append(f"  Symmetry argument: {r['symmetry_argument']}")
    out.append("  Proved:" + _lst(s.get("proved", [])))
    out.append("  Assumed:" + _lst(s.get("assumed", [])))
    out.append("  Checked partially:" + _lst(s.get("checked_partially", [])))
    out.append("  Evidence:" + _lst(r["evidence"]))
    out.append("  Red flags:" + _lst(r["red_flags"]))
    out.append("  Cross-checks:" + _lst(r["cross_checks"]))
    out.append(f"  Access / confidence: {r['source_access']} / {r['confidence']}")
    out.append(f"  Notes: {r['notes']}")
    return "\n".join(out) + "\n"


def write_record(r, outdir, progress):
    missing = [k for k in KEYS if k not in r]
    if missing:
        raise ValueError(f"{r.get('id')}: missing keys {missing}")
    for k in ("proved", "assumed", "checked_partially"):
        r["symmetry_status"].setdefault(k, [])
    r = {k: r[k] for k in KEYS}
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, r["id"] + ".json"), "w") as f:
        json.dump(r, f, ensure_ascii=False, indent=1)
    with open(os.path.join(outdir, r["id"] + ".txt"), "w") as f:
        f.write(rec_txt(r))
    with open(progress, "a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M} done {r['id']} ({r['arxiv']}) access={r['source_access']}\n")
