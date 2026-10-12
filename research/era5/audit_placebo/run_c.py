"""Part C: machine-readable test register. usage: python3 -I run_c.py OUTDIR. Plan: PREREGISTRATION.md (23f1b5d)."""
import sys, os, re, glob
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
OUT = sys.argv[1]
PC = re.compile(r"^(p|pval|p_value|p_perm|p_boot|p_two|perm_p|pperm|p_season|p_clu\w*|p_raw)$", re.I)
QC = re.compile(r"^(q|q_bh|bh_q|qval|q_fam\w*|q_placebo)$", re.I)
rows = []
for f in sorted(glob.glob(os.path.join(REPO, "research/era5/**/*.csv"), recursive=True)):
    rel = os.path.relpath(f, REPO)
    if "/audit_placebo/" in rel:
        continue
    try:
        d = pd.read_csv(f, nrows=200000)
    except Exception:
        continue
    pc = [c for c in d.columns if PC.match(str(c))]
    qc = [c for c in d.columns if QC.match(str(c))]
    if not pc and not qc:
        continue
    n = int(d[pc[0]].notna().sum()) if pc else int(d[qc[0]].notna().sum())
    npass = int((d[qc[0]] < 0.05).sum()) if qc else None
    rows.append(dict(file=rel, thread=rel.split("/")[2], n_tests=n, has_q=bool(qc), n_q_pass=npass, n_p_lt05=int((d[pc[0]] < 0.05).sum()) if pc else None))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "partC_register.csv"), index=False)
g = R.groupby("thread").agg(files=("file", "count"), tests=("n_tests", "sum"), q_pass=("n_q_pass", "sum")).reset_index()
g.to_csv(os.path.join(OUT, "partC_by_thread.csv"), index=False)
tot = int(R.n_tests.sum())
open(os.path.join(OUT, "partC_summary.txt"), "w").write(
    f"result files with a p or q column: {len(R)}\nthreads (result directories): {R.thread.nunique()}\n"
    f"rows with a p or q value (tests in machine-readable form, lower bound; includes sensitivity variants and sub-family rows): {tot}\n"
    f"rows marked q < 0.05 by their own file: {int(R.n_q_pass.fillna(0).sum())}\n"
    f"rows in files with a q column: {int(R[R.has_q].n_tests.sum())}; rows in files with p only: {int(R[~R.has_q].n_tests.sum())}\n"
    f"expected nominal passes at 0.05 if every row were a global null (crude): {0.05*tot:.0f}\n")
print(open(os.path.join(OUT, "partC_summary.txt")).read()); print(g.sort_values("tests", ascending=False).head(15))
