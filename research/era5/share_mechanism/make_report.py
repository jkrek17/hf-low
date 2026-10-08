"""Readable tables from results/part1_tests.csv and part2_tests.csv -> results/summary.txt. usage: make_report.py RESULTS_DIR"""
import sys, pandas as pd
R = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 500)
t1 = pd.read_csv(f"{R}/part1_tests.csv"); t2 = pd.read_csv(f"{R}/part2_tests.csv")
o = []
o.append("PART 1 (RA-17), pipeline A, ERA5 proxy. rr per SD of the pattern index; f = share of the log share effect\n")
c = ["variant", "basin", "test", "rr", "rr_lo", "rr_hi", "p", "q", "mde_rr", "f", "f_lo", "f_hi", "A4_rr", "A4_lo", "A4_hi"]
o.append(t1[c].round(3).to_string(index=False))
o.append("\nVerdicts:\n" + t1.groupby(["variant", "basin"]).verdict.first().to_string())
o.append("\nPART 2 (RA-18/RA-6), all variants\n")
o.append(t2.round(3).to_string(index=False))
open(f"{R}/summary.txt", "w").write("\n".join(o))
