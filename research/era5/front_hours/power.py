"""Planted-effect power for P1 (added after the outcomes were joined; uses only the marginal distribution of H, not its relation to the front).
Latent z ~ N(rho*x + sqrt(1-rho^2)*e), x standard normal (the front), y = the observed H values reordered by rank of z. 3,000 sims per rho,
Spearman test at alpha 0.05 (t approximation, n = 19) and at 0.05/11. usage: python3 -I power.py results/season_table.csv"""
import sys
import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(7)
T = pd.read_csv(sys.argv[1], index_col=0).loc[2007:2025]
H = np.sort(T.hours.values); n = len(H)
def power(rho, alpha, N=3000):
    tc = stats.t.ppf(1 - alpha / 2, n - 2); hit = 0
    for _ in range(N):
        x = rng.standard_normal(n); z = rho * x + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
        y = H[np.argsort(np.argsort(z))]
        r = stats.spearmanr(x, y)[0]
        hit += abs(r) * np.sqrt((n - 2) / (1 - r * r)) > tc
    return hit / N
out = ["rho_latent  power(alpha .05)  power(alpha .05/11)"]
for rho in (0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8):
    out.append(f"{rho:.2f}        {power(rho, .05):.2f}              {power(rho, .05/11):.2f}")
print("\n".join(out))
