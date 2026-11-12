"""
CYP2C19 / clopidogrel analysis
Input : cyp2c19_final.tsv  (made in the previous steps)
Output: pop_table.tsv, tests.tsv, figure1.png
"""
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # lets matplotlib save images without a screen
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# ---------- 1. Load data and keep only classifiable people ----------
df = pd.read_csv("cyp2c19_final.tsv", sep="\t")
print("People in file:", len(df))

n_unclass = (df["group"] == "Not classifiable").sum()
df = df[df["group"] != "Not classifiable"].copy()
df["alt"] = (df["group"] == "Alternative advised").astype(int)
print("Not classifiable (excluded):", n_unclass)
print("Classifiable people used:", len(df))


# ---------- 2. Helper functions ----------
def wilson_ci(k, n, z=1.96):
    """95% confidence interval for a proportion k/n (Wilson method)."""
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return centre - half, centre + half


def odds_ratio_ci(a, b, c, d, z=1.96):
    """Odds ratio and 95% CI for the table [[a, b], [c, d]] (Woolf method).
    Adds 0.5 to every cell if any cell is zero."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return or_, math.exp(math.log(or_) - z * se), math.exp(math.log(or_) + z * se)


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values (multiple-testing correction)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        adj[i] = prev
    return adj


# ---------- 3. Per-population table ----------
rows = []
for pop, g in df.groupby("pop"):
    n = len(g)
    k = int(g["alt"].sum())
    lo, hi = wilson_ci(k, n)
    rows.append(
        {
            "pop": pop,
            "super_pop": g["super_pop"].iloc[0],
            "n_classifiable": n,
            "n_alt": k,
            "pct_alt": round(100 * k / n, 1),
            "ci_low": round(100 * lo, 1),
            "ci_high": round(100 * hi, 1),
        }
    )
pop_table = pd.DataFrame(rows)
order = {"SAS": 0, "EUR": 1, "EAS": 2, "AFR": 3}
pop_table["o"] = pop_table["super_pop"].map(order)
pop_table = pop_table.sort_values(["o", "pct_alt"], ascending=[True, False]).drop(columns="o")
pop_table.to_csv("pop_table.tsv", sep="\t", index=False)
print("\n=== % advised an alternative to clopidogrel, by population ===")
print(pop_table.to_string(index=False))


# ---------- 4. Statistical tests (Fisher's exact test) ----------
def counts(mask):
    sub = df[mask]
    return int(sub["alt"].sum()), int(len(sub) - sub["alt"].sum())


tests = []
sas = df["super_pop"] == "SAS"

# (a) all South Asians pooled vs each other super-population
for other in ["EUR", "EAS", "AFR"]:
    a, b = counts(sas)
    c, d = counts(df["super_pop"] == other)
    odds, lo, hi = odds_ratio_ci(a, b, c, d)
    _, p = fisher_exact([[a, b], [c, d]])
    tests.append(("SAS (pooled)", other + " (pooled)", odds, lo, hi, p))

# (b) each South Asian population vs all Europeans pooled
eur = df["super_pop"] == "EUR"
for pop in sorted(df.loc[sas, "pop"].unique()):
    a, b = counts(df["pop"] == pop)
    c, d = counts(eur)
    odds, lo, hi = odds_ratio_ci(a, b, c, d)
    _, p = fisher_exact([[a, b], [c, d]])
    tests.append((pop, "EUR (pooled)", odds, lo, hi, p))

adj = bh_adjust([t[5] for t in tests])
tests_df = pd.DataFrame(
    [
        {
            "group_1": t[0],
            "group_2": t[1],
            "odds_ratio": round(t[2], 2),
            "ci_low": round(t[3], 2),
            "ci_high": round(t[4], 2),
            "p_value": f"{t[5]:.2e}",
            "p_adj_BH": f"{q:.2e}",
        }
        for t, q in zip(tests, adj)
    ]
)
tests_df.to_csv("tests.tsv", sep="\t", index=False)
print("\n=== Fisher's exact tests (odds of needing an alternative) ===")
print(tests_df.to_string(index=False))


# ---------- 5. Figure 1 ----------
colours = {"SAS": "#c0392b", "EUR": "#2980b9", "EAS": "#27ae60", "AFR": "#8e44ad"}
plot_df = pop_table.iloc[::-1].reset_index(drop=True)  # reverse so SAS is on top

fig, ax = plt.subplots(figsize=(8, 9))
for i, r in plot_df.iterrows():
    ax.barh(i, r["pct_alt"], color=colours[r["super_pop"]], alpha=0.9)
    ax.errorbar(
        r["pct_alt"], i,
        xerr=[[r["pct_alt"] - r["ci_low"]], [r["ci_high"] - r["pct_alt"]]],
        color="black", capsize=2, linewidth=1,
    )
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df["pop"])
ax.set_xlabel("% of people for whom CPIC advises an alternative to clopidogrel")
ax.set_title("CYP2C19 clopidogrel risk by population\n(1000 Genomes high-coverage, 95% CI)")
ax.set_xlim(0, 100)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colours.values()]
ax.legend(handles, ["South Asian", "European", "East Asian", "African"], loc="lower right")
plt.tight_layout()
plt.savefig("figure1.png", dpi=200)
print("\nSaved: pop_table.tsv, tests.tsv, figure1.png")
