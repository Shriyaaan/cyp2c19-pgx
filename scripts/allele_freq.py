"""
CYP2C19 allele frequencies per population
Input : cyp2c19_final.tsv
Output: allele_freq.tsv, figure2.png
"""
import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("cyp2c19_final.tsv", sep="\t")

# keep only rows with a real diplotype like *1/*2 (drops "no call")
ok = df["diplotype"].astype(str).str.match(r"^\*\d+/\*\d+$")
print("People with a called diplotype:", int(ok.sum()), "of", len(df))
df = df[ok].copy()

ALLELES = ["*2", "*3", "*17"]
rows = []
for pop, g in df.groupby("pop"):
    alleles = g["diplotype"].str.split("/").explode()   # 2 alleles per person
    total = len(alleles)                                  # = 2 x people
    row = {"pop": pop, "super_pop": g["super_pop"].iloc[0], "n_people": len(g)}
    for a in ALLELES:
        row[a] = round(100 * (alleles == a).sum() / total, 1)
    row["other"] = round(100 * (~alleles.isin(["*1"] + ALLELES)).sum() / total, 1)
    rows.append(row)

out = pd.DataFrame(rows)
order = {"SAS": 0, "EUR": 1, "EAS": 2, "AFR": 3}
out["o"] = out["super_pop"].map(order)
out = out.sort_values(["o", "*2"], ascending=[True, False]).drop(columns="o")
out.to_csv("allele_freq.tsv", sep="\t", index=False)
print("\n=== Allele frequency (% of alleles) ===")
print(out.to_string(index=False))

# ---------- Figure 2: heatmap ----------
fig, ax = plt.subplots(figsize=(6, 9))
data = out[ALLELES].values
im = ax.imshow(data, cmap="Reds", aspect="auto")
ax.set_xticks(range(len(ALLELES)))
ax.set_xticklabels(["*2\n(no function)", "*3\n(no function)", "*17\n(increased)"])
ax.set_yticks(range(len(out)))
ax.set_yticklabels([f"{p} ({s})" for p, s in zip(out["pop"], out["super_pop"])])
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        v = data[i, j]
        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                color="white" if v > data.max() * 0.6 else "black", fontsize=8)
ax.set_title("CYP2C19 allele frequency (% of alleles)")
plt.colorbar(im, ax=ax, shrink=0.6, label="%")
plt.tight_layout()
plt.savefig("figure2.png", dpi=200)
print("\nSaved: allele_freq.tsv, figure2.png")
