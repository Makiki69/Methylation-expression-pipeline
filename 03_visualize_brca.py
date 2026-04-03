import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

meth = pd.read_csv("processed/methylation_gene_level.csv", index_col=0)
expr = pd.read_csv("processed/expression_gene_level.csv", index_col=0)

# Colored scatterplot
merged = pd.DataFrame({
    "methylation": meth.mean(axis=1),
    "expression": expr.mean(axis=1)
})
#KDE contour
sns.kdeplot(
    x=merged["methylation"],
    y=merged["expression"],
    fill=True,
    cmap="viridis"
)
plt.title("Gene-Level Methylation vs Expression (BRCA)")
plt.savefig("processedmeth_vs_expr_density.png", dpi=300)
plt.close()
