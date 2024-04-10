### Overview:
Goal is to:
- house the original modules discovery method as well as the improved version
- allow all methods to work on adata objects as well as on the data matrices directly


### Breakdown:
<!-- TODO: Add the corresponding code statement or options for each -->
1. Normalize (outside of modules)
2. Perform SVD (or get pre-computed SVD from adata object)
    - Optional filter components with strong batch effects
    - Optional RUV/SVA on the PC space
3. PC adjustment
4. Estimating the correlation
5. Build adjacency matrix
6. Filter k-NN (thresholding)
7. Perform community detection (NMF, BigClam, Leiden)
8. Downstream tasks
9. Benchmarking

### Usage
```python
import scdemon as sm
from scdemon import plotting as pl

mo = sm.run_modules(adata)
pl.plot_genes(mo, 'base')

```


### Old todolist / desirables:
Done
- Fixed graph layout (nogrid works)
- Enrichments on graph / other attributes on graph
- P-values on corr for cutoff
- Quantiles for corr. cutoff
- Expr percent dpdt. cutoffs
- Writing module list

Todo:
- Saving enrichments so we don't recompute them
- Heatmap for the enrichments (make profiles, pass to scanpy dotplot)
- Compute gene properties
- Plot umap with any coloring (gene properties)
- Handle both X and adata
- Plot multiple graphs from the same adata
- Better docstrings
- Plot multiple graphs on the same layout (for subsetting!)
- Cell state discovery + annotation from gene expression modules.
