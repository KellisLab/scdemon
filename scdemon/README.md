### Overview:
Goal is to:
- house the original modules discovery method as well as the improved version
- allow all methods to work on adata objects as well as on the data matrices directly


### Breakdown:
**Make object:**
```python
mod = sm.modules_core(adata, suffix='tagline',
        estimate_sd=False, svd_k=max_k,
        filter_expr=0.05, z=4.5, 
        calc_raw=False)
```

**Setup:**
1. Filter dataset by genes (normalize outside)
2. Perform SVD (or get pre-computed SVD from adata object)
3. Flag PCs correlated with specific covariates
4. Adjust PCs using `robust_prepare`
```python
mod.setup()  # Setup the object
```

**Make graph:**
```python
mod.make_graph(graph_id='base',
        multigraph=False,
        filter_covariate=None)
```
In `framework`:
3. Select PCs to use
4. Estimate the correlation

Within the specific `graph` object:
5. Build adjacency matrix
6. Filter k-NN (thresholding)
    - Using bivariate splines
    - Using `robust_SE_t`
7. Perform community detection
    - Leiden
    - Hierarchical clustering
    - NMF / Factor Analysis
    - BigClam
    - Stoch. block model

Outside of module or graph objects:
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
