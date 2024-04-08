### Overview:
Goal is to:
- house the original modules discovery method as well as the improved version
- allow all methods to work on adata objects as well as on the data matrices directly


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
