# plaque_atlas_signature_matrix.py

import scanpy as sc
import pandas as pd
import numpy as np
import os

# =====================================================
# 1️⃣ Load AnnData
# =====================================================
adata = sc.read_h5ad("3bab1c0b-d3e3-4a01-840f-d49a8284d989.h5ad")

# Inspect structure
print(adata)

# =====================================================
# 2️⃣ Choose cell-type annotation
# =====================================================
# You can use 'cell_type_level1' (broad) or 'cell_type_level2' (fine-grained)
adata.obs['ct'] = adata.obs['cell_type_level1']
cell_types = sorted(adata.obs['ct'].unique().tolist())
print(f"Using {len(cell_types)} cell types at Level 1 resolution")

# =====================================================
# 3️⃣ Use pre-normalized expression
# =====================================================
# Plaque Atlas already includes log-normalized data
adata.X = adata.layers["log_normalized_individually"]

# =====================================================
# 4️⃣ Compute pseudo-bulk mean expression per cell type
# =====================================================
genes = adata.var_names
pseudo_bulk = pd.DataFrame(index=genes)

for ct in cell_types:
    cells = adata.obs_names[adata.obs['ct'] == ct]
    subset = adata[cells, :]
    pseudo_bulk[ct] = np.asarray(subset.X.mean(axis=0)).ravel()

# Convert from log1p to linear scale
pseudo_bulk_linear = np.expm1(pseudo_bulk)

# =====================================================
# 5️⃣ Filter for cell-type–specific genes
# =====================================================
# Compute fold-change between top and second-highest expression
max1 = pseudo_bulk_linear.max(axis=1)
max2 = pseudo_bulk_linear.apply(lambda x: np.partition(x.values, -2)[-2], axis=1)
fc = (max1 + 1e-9) / (max2 + 1e-9)

# Keep genes with FC ≥ 2 (specific markers)
marker_genes = fc[fc >= 2].index
signature = pseudo_bulk_linear.loc[marker_genes]

print(f"Selected {len(marker_genes)} marker genes")

# =====================================================
# 6️⃣ Normalize each column (cell-type profile sums to 1)
# =====================================================
signature_norm = signature.div(signature.sum(axis=0), axis=1)

# =====================================================
# 7️⃣ (Optional) Convert Ensembl IDs → Gene Symbols
# =====================================================
if 'feature_name' in adata.var.columns:
    mapping = adata.var[['feature_name']].copy()
    mapping.index = adata.var_names
    mapping = mapping.dropna()
    signature_norm = signature_norm.rename(index=mapping['feature_name'].to_dict())
    signature_norm = signature_norm[~signature_norm.index.duplicated(keep='first')]
    print("Mapped Ensembl IDs to gene symbols using 'feature_name' column.")

# =====================================================
# 8️⃣ Save signature matrix
# =====================================================
output_file = "plaque_atlas_matrix.csv"
signature_norm.to_csv(output_file)
print(f"✅ Signature matrix saved to: {output_file}")

