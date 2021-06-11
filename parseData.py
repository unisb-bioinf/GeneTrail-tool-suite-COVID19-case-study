# This script is used to parse the matrix and additional
# metadata from the .h5ad file provided by Wilk et al.
#

import anndata
import scipy
import scipy.io as scio
import pandas


data = anndata.read_h5ad("Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad", as_sparse=["X"])

with open("matrix.features", "w") as writer:
    writer.write("\n".join(list(data.var.index)))

with open("matrix.barcodes", "w") as writer:
    writer.write("\n".join(list(data.obs.index)))

c = data.obs
pandas.DataFrame(columns=c.columns).to_csv("matrix.metadata", index=False, sep="\t")
c.to_csv("matrix.metadata", header=None, mode='a', sep="\t")

x = scipy.sparse.csr_matrix(data.to_df("matrix").values.T)
scio.mmwrite("matrix.mtx", x)

