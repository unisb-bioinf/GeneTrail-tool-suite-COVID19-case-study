# This script creates pseudobulks for CD14+ monocytes. To this end,
# it operates on the raw count matrix and metadata provided by
# Wilk et al.
#

require(Matrix)
require(muscat)
require(SingleCellExperiment)

data = readMM("matrix.mtx")
metadata = read.table("matrix.metadata", sep="\t", quote="\"")

data = data[,metadata$cell_type_fine == "CD14 Monocyte"]
colnames(data) = rownames(metadata)[metadata$cell_type_fine == "CD14 Monocyte"]
rownames(data) = r[,1]

# This file is used for the single-cell analysis
write.table(as.data.frame(as.matrix(data)), "matrix_count_CD14_Monocytes_all.tsv", sep="\t", quote=F)

sce = SingleCellExperiment(list(abc=data))
colData(sce)$sample_id = metadata$Donor_full[metadata$cell_type_fine == "CD14 Monocyte"]
colData(sce)$cluster_id = metadata$cell_type_fine[metadata$cell_type_fine == "CD14 Monocyte"]
colData(sce)$group_id = metadata$Ventilated[metadata$cell_type_fine == "CD14 Monocyte"]

r = read.table("matrix.features", stringsAsFactors=F)

pb = aggregateData(sce)

rownames(pb) = r[,1]
colnames(pb) = c("C1_A", "C1_B", "C2", "C3", "C4", "C5", "C6", "C7", "H1", "H2", "H3", "H4", "H5", "H6")

x = assay(pb, 1)
c = colSums(x)
m = median(c)
for(i in 1:ncol(x)){
    x[,i] = x[,i]/c[i]
}
x = x * m

# This file is used for the REGGAE analysis
write.table(x, file="pseudobulk_CD14_Monocyte_All.normalized.tsv", sep="\t", quote=F)

