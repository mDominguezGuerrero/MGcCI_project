filter_cells <- function(Seurat) {
  ## Get filtering parameters
  count.max <- round(mean(Seurat$nCount_RNA) + 2 * sd(Seurat$nCount_RNA), digits = -2)
  count.min <- round(mean(Seurat$nCount_RNA) - 2 * sd(Seurat$nCount_RNA), digits = -2)
  feat.max <- round(mean(Seurat$nFeature_RNA) + 2 * sd(Seurat$nFeature_RNA), digits = -2)
  feat.min <- round(mean(Seurat$nFeature_RNA) - 2 * sd(Seurat$nFeature_RNA), digits = -2)
  ## Set minimum parameters to 0 if negative value
  if (count.min < 0){
    count.min <- 0
  } else {
    count.min <- count.min
  }
  
  if (feat.min < 0){
    feat.min <- 0
  } else {
    feat.min <- feat.min
  }
  ## Filter cells
  Seurat <- subset(Seurat, 
                   subset = nFeature_RNA > feat.min & nFeature_RNA < feat.max & nCount_RNA < count.max & nCount_RNA > count.min)
  
}
