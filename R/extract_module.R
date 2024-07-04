library(AnnotationDbi)
library(org.Mm.eg.db)

extract_module <- function(module_df, module_number) {
  genes.Module <- module_df %>%
    filter(module == module_number)
  gene.names <- mapIds(org.Mm.eg.db, genes.Module$id, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  gene.names <- na.omit(gene.names)
  gene.names <- c("SymbolName", gene.names)
  
  write(gene.names, file = paste0("namesModule", module_number, ".csv"))
}
