assign_group <- function(id) {
  output = character()
  for (i in id) {
    if (i %in% cells_id_sub_1) {
      output[i] = "cluster_2.1"
    } else if (i %in% cells_id_sub_2) {
      output[i] = "cluster_2.2"
    } else {
      output[i] = "cluster_1"
    }
  }
  as.factor(output)
}
