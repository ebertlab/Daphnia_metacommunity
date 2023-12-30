dm = readRDS("./distances_matrix.RDS")
mdm = dm
groups = readRDS("./islandsGroups.RDS")

for(i in 1:4){
  for(j in (i+1):5){
    mdm[groups[[i]], groups[[j]]] = 10e15
    mdm[groups[[j]], groups[[i]]] = 10e15
  }
}

write.table(mdm, file = "./distance_matrix_groups", row.names = F, col.names = F)
