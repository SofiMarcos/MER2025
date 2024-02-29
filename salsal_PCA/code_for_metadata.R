# Script to add metadata to the PCA table
pop <- read.table("Salmo_salar_localities.csv", sep = ";", header = TRUE)

pop$Sample == rownames(df2)

df2$subpopulation <- pop$Subpopulation
df2$population <- pop$Population

