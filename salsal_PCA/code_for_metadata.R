# Script to add metadata to the PCA table
# Load metadata
pop <- read.table("Salmo_salar_localities.csv", sep = ";", header = TRUE)

# Check if sample order in the metadata file is the same as in the pca dataframe
pop$Sample == rownames(df2)

# Add two columns to the pca dataframe
df2$subpopulation <- pop$Subpopulation
df2$population <- pop$Population

