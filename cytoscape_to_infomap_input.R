library(dplyr)
library(tidyr)

df1 <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/DEG_All_SLOWP_network04_subset_node.csv", header=TRUE)
nodes_df <- df1 %>% select(name, everything())

df2 <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/DEG_All_SLOWP_network04_subset_edge.csv", header=TRUE)
df2$name <- gsub(" \\(interacts with\\) ", ".", df2$name)
edges_df_prov <- df2 %>% separate(name, c('Source', 'Target'))
edges_df <- edges_df_prov %>% select(Source, Target, everything())

df_All_SLOWP <- edges_df[ , c('Source','Target','irp_score')]

#Input the number of nodes in the network for both IDs here:
ID_Source <- seq(dim(df1)[1])
ID_Target <- seq(dim(df1)[1])
Source <- nodes_df[, "name"]
Source <- sort(Source)
Target <- nodes_df[, "name"]
Target <- sort(Target)
df_Source_ID <- data.frame(Source, ID_Source)
df_Target_ID <- data.frame(Target, ID_Target)

df_Source_ID <- merge(x=df_All_SLOWP,y=df_Source_ID,by="Source")
df_Source_Target_ID <- merge(x=df_Source_ID,y=df_Target_ID,by="Target")
df_IDs_score <- df_Source_Target_ID[, c("ID_Source", "ID_Target", "irp_score")]
#Stuff to get a dataframe with the link LOC-ID called df_LOC_ID

write.table(df_IDs_score, "C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/df_IDs_score.txt", sep=" ", quote=FALSE, row.names=FALSE)
#Do the InfoMap run with the df_IDs_score.txt and dont forget to comment the first line.

write.table(df_Source_Target_ID, "C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/df_Source_Target_ID.txt", sep=" ", quote=FALSE, row.names=FALSE)

#Go to Cytoscape here and add the ID_Source, then the ID_Target (both from the df_Source_Target_ID), the key is the LOC.
#Then, export the node table of the network as (DEG_Testing.csv) and resume the script.
df3 <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/DEG_testing.csv", header=TRUE)
df_LOC_ID <- df3[ , c('name', 'ID_Source', 'ID_Target')]
ID_Source <- df_LOC_ID$ID_Source
ID_Target <- df_LOC_ID$ID_Target
LOC_ID <- coalesce(ID_Source, ID_Target)
df_LOC_ID <- cbind(df_LOC_ID, LOC_ID)
df_LOC_ID <- subset(df_LOC_ID, select = -c(ID_Source, ID_Target))
write.table(df_LOC_ID, "C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/df_LOC_ID.txt", sep=" ", quote=FALSE, row.names=FALSE)

#Treating the InfoMap output for Cytoscape
df_Comunities <- read.csv(file="C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/df_IDs_score.clu", header=FALSE, sep=" ", comment.char = "#")
colnames(df_Comunities) <- c('LOC_ID', 'Module', 'Flow')
df_ID_Comunity <- merge(x=df_LOC_ID,y=df_Comunities,by="LOC_ID")
df_LOC_Comunity <- subset(df_ID_Comunity, select = -c(LOC_ID))
write.table(df_LOC_Comunity, "C:/Users/Hugo Rodrigues/Documents/Hugo_Networks/df_LOC_Comunity.txt", sep=" ", quote=FALSE, row.names=FALSE)
