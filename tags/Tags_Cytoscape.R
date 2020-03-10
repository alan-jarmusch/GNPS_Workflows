# INFO --------------------------------------------------------------------
# Tag Cytoscape - Molecular Networking or Feature based Molecular Networking
# Author: Alan K. Jarmusch
# Date: 2020-03-06
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

# COMMON STEPS ------------------------------------------------------------

# Make a copy and download Master Tag Sheet (.tsv) from Google Sheets (https://docs.google.com/spreadsheets/d/1zgSpcgsSxRIgWdHH9khiEMmt2rvNL15raRFP8CxNnFE/edit?usp=sharing).
# INPUT
df_tag_master <- fread("GNPS Tag Template - MASTER - Master.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Open your molecular network in Cytoscape (using graphML) and export the node table
df <- fread("node_file.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
df$INCHI <- gsub('"','',df$INCHI)

# CYTOSCAPE ---------------------------------------------------------------
# Run the following code
## TAGS
tags <- cbind(df_tag_master[,"Inchi_whole string"],separate(df_tag_master[,"TAGS"], TAGS, paste0("X",1:max(sapply(strsplit(df_tag_master$TAGS,"\\|"),length))), sep="\\|"))
tags <- subset(tags, tags$`Inchi_whole string` != "")
tags_reformat <- gather(tags, key="tags", value="terms", 2:length(tags))
tags_reformat <- na.omit(tags_reformat)
tags <- tags_reformat[,-c(2)] %>% group_by(`Inchi_whole string`) %>% distinct()
tags$number <- 1
tags <- spread(tags, key=terms, value=number, fill = FALSE)
tags <- subset(tags, select = -c(V1))
colnames(tags)[-1] <- paste0("TAG_", colnames(tags[,-1]))
colnames(tags)[1] <- "INCHI"

## UBERON
tags_UBERON <- cbind(df_tag_master[,"Inchi_whole string"],separate(df_tag_master[,"UBERONBodyPartName"], UBERONBodyPartName, paste0("X",1:max(sapply(strsplit(df_tag_master$UBERONBodyPartName,"\\|"),length))), sep="\\|"))
tags_UBERON <- subset(tags_UBERON, tags_UBERON$`Inchi_whole string` != "")
tags_reformat <- gather(tags_UBERON, key="tags", value="terms", 2:length(tags_UBERON))
tags_reformat <- na.omit(tags_reformat)
tags_UBERON <- tags_reformat[,-c(2)] %>% group_by(`Inchi_whole string`) %>% distinct()
tags_UBERON$number <- 1
tags_UBERON <- spread(tags_UBERON, key=terms, value=number, fill = FALSE)
tags_UBERON <- subset(tags_UBERON, select = -c(V1))
colnames(tags_UBERON)[-1] <- paste0("TAG_", colnames(tags_UBERON[,-1]))
colnames(tags_UBERON)[1] <- "INCHI"

## LIFESTYLES
tags_Lifestyle <- cbind(df_tag_master[,"Inchi_whole string"],separate(df_tag_master[,"Lifestyle_Tag"], Lifestyle_Tag, paste0("X",1:max(sapply(strsplit(df_tag_master$Lifestyle_Tag,"\\|"),length))), sep="\\|"))
tags_Lifestyle <- subset(tags_Lifestyle, tags_Lifestyle$`Inchi_whole string` != "")
tags_reformat <- gather(tags_Lifestyle, key="tags", value="terms", 2:length(tags_Lifestyle))
tags_reformat <- na.omit(tags_reformat)
tags_Lifestyle <- tags_reformat[,-c(2)] %>% group_by(`Inchi_whole string`) %>% distinct()
tags_Lifestyle$number <- 1
tags_Lifestyle <- spread(tags_Lifestyle, key=terms, value=number, fill = FALSE)
tags_Lifestyle <- subset(tags_Lifestyle, select = -c(V1))
colnames(tags_Lifestyle)[-1] <- paste0("TAG_", colnames(tags_Lifestyle[,-1]))
colnames(tags_Lifestyle)[1] <- "INCHI"

## NCBI Taxonomy
NCBI_tax <- cbind(df_tag_master[,"Inchi_whole string"],separate(df_tag_master[,"NCBITaxonomy"], NCBITaxonomy, paste0("X",1:max(sapply(strsplit(df_tag_master$NCBITaxonomy,":"),length))), sep=":"))
NCBI_tax <- subset(NCBI_tax, NCBI_tax$`Inchi_whole string` != "")
tags_reformat <- gather(NCBI_tax, key="NCBITaxonomy", value="terms", 2:length(NCBI_tax))
tags_reformat <- na.omit(tags_reformat)
NCBI_tax <- tags_reformat[,-c(2)] %>% group_by(`Inchi_whole string`) %>% distinct()
NCBI_tax$number <- 1
NCBI_tax <- spread(NCBI_tax, key=terms, value=number, fill = FALSE)
NCBI_tax <- subset(NCBI_tax, select = -c(V1))
colnames(NCBI_tax)[-1] <- paste0("NCBITaxonomy_", colnames(NCBI_tax[,-1]))
colnames(NCBI_tax)[1] <- "INCHI"

# COMBINE and OUTPUT CYTOSCAPE FILE ---------------------------------------
# Output is an update node table that can be imported into Cytoscape
output <- df %>% left_join(tags, by="INCHI") %>% left_join(NCBI_tax, by="INCHI") %>% left_join(tags_UBERON, by="INCHI") %>% left_join(tags_Lifestyle, by="INCHI")
output[,c((length(df)+1):length(output))] <- replace(output[,c((length(df)+1):length(output))], is.na(output[,c((length(df)+1):length(output))]), 0.1)
output[is.na(output)] <- "N/A"
write.csv(output, "node_file_updated.csv", row.names=FALSE)