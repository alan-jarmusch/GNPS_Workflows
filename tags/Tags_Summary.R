# INFO --------------------------------------------------------------------
# Tag Summary - Molecular Networking or Feature based Molecular Networking
# Author: Alan K. Jarmusch
# Date: 2020-03-06
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# COMMON STEPS ----------------------------------------------------------

# Make a copy and download Master Tag Sheet (.tsv) from Google Sheets (https://docs.google.com/spreadsheets/d/1zgSpcgsSxRIgWdHH9khiEMmt2rvNL15raRFP8CxNnFE/edit?usp=sharing).
# INPUT
df_tag_master <- fread("GNPS Tag Template - MASTER - Master.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# INPUT:"View all library hits" table (downloaded from GNPS) 
df_hits <- fread("METABOLOMICS-SNETS-V2-e847e3b6-view_all_annotations_DB-main.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
df_hits$INCHI <- gsub('"','',df_hits$INCHI)

# Distinct INCHI
df_hits <- df_hits %>% group_by(INCHI) %>% distinct()
df_hits <- subset(df_hits, df_hits$INCHI != "" & df_hits$INCHI != "n/a" & df_hits$INCHI != "N/A")


# MOLECULAR NETWORKING - TAGS ---------------------------------------------
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
df_hit_tags <- merge(df_hits, tags, by="INCHI")

if (length(df_hit_tags) > length(df_hits)) {
  df_hit_tags_summary <- df_hit_tags[, c((length(df_hits)+1):length(df_hit_tags))] %>% summarise_all(funs(sum))
  df_hit_tags_table <- gather(df_hit_tags_summary, tag, number, 1:length(df_hit_tags_summary))
  df_hit_tags_table <- subset(df_hit_tags_table, df_hit_tags_table$number != "0")
  df_hit_tags_table$tag <- gsub("TAG_", "", df_hit_tags_table$tag)
  
  plot_num_annotations <- ggplot(data=df_hit_tags_table, aes(x= reorder(as.factor(tag),number), y=as.numeric(number)))+
    stat_summary_bin(width= 0.75, fill="grey50", fun.y = "sum", geom = "bar")+
    geom_text(aes(label=as.numeric(number)), hjust=-0.5, size=1.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(df_hit_tags_table$number))+max(as.numeric(df_hit_tags_table$number))*0.05))+
    coord_flip()+
    theme_minimal()+
    theme(panel.grid.major=element_line(colour ="grey75",size=0.25, linetype="dashed"),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="none",
          legend.title = element_text(colour="black", size=6),
          legend.text = element_text(colour="black", size=6)) +
    labs(x="Tags", y="Number")
  print(plot_num_annotations)
  ggsave(file="GNPS_Summary_TAG.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
}  else {
  print("no tags")
}


# MOLECULAR NETWORKING - UBERON -------------------------------------------
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
df_hit_tags_UBERON <- merge(df_hits, tags_UBERON, by="INCHI")
if (length(df_hit_tags_UBERON) > length(df_hits)) {
  df_hit_tags_UBERON_summary <- df_hit_tags_UBERON[, c((length(df_hits)+1):length(df_hit_tags_UBERON))] %>% summarise_all(funs(sum))
  df_hit_tags_UBERON_table <- gather(df_hit_tags_UBERON_summary, tags_UBERON, number, 1:length(df_hit_tags_UBERON_summary))
  df_hit_tags_UBERON_table <- subset(df_hit_tags_UBERON_table, df_hit_tags_UBERON_table$number != "0")
  df_hit_tags_UBERON_table$tags_UBERON <- gsub("TAG_", "", df_hit_tags_UBERON_table$tags_UBERON)
  
  plot_num_annotations <- ggplot(data=df_hit_tags_UBERON_table, aes(x= reorder(as.factor(tags_UBERON),number), y=as.numeric(number)))+
    stat_summary_bin(width= 0.75, fill="grey50", fun.y = "sum", geom = "bar")+
    geom_text(aes(label=as.numeric(number)), hjust=-0.5, size=1.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(df_hit_tags_UBERON_table$number))+max(as.numeric(df_hit_tags_UBERON_table$number))*0.05))+
    coord_flip()+
    theme_minimal()+
    theme(panel.grid.major=element_line(colour ="grey75",size=0.25, linetype="dashed"),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="none",
          legend.title = element_text(colour="black", size=6),
          legend.text = element_text(colour="black", size=6)) +
    labs(x="tags_UBERON", y="Number")
  print(plot_num_annotations)
  ggsave(file="GNPS_Summary_UBERON.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
}  else {
  print("no tags")
}

# MOLECULAR NETWORKING - LIFESTYLE ----------------------------------------
## LIFESTYLE
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
df_hit_tags_Lifestyle <- merge(df_hits, tags_Lifestyle, by="INCHI")
if (length(df_hit_tags_Lifestyle) > 40) {
  df_hit_tags_Lifestyle_summary <- df_hit_tags_Lifestyle[, c((length(df_hits)+1):length(df_hit_tags_Lifestyle))] %>% summarise_all(funs(sum))
  df_hit_tags_Lifestyle_table <- gather(df_hit_tags_Lifestyle_summary, Lifestyle_Tag, number, 1:length(df_hit_tags_Lifestyle_summary))
  df_hit_tags_Lifestyle_table <- subset(df_hit_tags_Lifestyle_table, df_hit_tags_Lifestyle_table$number != "0")
  df_hit_tags_Lifestyle_table$Lifestyle_Tag <- gsub("TAG_", "", df_hit_tags_Lifestyle_table$Lifestyle_Tag)
  
  plot_num_annotations <- ggplot(data=df_hit_tags_Lifestyle_table, aes(x= reorder(as.factor(Lifestyle_Tag),number), y=as.numeric(number)))+
    stat_summary_bin(width= 0.75, fill="grey50", fun.y = "sum", geom = "bar")+
    geom_text(aes(label=as.numeric(number)), hjust=-0.5, size=1.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(df_hit_tags_Lifestyle_table$number))+max(as.numeric(df_hit_tags_Lifestyle_table$number))*0.05))+
    coord_flip()+
    theme_minimal()+
    theme(panel.grid.major=element_line(colour ="grey75",size=0.25, linetype="dashed"),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="none",
          legend.title = element_text(colour="black", size=6),
          legend.text = element_text(colour="black", size=6)) +
    labs(x="Lifestyle_Tag", y="Number")
  print(plot_num_annotations)
  ggsave(file="GNPS_Summary_Lifestyle.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
}  else {
  print("no tags")
}

# MOLECULAR NETWORKING - NCBI Taxonomy ------------------------------------
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
df_hit_NCBI_tax <- merge(df_hits, NCBI_tax, by="INCHI")
if (length(df_hit_NCBI_tax) > length(df_hits)) {
  df_hit_NCBI_tax_summary <- as.numeric(df_hit_NCBI_tax[, c((length(df_hits)+1):length(df_hit_NCBI_tax))]) %>% summarise_all(funs(sum))
  df_hit_NCBI_tax_table <- gather(df_hit_NCBI_tax_summary, NCBI_tax, number, 1:length(df_hit_NCBI_tax_summary))
  df_hit_NCBI_tax_table <- subset(df_hit_NCBI_tax_table, df_hit_NCBI_tax_table$number != "0")
  
  plot_num_annotations <- ggplot(data=df_hit_NCBI_tax_table, aes(x= reorder(as.factor(NCBI_tax),number), y=as.numeric(number)))+
    stat_summary_bin(width= 0.75, fill="grey50", fun.y = "sum", geom = "bar")+
    geom_text(aes(label=as.numeric(number)), hjust=-0.5, size=1.5)+
    scale_y_continuous(limits=c(0,max(as.numeric(df_hit_NCBI_tax_table$number))+max(as.numeric(df_hit_NCBI_tax_table$number))*0.05))+
    coord_flip()+
    theme_minimal()+
    theme(panel.grid.major=element_line(colour ="grey75",size=0.25, linetype="dashed"),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="none",
          legend.title = element_text(colour="black", size=6),
          legend.text = element_text(colour="black", size=6)) +
    labs(x="NCBI_tax", y="Number")
  print(plot_num_annotations)
  ggsave(file="GNPS_Summary_NCBI.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)
}  else {
  print("no tags")
}
