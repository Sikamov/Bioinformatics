########### library
#library(ggrepel)
library(magrittr)
library(dplyr)   
library(ggplot2)
#library(tximport)
library(tidyr)
library(ggbreak) 
#library(biomaRt)
library(rtracklayer)
library(openxlsx)
library(FactoMineR)
library(tidyverse)
library(org.EcK12.eg.db)
library(BiocManager)
library(AnnotationDbi)
library(KEGGREST)
library(VennDiagram)
library(readxl)
library(biomaRt)
library(httr)
library(tidyverse)
library(Biostrings)
library(UpSetR)
library(ComplexUpset)
library(read.xlsx)
library(missMDA)
#BiocManager::install("missMDA", force = TRUE)





##########directory
dir <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/H34/"
annotation<- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Mh-H34_filt_bold.tsv"
names <- c("33", "43", "11", "1862", "12", "X-37", "7", "45", "40")



setwd(dir)
#png(paste0(dir,"vulcanoplots/vulcanoplots.png"),width=1000,height=1000,res=200)
#par(mfrow = c(3, 3))
wb <- createWorkbook()
i <-0
Protein <- NA
repeats <- NA
ID <- NA
df <- data.frame(ID)
df_PCA <- data.frame(ID)



#############################################################################################################

for (name in names) {
    setwd(dir)
   
    #######data
    data_all <- read.csv(paste0(dir, "mycoplasma_design_Mh_", name, "_openms_design_msstats_in_comparisons.csv")
                         , sep = "\t")
    colnames(data_all)[colnames(data_all) == "Protein"] <- "ID"
    data_all <- data_all %>% 
      mutate(
        Expression = case_when(log2FC >= 1 & adj.pvalue <= 0.05 & log2FC != Inf ~ "Up-regulated",
                               log2FC <= -1 & adj.pvalue <= 0.05 & log2FC != -Inf ~ "Down-regulated",
                               log2FC == -Inf ~ "Unique-down-regulated",
                               log2FC == Inf ~ "Unique-up-regulated",
                               TRUE ~ "Unchanged")
      )
    data <- data_all[!is.infinite(data_all$log2FC), ]
    n <- data %>% count(Expression) 
    diff <- data[data$adj.pvalue<=0.05 & abs(data$log2FC)>=1, ]
    diff <- diff[complete.cases(diff$log2FC), ]
  
    
    ###df
    df1 <- data_all[ ,c("ID","log2FC")]
    df <- merge(df, df1, by = "ID", all.y = TRUE, all.x = TRUE)
    colnames(df)[colnames(df) == "log2FC"] <- paste0("Mh",name,"/H34")
    
    
    ###df_PCA
    data_PCA<- read.csv(paste0(dir,"PCA/", "mycoplasma_design_Mh_", name, "_openms_design_msstats_in.csv"), stringsAsFactors = FALSE, sep = ",", dec = ".")
                              
    data_PCA <- data_PCA[c("ProteinName", "Intensity", "Reference")] #, "PeptideSequence"
    #data_PCA <- mutate(data_PCA, Reference = as.character(Reference), Intensity = as.character(Intensity))
    data_PCA <- data_PCA%>%
      group_by(ProteinName, Reference) %>%
      mutate(Intensity = sum(Intensity)) %>%
      distinct(ProteinName, Reference, .keep_all = TRUE) %>%
      pivot_wider(names_from = Reference, values_from = Intensity) %>%
      mutate(across(starts_with("Intensity"), ~ as.numeric(.))) 
 
    colnames(data_PCA)[colnames(data_PCA) == "ProteinName"] <- "ID"
    colnames(data_PCA) <- gsub("OB0(\\d{3}).mzML", "\\1", colnames(data_PCA))
    if (name != "40") {
      data_PCA <- data_PCA[, !(names(data_PCA) %in% c("707","706","708", "PeptideSequence"))]
    }
    df_PCA <- merge(data_PCA, df_PCA, by = "ID", all.y = TRUE, all.x = TRUE)
    
    
    ######annotations_BAKTA
    db_bakta <- read.table(annotation, header = TRUE, quote="", sep = "\t")
    extract_value <- function(string, prefix) {
      parts <- unlist(strsplit(string, ", "))
      value <- NA
      for (part in parts) {
        if (grepl(prefix, part)) {
          value <- sub(paste0("^", prefix), "", part)
          break
        }
      }
      return(value)
    }
    db_bakta$RefSeq <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "RefSeq:"))
    db_bakta$SO <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "SO:"))
    db_bakta$UniParc <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "UniParc:"))
    db_bakta$UniRef <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "UniRef:"))
    db_bakta$UserProtein <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "UserProtein:"))
    db_bakta$GO <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "GO:"))
    db_bakta$KEGG <- sapply(db_bakta$DbXrefs, function(x) extract_value(x, "KEGG:"))
    
    
    lines <- readLines("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/GCF_000085865.1_ASM8586v1_genomic.faa")
    locus_tag <- vector()
    old_locus_tag <- vector()
    
    for (k in 1:length(lines)) {
      if (grepl("/locus_tag", lines[k])) {
        locus_tag_value <- sub('^.*"([^"]+)".*', '\\1', lines[k])
        locus_tag <- c(locus_tag, locus_tag_value)
        if (grepl("/old_locus_tag", lines[k+1])) {
          old_locus_tag_value <- sub('^.*"([^"]+)".*', '\\1', lines[k+1])
          old_locus_tag <- c(old_locus_tag, old_locus_tag_value)
        }
        else {
          old_locus_tag_value <- NA
          old_locus_tag <- c(old_locus_tag, old_locus_tag_value)
        }
      }
    }
    
    old_new <- data.frame(locus_tag = locus_tag, old_locus_tag = old_locus_tag)
    old_new  <- old_new [!duplicated(old_new ),]
    colnames(old_new)[colnames(old_new) == "locus_tag"] <- "UserProtein"
    
    db_bakta <- merge(db_bakta, old_new, by = "UserProtein", all.x = TRUE)
    db_bakta <- merge(db_bakta, blast, by = "Locus.Tag", all.x = TRUE)
    colnames(db_bakta)[colnames(db_bakta) == "old_locus_tag"] <- "OldUserProtein"
    
    for (j in 1:nrow(db_bakta)) {
      if (is.na(db_bakta$Gene[j])||db_bakta$Gene[j] == "") {
        if (grepl("Hypothetical protein", db_bakta$Product[j])||grepl("hypothetical protein", db_bakta$Product[j])||is.na(db_bakta$Product[j])||db_bakta$Product[j] == "") {
          if (grepl("Hypothetical protein", db_bakta$BLAST[j])||grepl("hypothetical protein", db_bakta$BLAST[j])||is.na(db_bakta$BLAST[j])||db_bakta$BLAST[j] == "") {
            if (is.na(db_bakta$UniParc[j])||db_bakta$UniParc[j] == "") {
              db_bakta$Protein[j] <- db_bakta$Product[j]
            } else { db_bakta$Protein[j] <- db_bakta$UniParc[j]
            }
          } else {
            db_bakta$Protein[j] <- db_bakta$BLAST[j]
          }
        } else {
          db_bakta$Protein[j] <- db_bakta$Product[j]
        }
      } else {
        db_bakta$Protein[j] <- db_bakta$Gene[j]
      }
    }
    
    my_txdf <- db_bakta
    colnames(my_txdf)[colnames(my_txdf) == "Locus.Tag"] <- "ID"
    # my_txdf <- as.data.frame(my_txdf[,c("ID","UserProtein","Gene","Product","RefSeq","GO","KEGG","old_locus_tag")])
    # my_txdf <- my_txdf[!duplicated(my_txdf),]
    # colnames(my_txdf)[colnames(my_txdf) == "old_locus_tag"] <- "OldUserProtein"
    # my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
    # #my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]

    # for (j in 1:nrow(my_txdf)) {
    #   if (is.na(my_txdf$Gene[j])||my_txdf$Gene[j] == "") {
    #     if (grepl("Hypothetical protein", my_txdf$Product[j])||grepl("hypothetical protein", my_txdf$Product[j])||is.na(my_txdf$Product[j])||my_txdf$Product[j] == "") {
    #       my_txdf$Protein[j] <- my_txdf$OldUserProtein[j]
    #     } else {
    #       my_txdf$Protein[j] <- my_txdf$Product[j]
    #     }
    #   } else {
    #     my_txdf$Protein[j] <- my_txdf$Gene[j]
    #   }
    # }
    
    
    ######annotations_uniprot
    # my_txdf <- db_uniprot
    # colnames(my_txdf)[colnames(my_txdf) == "Entry"] <- "ID"
    # my_txdf <- as.data.frame(my_txdf[,c("ID","Protein.names","gene")])
    # my_txdf <- my_txdf[!duplicated(my_txdf),]
    # colnames(my_txdf)[colnames(my_txdf) == "Protein.names"] <- "product"
    # my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
    # #my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]
    # my_txdf$Protein <- my_txdf$gene
    # 
    # 
    # 
    diff <- merge(diff, my_txdf, by = "ID", all.x = TRUE)
    data_all <- merge(data_all, my_txdf, by = "ID", all.x = TRUE)

   
    write.table(as.data.frame(my_txdf), file = "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Annotation.txt", row.names = TRUE, col.names = TRUE)
    
    
    #########xlsx statistics
    addWorksheet(wb, sheetName = paste0("Mh",name,"_H34_diff"))
    i<-i+1
    writeData(wb, sheet = i, x = diff)
    addWorksheet(wb, sheetName = paste0("Mh",name,"_H34_all"))
    i<-i+1
    writeData(wb, sheet = i, x = data_all)
    
    if (i==2*length(names)) {
      df_common<- df[complete.cases(df), ]
      df <- merge(df, my_txdf, by = "ID", all.x = TRUE)
      df_common<- merge(df_common, my_txdf, by = "ID", all.x = TRUE)
      addWorksheet(wb, sheetName = "log2FC_all")
      writeData(wb, sheet = 2*length(names)+1, x = df)
    }
    
    saveWorkbook(wb, file = paste0(dir,"statistics_mycoplasma_proteomes.xlsx"), overwrite = TRUE)
    
    
    
  
    #########vulcanoplots
    dir.create(paste0(dir,"vulcanoplots/"))
    png(paste0(dir,"vulcanoplots/vulcanoplot_Mh",name,"_H34.png"),width=1000,height=1000,res=200)
    vulcanoplot<- ggplot(data, aes(log2FC, -log(adj.pvalue,10))) +
      geom_point(aes(color = Expression), size = 1) +
      xlab(expression("log"[2]*"FC")) + 
      ylab(expression("-log"[10]*"FDR")) +
      scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
      guides(colour = guide_legend(override.aes = list(size=1.5))) + 
      ggtitle(paste0("Mh-",name,"/H-34")) + theme_bw() +
      theme(legend.title = element_blank(),
            legend.text = element_text(size=10),
            legend.position = "top",
            panel.grid.major = element_line(colour=NA),
            panel.grid.minor = element_line(colour=NA))  +
      geom_hline(yintercept=-log(0.05,10), linetype="dashed", color = "gray50")  +
      geom_vline(xintercept=-1, linetype="dashed", color = "dodgerblue3") +
      geom_vline(xintercept=1, linetype="dashed", color = "firebrick3") +
      annotate(geom='label', x =c(-1.5,1.5), 
               y = max(-log(data$adj.pvalue,10),na.rm = TRUE, finite = TRUE), 
               label = c(n$n[1],n$n[3]), col=c("dodgerblue3","firebrick3"), size=3) 
    
    if (name=="") {
      vulcanoplot <- vulcanoplot + scale_x_break(breaks = c(-45, -10)) +
        scale_x_continuous(breaks = seq(-50, 10, by = 5), limits = c(-50, 10)) +
        theme(axis.title.x.top = element_blank(), axis.text.x.top = element_blank(),
              axis.ticks.x.top = element_blank())
    }
      
    print(vulcanoplot)
    dev.off()
    
    
    
    
    
    #########barplots
    dir.create(paste0(dir,"barplots/"))
    png(paste0(dir,"barplots/barplots_Mh",name,"_H34.png"),width=3000,height=2000,res=200) 
    
    barplots<-ggplot(diff, aes(x=ID, y=log2FC, fill = Expression)) + 
      geom_bar(stat = "identity", color="black") + 
      scale_fill_manual(values = c("Up-regulated" = "firebrick3", "Down-regulated" = "dodgerblue3", "Unchanged" = "gray50")) +
      theme_bw() +
      xlab(expression("Proteins")) + 
      ylab(expression("log"[2]*"FC")) + 
      ggtitle(paste0("Mh-",name,"/H-34")) +
      theme(legend.title = element_blank(),
            legend.text = element_text(size=10),
            legend.position = "top",
            panel.grid.major = element_line(colour=NA),
            panel.grid.minor = element_line(colour=NA),
            axis.text.x = element_text(angle=90, size=7, hjust = 1),
            axis.text.y = element_text(angle = 90)) +
      scale_x_discrete(labels = diff$Protein) +
      geom_errorbar(aes(ymin = log2FC - SE, ymax = log2FC + SE), na.rm = FALSE, orientation = NA, width = 0.2, color = "black")

      print(barplots)
      dev.off()
      

}




##########selected barplots
dir.create(paste0(dir,"GO_barplots/"))
importante_proteins <- read_excel("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/proteins_for_barplots.xlsx")
importante_proteins <- importante_proteins[rowSums(is.na(importante_proteins)) != ncol(importante_proteins), ]
my_txdf$N <- 1:nrow(my_txdf)

importante_proteins <- merge(my_txdf[,c("OldUserProtein","Gene","BLAST","Protein","ID","Product","N")], importante_proteins[,c("Group","Protein","N")], by = "N", all.y = TRUE)

#for (j in 1:nrow(importante_proteins)) {
#  if (importante_proteins$X.1[j] == "") {
#    importante_proteins$OldUserProtein[j] <- importante_proteins$gene[j]
#  }
#  else {
#    importante_proteins$OldUserProtein[j] <- importante_proteins$X.1[j]
#  }
#  if (j>=74) {
#    importante_proteins$OldUserProtein[j] <- importante_proteins$go_component[j]
#  }
#}
#importante_proteins <- importante_proteins[!duplicated(importante_proteins$gene),]

go_names <- unique(importante_proteins$Group)
go_names  <- na.omit(go_names)
df_go <- merge(df[,c(paste0("Mh",names,"/H34"),"ID")], importante_proteins, by = "ID", all.y = TRUE)
data_long <- df_go %>% gather(key = "variable", value = "value", -all_of(c("ID", "Gene", "Product", "N", "Group","OldUserProtein","BLAST","Protein.x", "Protein.y")))
#игнорирование выброса
#data_long[263,"value"] <- NA

for (go in go_names) {

  # go_tittle <- gsub("[[:space:]]", "_", go)
  # go_tittle <- strsplit(go_tittle, "|", fixed = TRUE)[[1]]
  # go_tittle <- go_tittle[1]
  # go_tittle_graph <- gsub("_", " ", go_tittle)
  go_tittle_graph <- go
  
  png(paste0(dir,"GO_barplots/GO_barplots_",go,".png"),width=3000,height=1000,res=200)

  filtered_data <- data_long %>%
    filter(Group == go)

  filtered_data <- filtered_data %>%
    group_by(variable) %>%
    mutate(x_intercept = row_number() + 0.5)

  num_colors <- length(unique(data_long$variable))
  color_palette <- rainbow(num_colors)

  go_barplots<- ggplot(filtered_data, aes(x = Protein.y, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab(expression("Proteins")) +
    ylab(expression("log"[2]*"FC")) +
    theme_bw() +
    geom_vline(aes(xintercept = x_intercept), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 0), color = "black", size = 0.3) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.position = "top",
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          axis.text.x = element_text(angle = 90, size = 7)) +
    scale_x_discrete(labels = filtered_data$Protein.y) +
    ggtitle(go_tittle_graph) +
    scale_fill_manual(values = color_palette) 
    #+ geom_errorbar(aes(ymin = value - SE, ymax = value + SE), na.rm = FALSE, orientation = NA, width = 0.2, color = "black")

  print(go_barplots)
  dev.off()

}


##########PCA
df_PCA <- df_PCA[complete.cases(df_PCA$ID), ]
rownames(df_PCA) <- df_PCA$ID
df_PCA$ID <- NULL
samples <- colnames(df_PCA)
isol <- rep(c("H-34", names), each = 3)


#deleted - common proteins
df_PCA_deleted <- df_PCA[complete.cases(df_PCA), ]
#df_PCA[is.na(df_PCA )] <- 0

#df_PCA<-apply(df_PCA,1,as.numeric)
df_PCA_deleted<-apply(df_PCA_deleted,1,as.numeric)


#scale - norm
res.pca = PCA(scale(df_PCA_deleted), ncp=5,  graph=F)
pca<-as.data.frame(res.pca$var$coord)
pca$Samples<-samples
pca$Isol<-isol

v<-round(as.numeric(unlist(res.pca$eig[,1]/sum(as.numeric(unlist(res.pca$eig[,1]))))*100),1)


png(paste0(dir,"PCA/PCA_deleted_norm.png"),width=4000,height=3000,res=200)

pca_plot<-ggplot(pca, aes(x=Dim.1, y=Dim.2))  +
  geom_point(aes(colour=factor(Isol)), size=10) +
  ggtitle("") + xlab(paste0('PC1 (',v[1],'%)')) + ylab(paste0("PC2 (",v[2],"%)")) +
  geom_text(label=pca$Samples, size=4, color='white', fontface="bold") + theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "top",
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(colour=NA))

print(pca_plot)
dev.off()




                      
