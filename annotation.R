######Vien diagrams
ann_bakta <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Mh-H34.gff3"
ann_bakta_NCBI <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Mh-H34_NCBI_ATCC.tsv"
ann_OV <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/аннотация к последней фасте.txt"
ann_NCBI <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/GCF_000759385.1_ASM75938v1_genomic.gtf"
ann_Uniprot <-  "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Uniprot_ATCC.txt"
ann_NCBI_Uniprot <-  "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/NCBI_Uniprot.txt"
ann_ATCC_23114 <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/GCA_000085865.1_ASM8586v1_genomic.gff"
ann_bakta_PCM <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/Mh-H34_PCM_ATCC23114.tsv"
ann_pgap <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/annot.gff"
ann_MHO <- "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/MHO_AU.tsv"

gff_bakta <- import.gff(ann_bakta, format = "gff3")
gff_bakta <- as.data.frame(gff_bakta)
gff_bakta$Note <- as.character(gff_bakta$Note)
#gff <- subset(gff,gff$product!="hypothetical protein")
db_bakta <- as.data.frame(gff_bakta[,c("locus_tag","gene","product","Note")])
db_bakta  <- db_bakta [!duplicated(db_bakta ),]
colnames(db_bakta)[colnames(db_bakta) == "locus_tag"] <- "ID"
db_bakta <- db_bakta[complete.cases(db_bakta$product), ]
db_bakta <- db_bakta[complete.cases(db_bakta$ID), ]
db_bakta <- db_bakta[complete.cases(db_bakta$gene), ]


gff_NCBI <- import.gff(ann_NCBI, format = "gtf")
gff_NCBI <- as.data.frame(gff_NCBI)
#gff <- subset(gff,gff$product!="hypothetical protein")
db_NCBI <- as.data.frame(gff_NCBI[,c("locus_tag","gbkey","gene","product","protein_id")])
db_NCBI  <- db_NCBI [!duplicated(db_NCBI ),]
colnames(db_NCBI)[colnames(db_NCBI) == "protein_id"] <- "ID"
db_NCBI <- db_NCBI[complete.cases(db_NCBI$ID), ]
db_NCBI <- db_NCBI[complete.cases(db_NCBI$gene), ]

gff_ATCC <- import.gff(ann_ATCC_23114, format = "gff")
#gff <- subset(gff,gff$product!="hypothetical protein")
db_ATCC <- as.data.frame(gff_ATCC)

db_ATCC$UniProtKB <- sapply(db_ATCC$Dbxref, function(x) extract_value(x, "UniProtKB/TrEMBL:"))
db_ATCC$NCBI_GP <- sapply(db_ATCC$Dbxref, function(x) extract_value(x, "NCBI_GP:"))


colnames(db_ATCC)[colnames(db_ATCC) == "locus_tag"] <- "OldUserProtein"
db_ATCC  <- gff_ATCC [!duplicated(gff_ATCC ),]
colnames(db_ATCC)[colnames(db_NCBI) == "protein_id"] <- "ID"
db_ATCC <- db_NCBI[complete.cases(db_NCBI$ID), ]
db_ATCC <- db_NCBI[complete.cases(db_NCBI$gene), ]


gff_pgap <- import.gff(ann_pgap, format = "gff")
#gff <- subset(gff,gff$product!="hypothetical protein")
db_pgap <- as.data.frame(gff_pgap)
colnames(db_pgap)[colnames(db_pgap) == "locus_tag"] <- "OldUserProtein"
db_pgap  <- gff_pgap [!duplicated(gff_pgap ),]
colnames(db_ATCC)[colnames(db_NCBI) == "protein_id"] <- "ID"
db_ATCC <- db_NCBI[complete.cases(db_NCBI$ID), ]
db_ATCC <- db_NCBI[complete.cases(db_NCBI$gene), ]



db_uniprot <- read.table(ann_Uniprot, header = TRUE, quote="", sep = "\t")
max_words <- max(sapply(db_uniprot$Gene.Names, length))
db_uniprot <- db_uniprot %>%
  separate(Gene.Names, into = "gene", sep = "\\s+", fill = "right")
# db_uniprot <- as.data.frame(db_uniprot[,c("locus_tag","ftype", "gene","product")])
# db_uniprot <- db_uniprot[!duplicated(db_uniprot),]
# db_uniprot <- db_uniprot[complete.cases(db_uniprot$locus_tag), ]
# colnames(db_OV)[colnames(db_uniprot) == "locus_tag"] <- "ID"
# db_uniprot <- db_uniprot[complete.cases(db_uniprot$gene), ]

db_MHO <- read.table(ann_MHO, header = TRUE, quote="", sep = "\t")
db_MHO$RefSeq <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "RefSeq:"))
db_MHO$SO <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "SO:"))
db_MHO$UniParc <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "UniParc:"))
db_MHO$UniRef <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "UniRef:"))
db_MHO$UserProtein <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "UserProtein:"))
db_MHO$GO <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "GO:"))
db_MHO$KEGG <- sapply(db_MHO$DbXrefs, function(x) extract_value(x, "KEGG:"))

db_OV <- read.table(ann_OV, header = TRUE, quote="", sep = "\t")
db_OV <- as.data.frame(db_OV[,c("locus_tag","ftype", "gene","product")])
db_OV <- db_OV[!duplicated(db_OV),]
db_OV <- db_OV[complete.cases(db_OV$locus_tag), ]
colnames(db_OV)[colnames(db_OV) == "locus_tag"] <- "ID"
db_OV <- db_OV[complete.cases(db_OV$gene), ]




db_bakta_NCBI <- read.table(ann_bakta_PCM, header = TRUE, quote="", sep = "\t")

# Функция для извлечения значения после указанного префикса
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
# Создание новых колонок в db_bakta_NCBI
db_bakta_NCBI$RefSeq <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "RefSeq:"))
db_bakta_NCBI$SO <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "SO:"))
db_bakta_NCBI$UniParc <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "UniParc:"))
db_bakta_NCBI$UniRef <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "UniRef:"))
db_bakta_NCBI$UserProtein <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "UserProtein:"))
db_bakta_NCBI$GO <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "GO:"))
db_bakta_NCBI$KEGG <- sapply(db_bakta_NCBI$DbXrefs, function(x) extract_value(x, "KEGG:"))


###################################################ols<-new

lines <- readLines("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/GCF_000085865.1_ASM8586v1_genomic.faa")
locus_tag <- vector()
old_locus_tag <- vector()

for (i in 1:length(lines)) {
  if (grepl("/locus_tag", lines[i])) {
    locus_tag_value <- sub('^.*"([^"]+)".*', '\\1', lines[i])
    locus_tag <- c(locus_tag, locus_tag_value)
    if (grepl("/old_locus_tag", lines[i+1])) {
      old_locus_tag_value <- sub('^.*"([^"]+)".*', '\\1', lines[i+1])
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



#colnames(db_bakta_NCBI)[colnames(db_bakta_NCBI) == "RefSeq"] <- "From"
db_bakta_NCBI <- merge(db_bakta_NCBI, old_new, by = "UserProtein", all.x = TRUE)
db_bakta_NCBI <- merge(db_bakta_NCBI, blast, by = "Locus.Tag", all.x = TRUE)
colnames(db_bakta_NCBI)[colnames(db_bakta_NCBI) == "old_locus_tag"] <- "OldUserProtein"
#db_bakta_NCBI <- as.data.frame(db_bakta_NCBI[,c("From", "gene","Gene","Product","UserProtein","Gene.Names..ordered.locus.")])


for (j in 1:nrow(db_bakta_NCBI)) {
  if (is.na(db_bakta_NCBI$Gene[j])||db_bakta_NCBI$Gene[j] == "") {
    if (grepl("Hypothetical protein", db_bakta_NCBI$Product[j])||grepl("hypothetical protein", db_bakta_NCBI$Product[j])||is.na(db_bakta_NCBI$Product[j])||db_bakta_NCBI$Product[j] == "") {
      if (grepl("Hypothetical protein", db_bakta_NCBI$BLAST[j])||grepl("hypothetical protein", db_bakta_NCBI$BLAST[j])||is.na(db_bakta_NCBI$BLAST[j])||db_bakta_NCBI$BLAST[j] == "") {
        if (is.na(db_bakta_NCBI$OldUserProtein[j])||db_bakta_NCBI$OldUserProtein[j] == "") {
          db_bakta_NCBI$Protein[j] <- db_bakta_NCBI$Product[j]
        } else { db_bakta_NCBI$Protein[j] <- db_bakta_NCBI$OldUserProtein[j]
          }
      } else {
        db_bakta_NCBI$Protein[j] <- db_bakta_NCBI$BLAST[j]
        }
    } else {
      db_bakta_NCBI$Protein[j] <- db_bakta_NCBI$Product[j]
    }
  } else {
    db_bakta_NCBI$Protein[j] <- db_bakta_NCBI$Gene[j]
  }
}

# db_bakta_NCBI <- db_bakta_NCBI[!duplicated(db_bakta_NCBI),]
# db_bakta_NCBI <- db_bakta_NCBI[complete.cases(db_bakta_NCBI$ID), ]
# db_bakta_NCBI <- db_bakta_NCBI[complete.cases(db_bakta_NCBI$Gene), ]
write.table(db_bakta_NCBI$RefSeq, file = "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/ID_.txt", row.names = FALSE, col.names = FALSE)

#######################################################################################



db <- read.table(ann_NCBI_Uniprot, header = TRUE, quote="", sep = "\t")
for (j in 1:nrow(db)) {
  if (db$Gene.Names..primary.[j] == "") {
    if (db$Gene.Names..ordered.locus.[j] == "") {
      db$gene[j] <- db$Protein.names[j]
    } else {
      db$gene[j] <- db$Gene.Names..ordered.locus.[j]
    }
  } else {
    db$gene[j] <- db$Gene.Names..primary.[j]
  }
}
db <- as.data.frame(db[,c("From","Entry","gene","Protein.names", "Gene.Names..primary.","Gene.Names..ordered.locus.","Gene.Ontology..cellular.component.","Gene.Ontology..biological.process.")])


# colnames(db_uniprot)[colnames(db_uniprot) == "gene"] <- "gene_uniprot"
# db <- merge(db, db_uniprot[,c("gene_uniprot")], by = "Entry", all.x = TRUE)
# colnames(db_uniprot)[colnames(db_uniprot) == "gene_uniprot"] <- "gene"

importante_proteins <- read.table("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/важные белки.csv", header = TRUE, sep = ";")
for (j in 1:nrow(importante_proteins)) {
  if (importante_proteins$X.1[j] == "") {
    importante_proteins$OldUserProtein[j] <- importante_proteins$gene[j]
  }
  else {
    importante_proteins$OldUserProtein[j] <- importante_proteins$X.1[j]
  }
  if (j>=74) {
    importante_proteins$OldUserProtein[j] <- importante_proteins$go_component[j]
  }
}
importante_proteins <- importante_proteins[!duplicated(importante_proteins$gene),]

write.table(db_MHO$RefSeq, file = "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/RefSeq_MHO.txt", row.names = FALSE, col.names = FALSE)

selected_ID <- read.table("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/ID_selected.txt", header = TRUE, sep = "\t")





list_df <- list(selected_ID$Entry, MHO_entry$Entry)



list_df <- lapply(list_df, function(df) {
  df[is.na(df)] <- " "
  return(df)
})

dir.create("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/venn_diagram/")
venn.diagram(
  x = list_df,
  category.names = c("selected proteins","MHO"),
  filename = "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/venn_diagram/venn_diagram_MHO_Entry.png",
  output = TRUE
)










########################################################annotations old

#########annotation GTF
GTF <- import.gff(annotation, format = "gff3")
gff <- as.data.frame(GTF[,c("locus_tag","gene","Name")])
#gff <- subset(gff,gff$product!="hypothetical protein")
my_txdf <- as.data.frame(gff[,c("locus_tag","gene","Name")])
my_txdf <- my_txdf[!duplicated(my_txdf),]
colnames(my_txdf)[colnames(my_txdf) == "locus_tag"] <- "ID"
my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]

for (j in 1:nrow(my_txdf)) {
  if (is.na(my_txdf$gene[j])) {
    if (my_txdf$Name[j] == "hypothetical protein") {
      my_txdf$Protein[j] <- my_txdf$ID[j]
    } else {
      my_txdf$Protein[j] <- my_txdf$Name[j]
    }
  } else {
    my_txdf$Protein[j] <- my_txdf$gene[j]
  }
}



diff <- merge(diff, my_txdf, by = "ID", all.x = TRUE)
data_all <- merge(data_all, my_txdf, by = "ID", all.x = TRUE)



#########annotation txt
#    my_txdf <- read.table(file=annotation", header = TRUE, quote="", sep = "\t")
#    my_txdf <- my_txdf[!duplicated(my_txdf),]
#    my_txdf <- my_txdf[complete.cases(my_txdf$locus_tag), ]
#    colnames(my_txdf)[colnames(my_txdf) == "locus_tag"] <- "Protein"
#    diff <- merge(diff, my_txdf, by = "Protein", all.x = TRUE)
#    data_all <- merge(data_all, my_txdf, by = "Protein", all.x = TRUE)


######annotations_ATCC
my_txdf <- db_ATCC
colnames(my_txdf)[colnames(my_txdf) == "ID"] <- "name"
colnames(my_txdf)[colnames(my_txdf) == "seqnames"] <- "ID"
my_txdf <- as.data.frame(my_txdf[,c("ID","locus_tag","gene","product","Name")])
my_txdf <- my_txdf[!duplicated(my_txdf),]
colnames(my_txdf)[colnames(my_txdf) == "old_locus_tag"] <- "OldUserProtein"
my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
#my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]

for (j in 1:nrow(my_txdf)) {
  if (is.na(my_txdf$gene[j])||my_txdf$gene[j] == "") {
    if (grepl("Hypothetical protein", my_txdf$product[j])||grepl("hypothetical protein", my_txdf$product[j])||is.na(my_txdf$product[j])) {
      my_txdf$Protein[j] <- my_txdf$locus_tag[j]
    } else {
      my_txdf$Protein[j] <- my_txdf$product[j]
    }
  } else {
    my_txdf$Protein[j] <- my_txdf$gene[j]
  }
}


######annotations_ATCC_uniprot

#######data_uniprot
data_all <- read.csv(paste0(dir, "mycoplasma_design_Mh_", name, "_openms_design_msstats_in_comparisons.csv")
                     , sep = "\t")
data_all <- data_all %>%
  separate(Protein, into = c("var1", "ID", "var3"), sep = "\\|", remove = FALSE) %>%
  dplyr::select(-Protein)


my_txdf <- db_uniprot
colnames(my_txdf)[colnames(my_txdf) == "Entry"] <- "ID"
my_txdf <- as.data.frame(my_txdf[,c("ID","Protein.names","gene")])
my_txdf <- my_txdf[!duplicated(my_txdf),]
colnames(my_txdf)[colnames(my_txdf) == "Protein.names"] <- "product"
my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
#my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]




############BAKTA_working_variant
######annotations
my_txdf <- db_bakta_NCBI
colnames(my_txdf)[colnames(my_txdf) == "Locus.Tag"] <- "ID"
my_txdf <- as.data.frame(my_txdf[,c("ID","UserProtein","Gene","Product","RefSeq","GO","KEGG","old_locus_tag")])
my_txdf <- my_txdf[!duplicated(my_txdf),]
colnames(my_txdf)[colnames(my_txdf) == "old_locus_tag"] <- "OldUserProtein"
my_txdf <- my_txdf[complete.cases(my_txdf$ID), ]
#my_txdf <- my_txdf[complete.cases(my_txdf$Name), ]

for (j in 1:nrow(my_txdf)) {
  if (is.na(my_txdf$Gene[j])||my_txdf$Gene[j] == "") {
    if (grepl("Hypothetical protein", my_txdf$Product[j])||grepl("hypothetical protein", my_txdf$Product[j])||is.na(my_txdf$Product[j])||my_txdf$Product[j] == "") {
      my_txdf$Protein[j] <- my_txdf$OldUserProtein[j]
    } else {
      my_txdf$Protein[j] <- my_txdf$Product[j]
    }
  } else {
    my_txdf$Protein[j] <- my_txdf$Gene[j]
  }
}


                      