samples <- c("32", "29", "30", "34", "33", "28", "35", "114","36")

samples_inf <- data.frame(Sample = c("32", "29", "30", "34", "33", "28", "35", "114","36"), Isol = c("33", "43", "11", "1862", "12", "7", "45", "40","H34"))
list_df <- list()
for (sample in samples) {
  annotation <- paste0("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/",sample,".tsv")
  
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

  
  
  write.table(as.data.frame(my_txdf), file = paste0("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Annotation_",samples_inf$Isol[samples_inf$Sample == sample],".tsv"), row.names = TRUE, col.names = TRUE)
  list_df <- append(list_df, list(my_txdf$Protein))
}

names(list_df) <- c("33", "43", "11", "1862", "12", "7", "45", "40", "H34")
png("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Upsetplot_1.png",width=2000,height=1000,res=200)
upset(fromList(list_df), intersect = c("33", "43", "11", "1862", "12", "7", "45", "40", "H34"), set_sizes  = FALSE, min_size=3)
dev.off()


