# Чтение данных из файла
data_blast <- readLines("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/H34_vs_MHO/H34_vs_ATCC23114.txt")

# Инициализация пустых списков для хранения значений
query_id <- vector("character")
description <- vector("character")

# Обработка данных
for (g in 1:(length(data_blast) - 1)) {
  if (startsWith(data_blast[g], "Query #")) {
    query_id_val <- gsub("Query #[0-9]+: (.+?) .+", "\\1", data_blast[g])
    if (grepl("No significant similarity found.", data_blast[g + 2])) {
    next
    }
  
  } else if (startsWith(data_blast[g], "Description")) {
    
    if (grepl("Metamycoplasma hominis ATCC 23114", data_blast[g + 1]) && grepl("100%", data_blast[g + 1])) {
    description_val <- gsub("\\[Meta.*|Meta.*", "", data_blast[g + 1])
    
    } else if (grepl("Metamycoplasma hominis ATCC 23114", data_blast[g + 2]) && grepl("100%", data_blast[g + 2])) {
      description_val <- gsub("\\[Meta.*|Meta.*", "", data_blast[g + 2])
      
      } else if (grepl("Metamycoplasma hominis ATCC 23114", data_blast[g + 3]) && grepl("100%", data_blast[g + 3])) {
        description_val <- gsub("\\[Meta.*|Meta.*", "", data_blast[g + 3])
        
        } else {
          description_val <- gsub("\\[Meta.*|Meta.*", "", data_blast[g + 1])
        }
    
    query_id <- c(query_id, query_id_val)
    description <- c(description, description_val)
  }
}

# db_bakta_NCBI <- merge(db_bakta_NCBI, old_new, by = "UserProtein", all.x = TRUE)
# db_bakta_NCBI <- merge(db_bakta_NCBI, blast, by = "Locus.Tag", all.x = TRUE)

# Создание DataFrame
blast <- data.frame(Locus.Tag = query_id, BLAST = description)
