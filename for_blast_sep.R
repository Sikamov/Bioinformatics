# Чтение файла FASTA
fasta <- readAAStringSet("/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Mh-H34_filt_bold.faa")

# Вычислите границы для разделения на три части
total_length <- length(fasta)
one_third <- total_length %/% 3

# Первая часть
fasta_first_part <- fasta[1:one_third]

# Вторая часть
fasta_second_part <- fasta[(one_third + 1):(2 * one_third)]

# Третья часть
fasta_third_part <- fasta[(2 * one_third + 1):total_length]

# Сохраните все три части в отдельные файлы
writeXStringSet(fasta_first_part, "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Mh-H34_filt_bold_first_part.faa")
writeXStringSet(fasta_second_part, "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Mh-H34_filt_bold_second_part.faa")
writeXStringSet(fasta_third_part, "/Users/kirillsikamov/Desktop/МФТИ/Биоинформатика ФХМ/micoplasma/references/bakta_MHO/Mh-H34_filt_bold_third_part.faa")