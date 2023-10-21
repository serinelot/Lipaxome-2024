setwd("C:/Users/Serine/Documents/1- Globus")

file_quantif = list.files("C:/Users/Serine/Documents/1- Globus", pattern = ".tsv$")

data_merge <- list()

for (i in (1:length(file_quantif))) {
  
  file_name = file_quantif[i]
  data = read.table(file_name, skip = 1, header = TRUE, sep = '\t')
  data <- data[, -c(2,3,4,5,6)]
  
  column_name <- colnames(data)
  bam_column <- column_name[2]
  
  srr_id_list <- strsplit(bam_column, split = ".", fixed = TRUE)
  srr_id_str <- unlist(srr_id_list)
  srr_id <- srr_id_str[3]
  names(data)[names(data) == bam_column] <- srr_id
  
  data_merge[[i]] <- data
}

data_merge <- Reduce(x = data_merge, f = function(x,y) {
  merge(x,y, by = "Geneid", sort = F)
})

write.csv2(data_merge, file = "data_merge.csv", row.names = F)
