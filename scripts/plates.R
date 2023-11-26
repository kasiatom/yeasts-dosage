require(dplyr)
require(readr)

strains <-
  read_delim(
    "strain-info.tsv",
    delim = "\t",
    col_names = TRUE,
    trim_ws = TRUE
  )
strains <- strains %>%
  select(plate_position, standard_name) %>%
  mutate(plate_position = gsub("__", "_", plate_position))

chosen_strains <-
  read_delim(
    "selected-ids-dendogram.txt",
    delim = "\t",
    col_names = FALSE,
    trim_ws = TRUE
  )
colnames(chosen_strains) <- "standard_name"


chosen_strains <- chosen_strains %>%
  left_join(., strains, by = "standard_name") %>%
  mutate("plate" = as.numeric(gsub("_.*", "", plate_position))) %>%
  mutate("well" = gsub("^[0-9]*_", "", plate_position)) %>%
  mutate("row" = gsub("[0-9]*", "", well)) %>%
  mutate("column" = as.numeric(gsub("[A-Z]*", "", well))) %>%
  arrange(plate, column, row) %>%
  select(standard_name, plate_position)

plates <- c(rep("1", 96), rep("2", 8))
rows <- rep(c("A", "B", "C", "D", "E", "F", "G", "H"), 13)
columns <-
  c(
    rep("1", 8),
    rep("2", 8),
    rep("3", 8),
    rep("4", 8),
    rep("5", 8),
    rep("6", 8),
    rep("7", 8),
    rep("8", 8),
    rep("9", 8),
    rep("10", 8),
    rep("11", 8),
    rep("12", 8),
    rep("1", 8)
  )
plates <- cbind(plates, rows, columns) %>%
  head(98) %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(plates = paste(plates, "_", sep=""))%>%
  mutate("new_plate_position" = paste0(plates, rows, columns, collapse ="")) %>%
  select(new_plate_position)

chosen_strains <- cbind(chosen_strains, plates)
colnames(chosen_strains) <- c("strain", "original_plate_position", "new_plate-position")

write_tsv(chosen_strains, "selected-strains-dendogram-plate-info.tsv")

