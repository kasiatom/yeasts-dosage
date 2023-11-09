require(readr)
require(dplyr)
require(tidyr)
require(stringr)



data <- read_delim("1011Matrix_vep_splited2.tsv",
  delim = "\t", escape_double = FALSE,
  na = "NA", trim_ws = TRUE)

 data <- data %>% 
   rename("pos" = POS, "ref" = REF, "alt" = ALT) %>%
   mutate("IS_PTV" = (VEP_BIOTYPE == "protein_coding" & grepl("HIGH", VEP_IMPACT)))
 
selected_strains <- read_delim("selected-ids-dendogram.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE)

selected_strains <- selected_strains$X1 ## all if do not filter

 ptv <- data %>%
  filter(IS_PTV) %>%
  pivot_longer(cols = matches('^[A-Z][A-Z][A-Z]$', ignore.case = FALSE), names_to = "STRAINS", values_to = "GENOTYPE") %>%
  mutate("IS_HET" = GENOTYPE %in% c("0/1", "1|0", "0|1")) %>%
   mutate("IS_ALT" = grepl("1", GENOTYPE)) %>%
   filter(STRAINS %in% selected_strains)
 
 summary_per_variant <- ptv %>%
   filter(IS_HET) %>%
   group_by(CHROM, pos, ref, alt) %>%
   add_count() %>%
   select(CHROM, pos, ref, alt, VEP_WORST_Gene, VEP_WORST_SYMBOL, VEP_WORST_Consequence, STRAINS, n) %>%
   mutate(STRAINS = paste(STRAINS, collapse=";")) %>%
   ungroup() %>%
   unique() %>%
   rename("HET_PTV_N" = n, "HET_PTV_STRAINS" = STRAINS) %>%
   as.data.frame()
 
 write_tsv(summary_per_variant, "summary_per_variant_selected_strains_dend.tsv")
 
 summary_per_gene <- ptv %>%
   filter(IS_HET) %>%
   rowwise() %>%
   mutate("VAR_ID" = paste(CHROM, "-", pos, "-", ref, "-", alt, sep="")) %>%
   ungroup() %>%
   select(VAR_ID, VEP_WORST_Gene, VEP_WORST_SYMBOL, STRAINS) %>%
   group_by(VEP_WORST_Gene) %>%
   mutate(STRAINS = paste(unique(STRAINS), collapse=";"), VAR_ID= paste(unique(VAR_ID), collapse=";")) %>%
   unique() %>%
   mutate("HET_PTV_N" = str_count(STRAINS, ";")+1) %>%
   ungroup() %>%
   unique() %>%
   rename( "HET_PTV_STRAINS" = STRAINS, "HET_PTV_VAR_IDS" = VAR_ID) %>%
   as.data.frame()
 
 write_tsv(summary_per_gene, "summary_per_gene_selected_strains_dend.tsv")
 
 summary_per_strain <- ptv %>%
   filter(IS_HET) %>%
   select(STRAINS, VEP_WORST_Gene, VEP_WORST_SYMBOL) %>%
   unique() %>%
   group_by(STRAINS) %>%
   add_count() %>%
   mutate("HET_PTV_GENES" = paste(unique(VEP_WORST_Gene), collapse=";"),"HET_PTV_NAMES" = paste(unique(VEP_WORST_SYMBOL), collapse=";")) %>%
   select(-VEP_WORST_Gene, -VEP_WORST_SYMBOL) %>%
   ungroup() %>%
   unique() %>%
   as.data.frame()
 
 write_tsv(summary_per_strain, "summary_per_strain_selected_strains_dend.tsv")
 