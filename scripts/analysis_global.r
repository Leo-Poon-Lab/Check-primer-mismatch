library(Biostrings)
library(tidyverse)

data <- read_tsv("../data/primer_probe_data.tsv")
seqs <- readDNAStringSet("../data/mmsa_2021-06-03/2021-06-03_masked.fa")
ref_seq <- readDNAStringSet("../../2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
df_meta <- read_tsv("../data/metadata.tsv", col_types = cols(.default = "c"))

# seqs in different continent, especially Africa and Europe
# since 2021
Locations <- c("Africa", "Europe", "Asia", "North America", "Oceania", "South America")
df_meta <- df_meta %>% filter(`Collection date`>= "2021")

df_check <- mclapply(Locations, function(location_t){
	print(location_t)
	df_t <- df_meta %>% filter(grepl(location_t, Location))
	seq_t <- seqs[names(seqs) %in% df_t$`Accession ID`]
	
	df_check <- lapply(seq_len(nrow(data)), function(i){
		print(i)
		data_i <- data[i,]
		name_i <- data_i$name
		start_i <- data_i$start
		stop_i <- data_i$stop
		seq_sub <- subseq(seq_t, start_i, stop_i)
		ref_seq_sub <- subseq(ref_seq, start_i, stop_i)
		con_m_i <- consensusMatrix(seq_sub)
		prop_wt <- sapply(seq_len(width(ref_seq_sub)), function(j){    
			num_j <- con_m_i[as.character(subseq(ref_seq_sub, j, j)),j]
			return(num_j/sum(con_m_i[1:4,j]))
		})
		df_out <- tibble(Name = name_i, Position = start_i:stop_i, Ref_base = names(prop_wt), mut_rate_precentage = round((1-prop_wt)*100, 5))
		mut_seqs <- sapply(seq_len(nrow(df_out)), function(ii){
			df_ii <- df_out[ii,]
			if(df_ii$mut_rate_precentage==0){return("")}
			rst <- names(seq_sub)[(subseq(seq_sub, ii, ii) != df_ii$Ref_base) & (subseq(seq_sub, ii, ii) %in% c("A", "C", "G", "T"))]
			paste(rst, collapse = "｜")
		})
		df_out$mut_seqs <- mut_seqs
		return(df_out)
	})

	df_check <- bind_rows(df_check)
	df_check$Location <- location_t
	write_csv(df_check %>% arrange(desc(mut_rate_precentage)), paste0("../resutls/df_check_", location_t, ".csv"))
}, mc.cores = 6)

df_check <- bind_rows(df_check)
df_check_n_probe <- df_check %>% filter(Name == "Probe_N") %>% arrange(desc(mut_rate_precentage)) %>% select(Name:mut_rate_precentage, Location)
write_excel_csv(df_check_n_probe, "../resutls/df_check_n_probe.csv")
