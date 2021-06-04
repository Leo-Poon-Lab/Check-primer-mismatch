library(Biostrings)
library(tidyverse)

data <- read_tsv("../data/primer_probe_data.tsv")
seqs <- readDNAStringSet("../../2020-09-01_COVID_NGS_pipeline/results/latest_genome_caseid.fasta")
ref_seq <- readDNAStringSet("../../2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")

df_check <- lapply(seq_len(nrow(data)), function(i){
    print(i)
    data_i <- data[i,]
    name_i <- data_i$name
    start_i <- data_i$start
    stop_i <- data_i$stop
    seq_sub <- subseq(seqs, start_i, stop_i)
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
        paste(rst, collapse = ", ")
    })
    df_out$mut_seqs <- mut_seqs
    return(df_out)
})

df_check <- bind_rows(df_check)
write_csv(df_check, "../resutls/df_check.csv")
