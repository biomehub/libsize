
PROJECT_DIRECOTRY <- '~/workstation/projetos/lib_size'
setwd(PROJECT_DIRECOTRY)
source('code/R/general_functions.R')

norm_data_dir = 'data/normalized'
check_dir(norm_data_dir)

## calculate size factors for normalization
seqs_data <- read_tsv('data/raw/sequencing_metadata.tsv')


seqs_data <- seqs_data %>%
  # calculate nominal recovery (observed reads / available a priori)
  mutate(
    # calculate relative availability (sample covs a priori - availavle -. not equal to observed lib size)
    # size factor is the seq's 'proportion' to max relative availability among all seqs
    size_factor = Mean_Sample_Cov_Pool / max(Mean_Sample_Cov_Pool)
  ) %>%
  as.data.frame()

metadata <- read_tsv('data/raw/metadata.tsv')

metadata_updated <- left_join(metadata, seqs_data, by = 'seq_id') %>%
  select(-E_coli) %>%
  mutate(total_cfu = reduce(select(., L_monocytogenes:E_faecalis), `+`))

a = metadata_updated %>%
  select(seq_id, size_factor) %>%
  unique() %>%
  as.data.frame(row.names = 1:4)
b = seqs_data %>%
  select(seq_id, size_factor) %>%
  unique() %>%
  as.data.frame(row.names = 1:4)

test = all.equal(a, b)

stopifnot(test)

metadata_updated %>%
  write_tsv(
    path = "data/normalized/metadata.tsv"
  )

## merge data and normalize

data_list <- list(
  'seq_1' = list(
    'otu' = 'data/raw/otu_table_count_190728-202300_neorefdb16sv6.tsv',
    'tax' = 'data/raw/tax_table_190728-202300_neorefdb16sv6.tsv',
    'seq_depth' = 4178618
  ),
  'seq_2' = list(
    'otu' = 'data/raw/otu_table_count_190829-003158_neorefdb16sv6.tsv',
    'tax' = 'data/raw/tax_table_190829-003158_neorefdb16sv6.tsv',
    'seq_depth' = 3179771
  ),
  'seq_3' = list(
    'otu' = 'data/raw/otu_table_count_190829-003329_neorefdb16sv6.tsv',
    'tax' = 'data/raw/tax_table_190829-003329_neorefdb16sv6.tsv',
    'seq_depth' = 2801656
  ),
  'seq_4' = list(
    'otu' = 'data/raw/otu_table_count_190930-091404_neorefdb16sv6.tsv',
    'tax' = 'data/raw/tax_table_190930-091404_neorefdb16sv6.tsv',
    'seq_depth' = 1472198
  )
)


phylo_list <- map(data_list, function(data_pars) {
  metadata <- 'data/normalized/metadata.tsv'
  phylo <- load_phyloseq(
    otutable = data_pars$otu,
    taxtable = data_pars$tax,
    metadata = metadata,
    meta_smpls_only = TRUE,
    remove_zeros = FALSE,
    merge_dilutions = FALSE
  ) %>% phylo_lca(.glom = TRUE)
  phylo@sam_data$seq_depth <- data_pars$seq_depth
  return(phylo)
})

phylo_merged <- merge_phyloseq(
  phylo_list$seq_1,
  phylo_list$seq_2,
  phylo_list$seq_3,
  phylo_list$seq_4
)
taxa_to_include <- c(
  "Listeria monocytogenes", "Salmonella enterica", "Bacillus cereus",
  "Staphylococcus aureus", "Staphylococcus epidermidis", "Enterococcus faecalis"
)

phylo_merged <- prune_taxa(taxa_names(phylo_merged) %in% taxa_to_include, phylo_merged)
# phylo_merged <- prune_samples(sample_sums(phylo_merged) > 0, phylo_merged)

saveRDS(phylo_merged, 'data/raw/phyloseq_raw_data.rds')
phylo_merged@sam_data$lib_size_raw <- sample_sums(phylo_merged) 
smpl_data <- phylo_merged@sam_data %>% data.frame()
otu_data <- phylo_merged@otu_table@.Data
otu_normalized <- otu_data / smpl_data$size_factor
smpl_data$lib_size <- rowSums(otu_normalized)[rownames(smpl_data)] 
sample_data(phylo_merged) <- smpl_data

otu_table(phylo_merged) <- otu_table(
            otu_normalized,
            taxa_are_rows = FALSE
          )

saveRDS(
  phylo_merged, 
  file = str_glue('{norm_data_dir}/phyloseq.rds')
)

otu <- otu_table(phylo_merged) %>% data.frame()
tax <- tax_table(phylo_merged) %>% data.frame()
smpl <- sample_data(phylo_merged) %>% data.frame()

write_tsv(otu, str_glue('{norm_data_dir}/otu_table.tsv'))
write_tsv(tax, str_glue('{norm_data_dir}/tax_table.tsv'))
write_tsv(smpl, str_glue('{norm_data_dir}/sample_data.tsv'))
