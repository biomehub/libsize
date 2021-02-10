source("code/R/general_functions.R")
library(patchwork)
theme_set(theme_pubr())
check_dir("output/supp_figures")


# Supp Fig 1

phylo <- load_phyloseq(
  otutable = "data/raw/extraction_method/200117-142337/otu_table_count_200117-142337.tsv",
  taxtable = "data/raw/extraction_method/200117-142337/tax_table_200117-142337.tsv",
  metadata = "data/raw/extraction_method/200117-142337/200117-142337_metadata.tsv",
  remove_zeros = FALSE
) %>%
  phylo_lca()

df <- tidy_phylo_data(phylo)
p <- df %>%
  mutate(extraction_method = factor(
    as.character(extraction_method),
    levels = c("Beads", "PowerSoil",
               "QIAMP", "PowerSoilPRO")
  )) %>%
  ggplot(aes(cfu, .value + 1, fill = extraction_method)) +
  geom_point(pch = 21, size = 3) +
  yscale('log10', TRUE) +
  xscale('log10', TRUE) +
  facet_wrap(~extraction_method) +
  labs(x = "CFU", y = "# reads + 1", fill = NULL) +
  theme(strip.text = element_text(face = "bold", color = "gray20", size = 14),
        axis.title = element_text(face = "bold", color = "gray20", size = 12)) +
  geom_smooth(method = 'gam',
              aes(group = extraction_method,
                  fill = NULL),
              formula = y ~ splines::bs(x),
              show.legend = FALSE)

ggsave("output/supp_figures/Supp_fig_1.png", p,
        width = 10, height = 6.5, dpi = 400)


