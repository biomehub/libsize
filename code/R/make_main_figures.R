source("code/R/general_functions.R")
library(patchwork)
theme_set(theme_pubr())

#
# This script mostly joins the figures obtained in the supplementary materials. 
# Specific dependencies and guides to source code follow load/read/import commands.
#

figs_dir <- "output/main_text_figures/final_figures/"
check_dir(figs_dir)

# Fig 1:

"
Figure 1 was made by hand :)
"

# Fig 2
set.seed(13528)
fig2a <- readr::read_tsv("data/raw/dna_concentrations/dna_concentration_data.tsv") %>%
  rename(dna_concentration := `dna_concentration(ng/uL)`) %>% 
  ggplot(aes(dna_concentration, ngs_reads + 1)) +
  geom_boxplot(aes(group = factor(dna_concentration)), color = 'gray40', alpha = 0) + 
  geom_jitter(height = 0, width = .2, size = 2) +
  geom_smooth(method = 'lm', formula = y ~ poly(x, 3)) +
  yscale("log10", .format = TRUE) +
  xscale("log10", .format = TRUE) +
  labs(
    x = expression(bold(paste("DNA concentration (",n,g,"/",mu,l, ")"))),
    y = "# reads"
  ) +
  theme(axis.title = element_text(face = "bold", color = "gray20", size = 12))
set.seed(7236)
fig2b <- readr::read_tsv("data/raw/dna_concentrations/dna_concentration_data.tsv") %>%
  ggplot(aes(dna_copies, ngs_reads + 1)) +
  geom_boxplot(aes(group = factor(dna_copies)), color = 'gray40', alpha = 0) + 
  geom_jitter(height = 0, width = .2, size = 2) +
  geom_smooth(method = 'lm', formula = y ~ poly(x, 3)) +
  yscale("log10", .format = TRUE) +
  xscale("log10", .format =  TRUE) +
  labs(
    x = "DNA copies",
    y = "# reads"
  ) +
  theme(axis.title = element_text(face = "bold", color = "gray20", size = 12))

### depend on Supp Material 2
fig2c <- readRDS("output/main_text_figures/fig2/fig2c.rds")
fig2d <- readRDS("output/main_text_figures/fig2/fig2d.rds")

### depends on Supp Material 3
fig2e <- readRDS("output/main_text_figures/fig2/fig2e.rds")
fig2f <- readRDS("output/main_text_figures/fig2/fig2f.rds")

figure2 <- ggarrange(fig2a, fig2b, fig2c, 
                     fig2d + guides(fill = 'none'),
                     fig2e, fig2f,
                     nrow = 3, ncol = 2, labels = "AUTO",
                     common.legend = TRUE,
                     legend = "bottom")

ggsave(str_glue("{figs_dir}/figure2.png"), figure2,
       width = 9, height = 9, dpi = 400)

# Fig 3

# from total microbial load model (supp mat 2, last code chunk, line ~ 1005)
fig3a <- readRDS("output/main_text_figures/fig3/posterior_probablities.rds") +
  theme(legend.key.size = unit(.25, 'cm'),
        legend.text = element_text(size = 8))

fig3b <- readRDS("output/main_text_figures/fig3/conditional_expectations.rds")

fig3c <- readRDS("output/main_text_figures/fig3/posterior_check_bars.rds") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10)) +
  scale_x_discrete(labels = paste0("0.84e+0", 2:6),
                   limits = 1:5)

fig3d <- readRDS("output/main_text_figures/fig3/tail_probs_fit1.rds") + 
  theme(legend.key.size = unit(.25, 'cm'),
        legend.text = element_text(size = 8),
        legend.position = "top") +
  guides(color = "none") +
  labs(y = "Tail probabilities")

# from mixed effects model (supp mat 2, last code chunk, line ~ 1185)
fig3e <- readRDS("output/main_text_figures/fig3/mixed_model_post_pred_check.rds") + 
  scale_x_discrete(labels = scales::scientific(2 * 10^(2:5)),
                   limits = 1:4) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        strip.text = element_text(size = 11, face = "bold.italic", color = "gray20"))

figure3 <- ggarrange(
  ggarrange(fig3a, fig3c, fig3b, fig3d, nrow = 2, ncol=2, labels = c("A", "C", "B", "D") ),
  ggarrange(fig3e, labels = "E"),
  nrow = 2
)

ggsave("output/main_text_figures/final_figures/figure3.png", figure3,
       height = 14, width = 12, dpi = 400)

# Fig 4

fig4a <- readRDS("output/main_text_figures/fig4/total_microbial_load_kfold.rds") +
  theme(axis.text.x = element_text(size = 8))
fig4b <- readRDS("output/main_text_figures/fig4/total_microb_load_test_set_performance.rds") +
  theme(axis.text.x = element_text(size = 8),
        legend.position = "right")
fig4c <- readRDS("output/main_text_figures/fig4/mixed_model_kfold.rds") +
  theme(axis.text = element_text(size = 9)) + guides(fill = 'none') +
  scale_y_continuous(breaks = scales::pretty_breaks(5))
fig4d <- readRDS("output/main_text_figures/fig4/mixed_model_test_set_performance.rds") +
  theme(axis.text = element_text(size = 9),
        legend.text = element_text(size = 9, face = "bold.italic")) +
  scale_y_continuous(breaks = scales::pretty_breaks(5))
fig4d$layers[[1]] <- NULL

tot_micob_load <- ggarrange(fig4a, fig4b, ncol = 2, labels = "AUTO",
                            common.legend = TRUE, legend = "right")
mixed_effects <- ggarrange(fig4c, fig4d, ncol = 2, labels = c("C", "D"),
                           widths = c(.6, .4),
                           legend = "bottom", common.legend = TRUE)

figure4 <- ggarrange(
  tot_micob_load,
  mixed_effects,
  nrow = 2,
  heights = c(.35, .65)
)

ggsave("output/main_text_figures/final_figures/figure4.png", figure4,
       height = 14, width = 12.5, dpi = 400)

# Fig 5

p <- readRDS("output/supp_material/supp3_mixed_model/lgo_plot.rds") +
  theme(axis.text.x = element_text(size = 9, face = "bold", color = "gray40"),
        legend.position = "bottom",
        legend.text = element_text(size = 9, face = "bold.italic", color = "gray40")) +
  guides(fill = guide_legend(nrow = 2))
ggsave("output/main_text_figures/final_figures/figure5.png", p,
       width = 6.5, height = 9, dpi = 400)
