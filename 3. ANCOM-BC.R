setwd("C:/Users/sbald/Dropbox/Simone/DTU/10. semester/Speciale/Bioinformatic analysis/R/R scripts")

load('R workspaces/ps.RData')

# Libraries
library(ANCOMBC)
library(microbiome)
library(DT)

## ANCOM-BC analysis ##

# Aggregate to a taxa level (Phylum, Class, Order, Family, Genus)
data = aggregate_taxa(ps.pcs, "Genus") # with P. piscicida B39bio
data.noB39bio = aggregate_taxa(ps.noB39bio, "Genus") # without P. piscicida B39bio

# Run ancom-bc
ancom_output = ancombc(phyloseq = data, formula = "Setup_2", 
                       p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                       group = "Setup_2", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                       max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

ancom_output.noB39bio = ancombc(phyloseq = data.noB39bio, formula = "Setup_2", 
                              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                              group = "Setup_2", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)


res = ancom_output$res
res.noB39bio = ancom_output.noB39bio$res

res_global = ancom_output$res_global
res_global.noB39bio = ancom_output.noB39bio$res_global

# P-values (p_val)
# with B39bio
col_name = c("Piscicida-Control", "Sewater-Control")
tab_p = res$p_val
colnames(tab_p) = col_name
write.table(tab_p, file = "Files/p_values_genus.csv", sep = ";",
            row.names = TRUE, col.names = NA)
# without B39bio
tab_p.noB39bio = res.noB39bio$p_val
colnames(tab_p.noB39bio) = col_name
write.table(tab_p.noB39bio, file = "Files/p_values_genus_noB39bio.csv", sep = ";",
            row.names = TRUE, col.names = NA)

# Differentially abundant taxa
# with B39bio
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
write.table(tab_diff, file = "Files/diff_abn_genus.csv", sep = ";",
            row.names = TRUE, col.names = NA)
# without B39bio
tab_diff.noB39bio = res.noB39bio$diff_abn
colnames(tab_diff.noB39bio) = col_name
write.table(tab_diff.noB39bio, file = "Files/diff_abn_genus.noB39bio.csv", sep = ";",
            row.names = TRUE, col.names = NA)


## Log fold change plot (P. piscicida B39bio vs. Control) - with B39bio ##

df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
# naming columns of df_fig2 SD (standard deviation) in the end, so they dont have the same column names as df_fig1
colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
# Merge the two dataframes
df_fig = df_fig1 %>% left_join(df_fig2, by = "taxon_id") %>%
  transmute(taxon_id, Setup_2Piscicida, Setup_2PiscicidaSD) %>%
  filter(Setup_2Piscicida != 0) %>%
  filter(Setup_2Piscicida > 0.3 | Setup_2Piscicida < -0.3)
df_fig$taxon_id = factor(df_fig$taxon_id, levels = df_fig$taxon_id)

# Plot
p = ggplot(data = df_fig, 
                  aes(x = reorder(taxon_id, Setup_2Piscicida), y = Setup_2Piscicida)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Setup_2Piscicida - Setup_2PiscicidaSD, ymax = Setup_2Piscicida + Setup_2PiscicidaSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "B39bio vs. Control - Family") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12))
p


## Log fold change plot (P. piscicida B39bio vs. Control) - without B39bio ##

df_fig1.noB39bio = data.frame(res.noB39bio$beta * res.noB39bio$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_fig2.noB39bio = data.frame(res.noB39bio$se * res.noB39bio$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
# naming columns of df_fig2 SD (standard deviation) in the end, so they dont have the same column names as df_fig1
colnames(df_fig2.noB39bio)[-1] = paste0(colnames(df_fig2.noB39bio)[-1], "SD")
# Merge the two dataframes
df_fig.noB39bio = df_fig1.noB39bio %>% left_join(df_fig2.noB39bio, by = "taxon_id") %>%
  transmute(taxon_id, Setup_2Piscicida, Setup_2PiscicidaSD) %>%
  filter(Setup_2Piscicida != 0) %>%
  filter(Setup_2Piscicida > 0.3 | Setup_2Piscicida < -0.3)
df_fig.noB39bio$taxon_id = factor(df_fig.noB39bio$taxon_id, levels = df_fig.noB39bio$taxon_id)

# Plot
p.noB39bio = ggplot(data = df_fig.noB39bio, 
           aes(x = reorder(taxon_id, Setup_2Piscicida), y = Setup_2Piscicida)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Setup_2Piscicida - Setup_2PiscicidaSD, ymax = Setup_2Piscicida + Setup_2PiscicidaSD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "B39bio vs. Control - Family, noB39bio ASVs") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12))
p.noB39bio



## Bias adjusted abundances (with B39bio) - only run at phylum level ##

samp_frac = ancom_output$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(abundances(data) + 1) 
# Adjust the log observed abundances
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_obs_abn_adj[, 6:12], 2) %>% 
  datatable(caption = "Bias-adjusted log observed abundances")

log_obs_abn_adj_df <- as.data.frame(log_obs_abn_adj)
write.csv(log_obs_abn_adj_df, file = "Files/abundance_phylum.csv")
# File 'abundance_sd_phylum.csv' was made in excel
abundance_sd_phylum <- read.csv(file = "Files/abundance_sd_phylum.csv", sep = ";")
abundance_sd_phylum = data.frame(abundance_sd_phylum)

# Plotting bias-adjusted log observed abundances
ggplot(abundance_sd_phylum, aes(x = ï..taxon_id, y = abundance, fill = setup)) +
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = abundance - sd, ymax = abundance + sd),
                position = "dodge", color = "black") + 
  labs(x = NULL, y = "Bias-adjusted log observed abundances", 
       title = "") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12))
