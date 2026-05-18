############################################################
# Title: Data analysis Wolfwinkel 2 Pollen dataset
#
# Author: Wouter Driessen
# Date: 1 july 2025
#
# Description:
# This R script performs paleoecological and multivariate
# analyses on the Wolfwinkel 2 pollen dataset and compares
# it with regional Neotoma pollen records from the
# Weichselian, Early Holocene, Mid Holocene, and Late
# Holocene periods.
#
# The script:
# - Imports and preprocesses pollen datasets
# - Removes rare taxa and harmonizes taxonomic names
# - Performs Detrended Correspondence Analysis (DCA)
# - Creates DCA ordination plots
# - Performs CONISS cluster analysis and broken stick analysis
# - Integrates regional Neotoma pollen datasets
# - Harmonizes taxa between local and regional datasets
# - Compares arboreal and Poaceae pollen abundances
# - Produces boxplots and ecological visualizations
# - Examines correlations between charcoal abundance
#   and selected pollen taxa
# - Exports statistical results as CSV files
############################################################

install.packages("permute")
install.packages("lattice")
install.packages("vegan")
install.packages("analogue")
install.packages("rioja")

library(vegan)
library(analogue)
library(rioja)

setwd("~/Documents/1. Projects/BA project/BA_data_analysis")
data <- read.csv('C2_pol1.0.csv')
data <- data[1:27,]

 
# omitting rare taxa (only three appearences)
taxa <- chooseTaxa(data[-c(6, 18, 20, 21, 22), 8:27], n.occ = 5)

### Correcting Quercus 
artemisia_value <- taxa[8, 15]
taxa[8, 3] <- taxa[8, 3] + artemisia_value
taxa[8, 14] <- 0


### Wolfwinkel 2 Zonation 


# Run DCA
dca <- decorana(taxa)
sc.dca <- scores(dca, display = "sites", choices = 1:2)
dca

zone_vector <- c(rep(1, 5), rep(2, 3), rep(3, 14))
sample_labels <- c(1, 3, 5, 7, 9, 13, 14, 15, 17, 19, 21, 23, 25, 27, 29, 31, 34, 38, 39, 40, 41, 42)

dca_df <- data.frame(
  DCAaxis1 = sc.dca[, 1],
  DCAaxis2 = sc.dca[, 2],
  zone = factor(zone_vector),
  label = sample_labels
)

zone_colors <- c("1" = "indianred", "2" = "blue", "3" = "purple")

# Wolfwinkel 2 DCA Plot  
ggplot(dca_df, aes(x = DCAaxis1, y = DCAaxis2, color = zone, label = label)) +
  geom_point(size = 5) +
  geom_text_repel(size = 6, max.overlaps = Inf) +
  scale_color_manual(values = zone_colors) +
  labs(title = "DCA Wolfwinkel 2 Pollen Zones 1-3", x = "DCA1", y = "DCA2") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 22),  # increase x-axis title size
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),   # axis tick labels size
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 22),  # legend title size
    legend.text = element_text(size = 20)# increase y-axis title size
  )


# Wolfwinkel 2 CONISS dendogram
clust <- chclust(dist(taxa), method = "coniss")
clust$labels <- sample_labels

plot(clust, hang = -1, horiz = TRUE, x.rev = TRUE, cex = 1.5, xaxt = "n",
     main = "CONISS Dendrogram with Sample Labels")
# plot axis mamually
axis(side = 1, cex.axis = 1.5, font = 2)  # font=2 makes labels bold
mtext("Dissimilarity", side = 1, line = 3, cex = 2, font = 2)



# Wolfwinkel 2 broken Stick analysis 
bstick_data <- bstick(clust, plot = FALSE)

plot(bstick_data$dispersion, type = "b", col = "red", lwd = 2, pch = 19,
     xlab = "Number of Clusters", ylab = "Eigenvalue Dispersion",
     cex.axis = 2, cex.lab = 2, font.lab = 2,  # <-- here: bold axis labels
     main = "Determining Optimal Clusters Using the Broken Stick Model",
     cex.main = 2)

lines(bstick_data$bstick, type = "b", col = "black", lwd = 2, pch = 19)

legend("topright",
       legend = c("Observed eigenvalues (dispersion)", "Broken stick values"),
       col = c("red", "black"),
       lty = 1,
       lwd = 2,
       cex = 1.8)



# Regional analysys with neotoma data 
ww2_taxa <-  taxa

# Neotoma data 

# Read Neotoma datasets pulled by Neotoma_data_pull.R
neotoma_late_hol <- read.csv("neotoma_late_hol1.csv")
neotoma_mid_hol <- read.csv("neotoma_mid_hol")
neotoma_early_hol <- read.csv("neotoma_early_hol")
neotoma_weich <- read.csv("neotoma_weich")

# Select taxa
taxta_late_hol_names <- chooseTaxa(neotoma_late_hol, n.occ = 5)
taxa_late_hol <- chooseTaxa(neotoma_late_hol[,4:677], n.occ = 5)
taxa_late_hol <- taxa_late_hol[!apply(taxa_late_hol < 0, 1, any), ]
taxa_late_hol <- taxa_late_hol[rowSums(taxa_late_hol) > 0, ]

taxa_mid_hol_names <- chooseTaxa(neotoma_mid_hol, n.occ = 5)
taxa_mid_hol <- chooseTaxa(neotoma_mid_hol[,4:677], n.occ = 5)
taxa_mid_hol <- taxa_mid_hol[!apply(taxa_mid_hol < 0, 1, any), ]

taxa_early_hol_names <- chooseTaxa(neotoma_early_hol, n.occ = 5)
taxa_early_hol <- chooseTaxa(neotoma_early_hol[,4:677], n.occ = 5)
taxa_early_hol <- taxa_early_hol[!apply(taxa_early_hol < 0, 1, any), ]

taxa_weich_names <- chooseTaxa(neotoma_weich, n.occ = 5)
taxa_weich <- chooseTaxa(neotoma_weich[,4:677], n.occ = 5)
taxa_weich <- taxa_weich[!apply(taxa_weich < 0, 1, any), ]


# aligning ww2 taxon names 
rename_map <- c(
  "ulmus"               = "Ulmus",
  "betula"              = "Betula",
  "quercus"             = "Quercus",
  "pinus"               = "Pinus",
  "tilia"               = "Tilia",
  "alnus"               = "Alnus",
  "corylus"             = "Corylus",
  "salix"               = "Salix",
  "ilex"                = "Ilex",
  "pteropsida"          = "Pteropsida",
  "poaceae"             = "Poaceae",
  "cereal"              = "Poaceae (grainsize >37μm, annulus > 7)",
  "ericaceaea"          = "Ericaceae",
  "artimisia"           = "Artemisia",
  "asteraceae_tub"      = "Asteraceae tubuliflorae",
  "asteraceae_lig"      = "Asteraceae liguliflorae",
  "caryophyllaceae"     = "Caryophyllaceae",
  "cyperaceae"          = "Cyperaceae",
  "plantago_lanceolata" = "Plantago.lanceolata",
  "mentha"              = "Mentha.type"
)

for (old_name in names(rename_map)) {
  new_name <- rename_map[[old_name]]
  if (old_name %in% names(ww2_taxa)) {
    names(ww2_taxa)[names(ww2_taxa) == old_name] <- new_name
  }
}


ww2_taxa <- ww2_taxa %>%
  mutate(
    Asteraceae = `Asteraceae tubuliflorae` + `Asteraceae liguliflorae`,
    Poaceae = Poaceae + `Poaceae (grainsize >37μm, annulus > 7)`,
    Pteropsida = pteropsida_monolete + pteropsida_trilete) %>% 
  select(Asteraceae, Poaceae, Pteropsida, Ulmus, Betula, 
         Quercus, Pinus, Tilia, Alnus, Salix, Corylus, Ilex, Ericaceae, Artemisia, Caryophyllaceae)
combined_taxa_neotoma <- bind_rows(taxa_late_hol, taxa_mid_hol, taxa_early_hol, taxa_weich)
str(combined_taxa_neotoma)


# Harmonization

# Input data frames
ww2_taxa <- ww2_taxa
full_taxa <- combined_taxa_neotoma

# Step 1: Shared taxa
shared_taxa <- intersect(names(ww2_taxa), names(full_taxa))

# Step 2: Start harmonized frames with shared taxa
harmonized_ww2 <- ww2_taxa[, shared_taxa, drop = FALSE]
harmonized_full <- full_taxa[, shared_taxa, drop = FALSE]


# Step 3: Define additional taxon mappings
extra_mappings <- list(
  Asteraceae = c(
    "Anthemis.Aster.Cirsium.Serratula",
    "Anthemis.type",
    "Aster.type",
    "Achillea.type",
    "Asteroideae",
    "Asteroideae.undiff.",
    "Centaurea.cyanus",
    "Centaurea.jacea.type",
    "Cirsium",
    "Cirsium.Carduus",
    "Centaurea.nigra.type"
  ),
  
  Poaceae = c(
    "Poaceae",
    "Poaceae.cf..Hordeum.group",
    "Poaceae..Cerealia.",
    "Poaceae..Cerealia.type.",
    "Poaceae..Cerealia..excluding.Secale",
    "Poaceae..Cerealia..undiff.",
    "Poaceae.indet.",
    "Hordeum.type",
    "Avena.Triticum.type",
    "Secale.cereale",
    "Zea.mays",
    "Poaceae...35.µm.",
    "Poaceae.undiff."
  ),
  
  Pteropsida = c(
    "Pteropsida",
    "Dryopteris.type",
    "Dryopteris.Thelypteris",
    "Polypodium.vulgare.type",
    "Polypodium.vulgare",
    "Pteridium.aquilinum",
    "Thelypteris.palustris.type",
    "Thelypteris.palustris",
    "Ophioglossum.vulgatum",
    "Selaginella.selaginoides",
    "Selaginella",
    "Huperzia.selago"
  ),
  
  Ulmus = c(
    "Ulmus",
    "Ulmus.glabra.type",
    "Ulmus.glabra"
  ),
  
  Betula = c(
    "Betula",
    "Betula.pubescens.sensu.lato",
    "Betula.pubescens.type",
    "Betula.pendula.B..pubescens",
    "Betula.undiff.",
    "Betula.nana.type",
    "Betula.nana"
  ),
  
  Quercus = c(
    "Quercus",
    "Quercus.robur.type",
    "Quercus.robur.group"
  ),
  
  Pinus = c(
    "Pinus",
    "Pinus.undiff.",
    "Pinus.sylvestris",
    "Pinus.sylvestris.type",
    "Pinus.haploxylon"
  ),
  
  Tilia = c(
    "Tilia",
    "Tilia.cordata",
    "Tilia.undiff."
  ),
  
  Alnus = c(
    "Alnus",
    "Alnus.glutinosa.type",
    "Alnus.glutinosa",
    "Frangula.alnus"
  ),
  
  Salix = c(
    "Salix",
    "Salix.type",
    "Salix.caprea",
    "Salix.glauca.type",
    "Salix.pentandra.type",
    "Salix.herbacea.type",
    "Salix.herbacea",
    "Salix.polaris.type",
    "Salix.polaris",
    "Lythrum.salicaria",
    "Lythrum.portula"
  ),
  
  Corylus = c(
    "Corylus",
    "Corylus.avellana",
    "Betula.Alnus.Corylus.Myrica..type.114..HdV."
  ),
  
  Ilex = c(
    "Ilex",
    "Ilex.aquifolium"
  ),
  
  Ericaceae = c(
    "Ericaceae",
    "Ericales",
    "Ericales.undiff.",
    "Ericaceae.undiff.",
    "Vaccinium.myrtillus.type",
    "Calluna",
    "Empetrum",
    "Hippophaë"
  ),
  
  Artemisia = c(
    "Artemisia",
    "Artemisia.type",
    "Artemisia.vulgaris.type"
  ),
  
  Caryophyllaceae = c(
    "Caryophyllaceae",
    "Caryophyllaceae.undiff."
  )
)

# Step 4: Sum mapped taxa in full_taxa and inject into harmonized_ww2/full
for (taxon in names(extra_mappings)) {
  subtypes <- extra_mappings[[taxon]]
  existing_subtypes <- subtypes[subtypes %in% names(full_taxa)]
  
  if (length(existing_subtypes) > 0) {
    # Always (re)add to harmonized_full
    harmonized_full[[taxon]] <- rowSums(full_taxa[, existing_subtypes, drop = FALSE], na.rm = TRUE)
    
    # Copy from ww2_taxa if it exists there
    if (taxon %in% names(ww2_taxa)) {
      harmonized_ww2[[taxon]] <- ww2_taxa[[taxon]]
    } else {
      # If not present, fill with 0
      harmonized_ww2[[taxon]] <- 0
    }
  }
}

# Step 5: Final re-alignment
common_taxa_final <- intersect(names(harmonized_ww2), names(harmonized_full))
harmonized_ww2 <- harmonized_ww2[, common_taxa_final]
harmonized_full <- harmonized_full[, common_taxa_final]

# Joining data frames
combined_taxa <- bind_rows(harmonized_ww2,harmonized_full)

combined_taxa[is.na(combined_taxa)] <- 0

#omit corylus 
combined_taxa <- combined_taxa[,-10]




#DCA 
dca <- decorana(combined_taxa)
sc.dca <- scores(dca, display = "sites", choices = 1:2)

zone_vector <- c(rep(1, 5), rep(2, 3), rep(3, 14), rep("Late", 526), rep('Mid', 93), rep('Early', 80), rep('Weich', 288))

ww2_taxa_dca <- data.frame(
  DCAaxis1 = sc.dca[, 1],
  DCAaxis2 = sc.dca[, 2],
  zone = factor(zone_vector) # Assigning the zone_vector to zone column
)



# Your existing zone_vector (small zones + big groups)
zone_vector <- c(rep(1, 5), rep(2, 3), rep(3, 14), rep("Late", 526), rep('Mid', 93), rep('Early', 80), rep('Weich', 288))


zone_colors <- c(
  "1" = "indianred",
  "2" = "blue",
  "3" = "purple",
  "Late" = "#8C7B75",  # warm gray-brown
  "Mid" = "#0072B2",
  "Early" = "#D62728",  # bright red
  "Weich" = "#009E73"
)


zone_shapes <- c("1" = 17, "2" = 17, "3" = 17, 
                 "Late" = 1,  "Mid" = 1, "Early" = 1, "Weich" = 1) 
zone_sizes <- c("1" = 1.7, "2" = 1.7, "3" = 1.7, 
                "Late" = 0.8,  "Mid" = 0.8, "Early" = 0.8, "Weich" = 0.8)
transparent_zones <- c("Late", "Mid", "Early", "Weich")

zone_colors[transparent_zones] <- adjustcolor(zone_colors[transparent_zones], alpha.f = 0.8)
ww2_taxa_dca$zone <- factor(ww2_taxa_dca$zone, levels = c("1", "2", "3", "Late", "Mid", "Early", "Weich"))

# Filter the dataset first
filtered <- ww2_taxa_dca$zone %in% c("Late","Mid", "Early", "Weich")


# Define the specific labels for your 22 samples IN THE CORRECT ORDER
sample_labels <- c(1, 3, 5, 7, 9, 13, 14, 15, 17, 19, 21, 23, 25, 27, 29, 31, 34, 38, 39, 40, 41)


# Identify rows belonging to zones 1, 2, 3
zones_123 <- ww2_taxa_dca$zone %in% c("1", "2", "3")

# Subset your coordinates for zones 1-3
coords_z123 <- ww2_taxa_dca[zones_123, ]

# Your sample labels (should be length 22 if you have 22 rows)
labels_z123 <- sample_labels

# Check lengths match
if (length(labels_z123) != nrow(coords_z123)) {
  stop("Sample label count does not match number of samples in zones 1-3.")
}

# Map zone colors to text labels
label_colors <- zone_colors[as.character(coords_z123$zone)]

# Plot as usual
plot(
  ww2_taxa_dca$DCAaxis1[filtered],
  ww2_taxa_dca$DCAaxis2[filtered],
  col = zone_colors[ww2_taxa_dca$zone[filtered]],
  pch = zone_shapes[ww2_taxa_dca$zone[filtered]],
  cex = zone_sizes[ww2_taxa_dca$zone[filtered]],
  xlab = "DCA1",
  ylab = "DCA2",
  main = "DCA: Wolfwinkel 2 sequence vs. Regional Neotoma Pollen Records",
  xlim = range(ww2_taxa_dca$DCAaxis1, na.rm = TRUE),
  ylim = range(ww2_taxa_dca$DCAaxis2, na.rm = TRUE),
  cex.lab = 1.5,    # bigger axis labels
  cex.axis = 1.3,   # bigger axis tick numbers
  cex.main = 1.6    # bigger plot title
)

custom_legend_labels <- c("Wolfwinkel 2: Zone 3", "Wolfwinkel 2: Zone 2", "Wolfwinkel 2: Zone 1", "Late Holocene (Present to 4.2 Ka BP)",
                          "Middle Holocene (4.2 - 8.2 ka BP)", "Early Holocene (8.2 - 11.7 ka BP)", "Weichselian (11.7 - 115 ka BP)")

# Plot points for zones 1-3
points(
  coords_z123$DCAaxis1,
  coords_z123$DCAaxis2,
  col = zone_colors[coords_z123$zone],
  pch = zone_shapes[coords_z123$zone],
  cex = zone_sizes[coords_z123$zone]
)



legend(
  "topright",
  legend = custom_legend_labels,
  pch = zone_shapes[levels(ww2_taxa_dca$zone)],
  col = zone_colors[levels(ww2_taxa_dca$zone)],
  title = title("Legend", font = 2),
  text.font = 2,
  cex = 1
)


# Alternate label positions
positions <- rep(c(3, 3, 2, 4), length.out = nrow(coords_z123))

text(
  coords_z123$DCAaxis1,
  coords_z123$DCAaxis2,
  labels = labels_z123,
  pos = positions,
  col = label_colors,
  cex = 1.2
)


# Boxplot Poacheae VS Arboreal 


library(dplyr)
library(tidyr)
library(ggplot2)

custom_legend_labels <- c(
  "Wolfwinkel 2: Zone 3",
  "Weichselian (11.7 - 115 ka BP)",
  "Early-Mid Holocene (11.7 - 4.2 ka BP)"
)

plot_data <- combined_taxa %>%
  mutate(
    Arboreal = Ulmus + Betula + Quercus + Pinus + Tilia,
    # Shrubs omitted as requested
    period = c(
      rep("Zone 3", 22),
      rep("Late holocene", 526),
      rep("Mid holocene", 93),
      rep("Early  holocene", 80),
      rep("Weichselian", 288)
    )
  ) %>%
  filter(period != "Late holocene") %>%  # omit Late Holocene
  mutate(
    period = case_when(
      period %in% c("Early  holocene", "Mid holocene") ~ "Early-Mid Holocene",
      TRUE ~ period
    )
  ) %>%
  filter(
    (period != "Zone 3") |
      (period == "Zone 3" & row_number() > 8)  # Keep last 14 of Zone 3
  ) %>%
  select(period, Arboreal, Poaceae) %>%
  pivot_longer(
    cols = c(Arboreal, Poaceae),
    names_to = "Taxa",
    values_to = "Value"
  ) %>%
  mutate(
    period = factor(
      period,
      levels = c("Zone 3", "Weichselian", "Early-Mid Holocene"),
      labels = custom_legend_labels
    )
  )

# Highlight points sample 41 and 42
highlight_points <- plot_data %>%
  filter(period == custom_legend_labels[1]) %>%
  group_by(Taxa) %>%
  slice_tail(n = 2) %>%
  ungroup() %>%
  mutate(Highlight = factor("Sample 41 and 42", levels = c("Sample 41 and 42")))

# Calculate sample sizes per period & Taxa
sample_sizes <- plot_data %>%
  group_by(period, Taxa) %>%
  summarise(n = n(), .groups = "drop")

# Define custom colors for legend fill
custom_colors <- c(
  "Wolfwinkel 2: Zone 3" = "#1b9e77",
  "Weichselian (11.7 - 115 ka BP)" = "#d95f02",
  "Early-Mid Holocene (11.7 - 4.2 ka BP)" = "#7570b3"
)

# Plot
ggplot(plot_data, aes(x = period, y = Value, fill = period)) +
  geom_boxplot() +
  geom_point(
    data = highlight_points,
    aes(x = period, y = Value, shape = Highlight, color = Highlight),
    size = 3, stroke = 1.5, fill = "yellow"
  ) +
  geom_text(
    data = sample_sizes,
    aes(x = period, y = -Inf, label = paste0("n=", n)),
    vjust = -0.5,
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = custom_colors,
    name = "Periods"
  ) +
  scale_color_manual(
    name = "Highlighted Points",
    values = c("Sample 41 and 42" = "black")
  ) +
  scale_shape_manual(
    name = "Highlighted Points",
    values = c("Sample 41 and 42" = 21)
  ) +
  facet_wrap(~ Taxa, scales = "free_y") +
  labs(
    title = "Poaceae and Arboreal Pollen Comparison: Wolfwinkel 2 vs Neotoma by Period",
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold", family = "sans")
  )



### Corralations 
# Cahrcoal Volume_mm3 vector
Volume_mm3 <- c(
  3.409, 4.628, 2.802, 1.987, 4.201, 0.924, 0.349, 8.027, 1.766, 0,
  0.013, 6.519, 0.368, 0.601, 1.79, 0.049, 0.073, 2.076, 8.194, 1.825,
  2.471, 2.081, 1.953, 1.038, 3.884, 6.026, 0.172, 0.233, 0.454, 0.038,
  0.276, 0, 1.543, 4.116, 0, 0, 0, 0, 0, 0, 0, 0
)

# Interpolate Volume_mm3 to match ww2_taxa rows
interp_vol <- approx(
  x = seq_along(Volume_mm3),
  y = Volume_mm3,
  xout = seq(from = 1, to = length(Volume_mm3), length.out = nrow(ww2_taxa))
)$y

# Combine with taxa data
ww2_taxa_ext <- cbind(ww2_taxa, Interpolated_Volume_mm3 = interp_vol)

# Function to fit models and create summary table
fit_and_summarize <- function(data){
  model_quercus <- lm(Interpolated_Volume_mm3 ~ Quercus, data = data)
  model_ericaceae <- lm(Interpolated_Volume_mm3 ~ Ericaceae, data = data)
  model_poaceae <- lm(Interpolated_Volume_mm3 ~ Poaceae, data = data)
  
  data.frame(
    Predictor = c("Quercus", "Ericaceae", "Poaceae"),
    Correlation = c(
      round(cor(data$Quercus, data$Interpolated_Volume_mm3, use = "complete.obs"), 3),
      round(cor(data$Ericaceae, data$Interpolated_Volume_mm3, use = "complete.obs"), 3),
      round(cor(data$Poaceae, data$Interpolated_Volume_mm3, use = "complete.obs"), 3)
    ),
    Estimate = c(
      round(coef(model_quercus)["Quercus"], 3),
      round(coef(model_ericaceae)["Ericaceae"], 3),
      round(coef(model_poaceae)["Poaceae"], 3)
    ),
    P_value = c(
      summary(model_quercus)$coefficients["Quercus", "Pr(>|t|)"],
      summary(model_ericaceae)$coefficients["Ericaceae", "Pr(>|t|)"],
      summary(model_poaceae)$coefficients["Poaceae", "Pr(>|t|)"]
    )
  )
}

# Results for subset rows 11 to 17
subset_df <- ww2_taxa_ext[11:17, ]
results_subset <- fit_and_summarize(subset_df)
cat("Results for subset (rows 11 to 17):\n")
print(results_subset)

# Results for all rows
results_all <- fit_and_summarize(ww2_taxa_ext)
cat("\nResults for all rows:\n")
print(results_all)

# Save subset results as CSV
write.csv(results_subset, file = "results_subset.csv", row.names = FALSE)

# Results for all rows
results_all <- fit_and_summarize(ww2_taxa_ext)

# Save all results as CSV
write.csv(results_all, file = "results_all.csv", row.names = FALSE)



