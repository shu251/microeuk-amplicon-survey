# Environmental characterization of all vent sites
## Three sites in separate ocean regions.

library(tidyverse); library(cowplot)
library(phyloseq); library(decontam)
library(compositions); library(patchwork); library(ggpubr)
library(ggupset); library(gt); library(treemapify)
library(plotly); library(viridis); library(vegan)


metadata <- read.csv("data-input/samplelist-metadata.csv", na.strings = "")

# Filter to environmental data only.
# head(metadata)
# names(metadata)

env_deepsea <- metadata %>% 
  mutate_all(as.character) %>%
  filter(Sample_or_Control == "Sample") %>% 
  filter(!(SAMPLETYPE == "Incubation")) %>% 
  filter(!(SAMPLETYPE == "Microcolonizer")) %>% 
  select(SAMPLE, VENT, COORDINATES, SITE, SAMPLEID, DEPTH, SAMPLETYPE, YEAR, TEMP = starts_with("TEMP"), pH, PercSeawater = starts_with("Perc"), Mg = starts_with("Mg"), H2 = starts_with("H2."), H2S = starts_with("H2S"), CH4 = starts_with("CH4"), ProkConc) %>%  
  pivot_longer(cols = TEMP:ProkConc, names_to = "VARIABLE", values_to = "VALUE", values_drop_na = FALSE) %>% 
  distinct() %>% 
  group_by(SAMPLE, VENT, COORDINATES, SITE, DEPTH, SAMPLETYPE, YEAR, VARIABLE) 

# head(env_deepsea)
# unique(env_deepsea$VARIABLE)
# str(env_deepsea)

## Fix coordinates
unique(env_deepsea$SITE)
gr_axial_coord <- env_deepsea %>% 
  ungroup() %>% 
  filter((SITE == "Axial" | SITE == "GordaRidge")) %>% 
  select(VENT, COORDINATES, SITE, SAMPLETYPE) %>% 
  distinct() %>% 
  separate(COORDINATES, into = c("lat", "N", "long", "W"), sep = " ") %>%
  mutate(
    LONG_EW = as.numeric(formatC(as.numeric(long), digits = 4, format = "f")),
    LAT = as.numeric(formatC(as.numeric(lat), digits = 4, format = "f")),
  ) %>%
  mutate(LONG = case_when(
    W == "W" ~ (LONG_EW*-1),
    W == "E" ~ LONG_EW,
    W == "" ~ (LONG_EW*-1)
  )) %>%
  select(-lat, -N, -long, -W, -LONG_EW) %>%
  relocate(LONG)

mcr_coord <- env_deepsea %>% 
  ungroup() %>% 
  filter((SITE == "Piccard" | SITE == "VonDamm")) %>% 
  select(VENT, COORDINATES, SITE, SAMPLETYPE) %>% 
  distinct() %>% 
  separate(COORDINATES, into = c("lat", "long"), sep = ", ", convert = TRUE) %>%
  mutate(
    LONG = as.numeric(formatC(as.numeric(long), digits = 4, format = "f")),
    LAT = as.numeric(formatC(as.numeric(lat), digits = 4, format = "f")),
  ) %>% select(-lat, -long)

all_coords <- bind_rows(mcr_coord, gr_axial_coord)
# View(all_coords)

# Units for the variables are as follows:
#   - Temp = Celsius
# - Percent Seawater = %
# - Mg = mmol/kg (or mM)
# - H2 = µmol/L (or µM)
# - H2S = mmol/L (or mM)
# - CH4 = µmol/L (or µM)
# - ProkConc = cells/ml

## Visualize geochemistry

rm <- c("-", "", "nd", "bd", NA)

geochem_violin <- env_deepsea %>% 
  filter(VARIABLE != "ProkConc") %>% 
  filter(!(VALUE %in% rm)) %>% 
  mutate(VALUE = as.numeric(as.character(VALUE))) %>% 
  ggplot(aes(x = SAMPLETYPE, y = VALUE, fill = SITE, shape = SAMPLETYPE)) +
  geom_boxplot(alpha = 0.3, aes(group = SAMPLETYPE), fill = "grey", width = 0.3) +
  geom_jitter(size = 2, width = 0.2) +
  facet_wrap(VARIABLE ~ ., scales = "free") +
  scale_shape_manual(values = c(21, 23, 24)) +
  scale_fill_manual(values = c("#fdbb84", "#31a354", "#ef3b2c", "#02818a")) +
  guides(fill = guide_legend(override.aes = list(shape = 21) ),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.text = element_text(color="black", size=12),
        legend.title = element_blank(),
        axis.title = element_text(color="black", size=14),
        legend.text = element_text(color = "black", size = 14),
        plot.margin = margin(2, 1, 2, 1, "cm"),
        strip.background = element_blank()) +
  labs(x = "", y = "")
geochem_violin


geochem_2 <- env_deepsea %>% 
  filter(VARIABLE == "ProkConc") %>% 
  filter(!(VALUE %in% rm)) %>% 
  mutate(VALUE = as.numeric(as.character(VALUE))) %>% 
  ggplot(aes(x = SAMPLETYPE, y = VALUE, fill = SITE, shape = SAMPLETYPE)) +
  geom_boxplot(alpha = 0.3, aes(group = SAMPLETYPE), fill = "grey", width = 0.3) +
  geom_jitter(size = 2, width = 0.2) +
  facet_wrap(VARIABLE ~ ., scales = "free") +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 23, 24)) +
  scale_fill_manual(values = c("#fdbb84", "#31a354", "#ef3b2c", "#02818a")) +
  guides(fill = guide_legend(override.aes = list(shape = 21) ),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  theme(axis.text = element_text(color="black", size=12),
        legend.title = element_blank(),
        axis.title = element_text(color="black", size=14),
        legend.text = element_text(color = "black", size = 14),
        plot.margin = margin(2, 1, 2, 1, "cm"),
        strip.background = element_blank()) +
  labs(x = "", y = "")
geochem_2


### Set colors and symbols for dataset
sampletype_order <- c("Background", "Plume", "Vent")
sampletype_symbol<- c(21, 23, 24)

site_order <- c("Axial", "GordaRidge", "Piccard", "VonDamm")
site_color <- c("#fdbb84", "#31a354", "#ef3b2c", "#02818a")


## **Table S1** Geochemistry values
### Generate output table for all parameters measured

# colnames(metadata)
# geomchem_table <- metadata %>% 
#   mutate_all(as.character) %>%
#   filter(Sample_or_Control == "Sample") %>% 
#   filter(!(SAMPLETYPE == "Incubation")) %>% 
#   filter(!(SAMPLETYPE == "Microcolonizer")) %>% 
#   select(-COORDINATES) %>% 
#   left_join(all_coords) %>% 
#   select(SAMPLE, VENT, SITE, SAMPLEID, LAT, LONG, DEPTH, SAMPLETYPE, YEAR, TEMP = starts_with("TEMP"), pH, PercSeawater = starts_with("Perc"), Mg = starts_with("Mg"), H2 = starts_with("H2."), H2S = starts_with("H2S"), CH4 = starts_with("CH4"), ProkConc) %>% 
#   group_by(SITE, SAMPLETYPE, YEAR, VENT, LAT, LONG, DEPTH, TEMP, pH, PercSeawater, Mg, H2, H2S, CH4, ProkConc) %>% 
#   summarize(SampleIDs = paste(unique(SAMPLEID), collapse = ", ")) %>% 
#   ungroup() 
  
# ?distinct()
# geomchem_table
# head(metadata)
geomchem_table <- metadata %>% 
  mutate_all(as.character) %>%
  filter(Sample_or_Control == "Sample") %>% 
  filter(!(SAMPLETYPE == "Incubation")) %>% 
  filter(!(SAMPLETYPE == "Microcolonizer")) %>% 
  select(SAMPLE, VENT, COORDINATES, SITE, SAMPLEID, DEPTH, SAMPLETYPE, YEAR, TEMP = starts_with("TEMP"), pH, PercSeawater = starts_with("Perc"), Mg = starts_with("Mg"), H2 = starts_with("H2."), H2S = starts_with("H2S"), CH4 = starts_with("CH4"), ProkConc) %>%  
  pivot_longer(cols = TEMP:ProkConc, names_to = "VARIABLE", values_to = "VALUE") %>% 
  mutate(VALUE = as.numeric(as.character(VALUE))) %>% 
  # group_by(SAMPLE, VENT, SITE, SAMPLETYPE, YEAR, VARIABLE) %>% 
  # summarise(VALUE_MEAN = mean(VALUE),
  #           SampleIDs = paste(unique(SAMPLEID), collapse = ", "),
  #           DEPTH_mean = mean(as.numeric(DEPTH)),
  #           ) %>% 
  pivot_wider(names_from = VARIABLE, values_from = VALUE) %>%
    left_join(all_coords %>% distinct(VENT, SITE, SAMPLETYPE, .keep_all=TRUE))
# str(geomchem_table)
head(geomchem_table)
# 
# write_delim(geomchem_table, file = "tableS1-geochem-params.txt", delim = "\t")

# View(geomchem_table)
###
# Create gt table
###

# my_color_pal <- function(x) {
#   scales::col_numeric(
#     palette = paletteer::paletteer_d(
#       palette = "ggsci::red_material"
#     ) %>% as.character(),
#     domain = NULL
#   )(x)
# }


# - Percent Seawater = %
# - Mg = mmol/L (or mM)
# - H2 = µmol/L (or µM)
# - H2S = mmol/L (or mM)
# - CH4 = µmol/L (or µM)

table_all_geochem <- geomchem_table %>% 
  ungroup() %>% 
  # select(-SampleIDs) %>%
  mutate('Temperature (°C)' = as.numeric(as.character(TEMP)),
         'Microbial concentration (cells/mL)' = as.numeric(as.character(ProkConc)),
         'Percent Seawater' = (as.numeric(as.character(PercSeawater)))/100,
         pH = as.numeric(as.character(pH)),
         'CH4 µmol/L (µM)' = as.numeric(as.character(CH4)),
         'Mg mmol/L (mM)' = as.numeric(as.character(Mg)),
         'H2 µmol/L (µM)' = as.numeric(as.character(H2)),
         'H2S mmol/L (mM)' = as.numeric(as.character(H2S)),
         'Depth (m)' = as.numeric(as.character(DEPTH))
         ) %>%
  # replace_na("na") %>% 
  select(SITE, SAMPLETYPE, VENT, Year = YEAR, -TEMP, -ProkConc, -PercSeawater, -Mg, -H2S, -H2, -CH4, 'Temperature (°C)', 'Microbial concentration (cells/mL)', 
         'Percent Seawater', 'CH4 µmol/L (µM)', pH,
         'Mg mmol/L (mM)', 'H2 µmol/L (µM)', 'H2S mmol/L (mM)', 'Depth (m)', Latitude=LAT, Longitude=LONG) %>% 
  # unite(LOCATION, SITE, SAMPLETYPE, sep = " ", remove = FALSE) %>% 
  group_by(SITE, SAMPLETYPE) %>% 
  arrange(Latitude) %>%
  filter(VENT != "Quakeplume") %>%
  # replace_na("na") %>% 
  distinct() %>%  
  gt(
    groupname_col = c("SITE", "SAMPLETYPE"),
    rowname_col = c("VENT")
  ) %>% 
  row_group_order(
  groups = c("Axial - Vent", "Axial - Plume", "Axial - Background",
             "GordaRidge - Vent", "GordaRidge - Plume", "GordaRidge - Background",
             "VonDamm - Vent", "VonDamm - Plume", "VonDamm - Background",
             "Piccard - Vent", "Piccard - Plume", "Piccard - Background")
  ) %>%
  # data_color(columns = c(Temperature, 'Microbial concentration', 
  #                        'Percent Seawater', 'CH4 µmol/L', pH,
  #                        'Mg mmol/L', 'H2 µmol/L', 'H2S mmol/L'),
  #            colors = my_color_pal) %>% 
  fmt_number(columns = c('Depth (m)'), decimals = 0) %>%
  fmt_number(columns = c('Temperature (°C)'), decimals = 1) %>% 
  fmt_scientific(columns = 'Microbial concentration (cells/mL)') %>% 
  fmt_percent(columns = 'Percent Seawater', decimals = 0) %>% 
  fmt_number(columns = c(pH), n_sigfig = 2) %>%
  fmt_number(columns = c('Mg mmol/L (mM)'), n_sigfig = 3) %>% 
  fmt_number(columns = c('CH4 µmol/L (µM)'), n_sigfig = 2) %>% 
  fmt_number(columns = c('H2 µmol/L (µM)', 'H2S mmol/L (mM)'), n_sigfig = 2) %>%
  tab_header(
    title = md("**Vent fluid parameters**")
  ) %>% 
  tab_spanner(
    label = "Location",
    columns = c('Depth (m)', Latitude, Longitude)
  ) %>% 
  tab_style(
    style = list(
      cell_fill("black"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = cells_stub()) %>% 
  tab_source_note(md("*NA represents no data, bd indicates below detection*")) %>% 
  tab_options(table.width = px(500)) %>% 
  opt_table_lines()
table_all_geochem
## Include this html file in github output files too
gtsave(table_all_geochem, filename = "geochemistry_all_gt.html", path = "output-tables/")

# install.packages("webshot")
# webshot::install_phantomjs()
# gtsave(table_all_geochem, filename = "geochemistry_all_gt.pdf", path = "output-tables/")
# gtsave(table_all_geochem, filename = "geochemistry_all_gt.png", path = "output-tables/")
