# Environmental characterization of all vent sites
## Three sites in separate ocean regions.

metadata <- read.delim("data-input/samplelist-metadata.txt", na.strings = "")

# Filter to environmental data only.
# head(metadata)
# names(metadata)

env_deepsea <- metadata %>% 
  mutate_all(as.character) %>%
  filter(Sample_or_Control == "Sample") %>% 
  filter(!(SAMPLETYPE == "Incubation")) %>% 
  filter(!(SAMPLETYPE == "Microcolonizer")) %>% 
  select(VENT, COORDINATES, SITE, SAMPLEID, DEPTH, SAMPLETYPE, YEAR, TEMP = starts_with("TEMP"), pH, PercSeawater = starts_with("Perc"), Mg = starts_with("Mg"), H2 = starts_with("H2."), H2S = starts_with("H2S"), CH4 = starts_with("CH4"), ProkConc) %>%  
  pivot_longer(cols = TEMP:ProkConc, names_to = "VARIABLE", values_to = "VALUE", values_drop_na = FALSE) %>% 
  distinct() %>% 
  group_by(VENT, COORDINATES, SITE, DEPTH, SAMPLETYPE, YEAR, VARIABLE) 

# head(env_deepsea)
# unique(env_deepsea$VARIABLE)
# str(env_deepsea)


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
geomchem_table <- metadata %>% 
  mutate_all(as.character) %>%
  filter(Sample_or_Control == "Sample") %>% 
  filter(!(SAMPLETYPE == "Incubation")) %>% 
  filter(!(SAMPLETYPE == "Microcolonizer")) %>% 
  select(SAMPLE, VENT, COORDINATES, SITE, SAMPLEID, DEPTH, SAMPLETYPE, YEAR, TEMP = starts_with("TEMP"), pH, PercSeawater = starts_with("Perc"), Mg = starts_with("Mg"), H2 = starts_with("H2."), H2S = starts_with("H2S"), CH4 = starts_with("CH4"), ProkConc) %>% 
  group_by(SITE, SAMPLETYPE, YEAR, VENT, COORDINATES, DEPTH, TEMP, pH, PercSeawater, Mg, H2, H2S, CH4, ProkConc) %>% 
  summarize(SampleIDs = paste(unique(SAMPLEID), collapse = ", "))
# geomchem_table

write_delim(geomchem_table, file = "tableS1-geochem-params.txt", delim = "\t")
