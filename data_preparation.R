#' Supplementary code to the article: 
#' Klinkovská et al. Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes
#' 
#' Author: Klára Klinkovská, 2023-06-21
#' R version 4.3.0

library(data.table)
library(sf)
library(plyr)
library(tidyverse)
library(readxl)
library(ggspatial)
library(patchwork)


# data import, join grid cells, exclude unreliable records ----------------

# unreliable data for exclusion, keep only veroh == 1 - verified records
ndop.del <- read_xlsx("munzbergova.xlsx") %>% 
  filter(VEROH != 1)

# grid cells after expert evaluation with the occurrence of the given species - green
# and not revised yet - grey
# this excludes records that were marked as implausible
greygreen <- read_csv("gray_green_qudrants_aggregated.csv") %>% 
  rename("code" = "quadrant_name")

## data import
pladias <- read_csv("records.csv") %>% 
  mutate(coords2 = str_sub(coords, 11, 50), # extract coordinates
         long = str_split_i(coords2, " ", 1) %>% 
           str_remove("POINT\\("), 
         lat = str_split_i(coords2, " ", 2) %>% 
           str_remove("\\)"), 
         year_pub = as.numeric(str_extract(source, "[:digit:]{4}")), # year of sampling
                year2 = year_pub-2, year =  coalesce(year, year2), 
                half_dec = round_any(year-1, 5, floor)) %>% 
  filter(year > 1960 & year <= 2020 & (gps_coords_precision < 5000 | is.na(gps_coords_precision)) & 
  ((project == "NDOP" & (!original_id %in% ndop.del$ID_ND_NALEZ | description == "Přijatý záznam"))
   | project != "NDOP")) %>% # keep records from 1961-2020, location uncertainty < 5 km, exclude unreliable data
  dplyr::select(record_id, name_lat, gps_coords_precision, half_dec, year, source, nalezce, 
                lat, long, project, locality, description) 

# convert to sf object to add information about grid cell
pladias.sf <- pladias %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Grid of Central European flora mapping (5'x3')
grid <- read_sf("shps/pladias_qq4_CZ_2017_WGS84.shp")

# assign species records to grid cells
plad.grid <- st_join(pladias.sf, grid) %>% 
  as.data.frame() 

# select gray-green grid cells and add variable for habitat mapping
plad.grid1 <- plad.grid %>% 
  dplyr::select(record_id, name_lat, half_dec, year, source, nalezce, 
                code, project, locality, description) %>% 
  semi_join(greygreen) %>% 
  mutate(hab.map = ifelse(project == "NDOP" & (str_detect(locality, "[:digit:]{6}") | 
                            str_detect(source, 'aktualizace|Aktualizace|Mapování biotopů|mapování biotopů')), 
                          1, 0) %>% replace_na(0))

# export to csv
write_csv(plad.grid1, "pladias_grid_230621.csv")

# nomenclature modification --------------------------------------------
nom_mod <- read_csv2("pladias_tax_names6.csv") %>% 
  distinct()

df1 <- plad.grid1 %>% 
  left_join(nom_mod) %>% 
  mutate(name_lat = coalesce(name_lat_mod, name_lat)) %>% 
  filter(name_lat != "XX" & name_lat != "XX_species_agg")

# replace NA in author with project/publication where available
df1$nalezce <- coalesce(df1$nalezce, df1$source)

# modify author names to surname
auth <- read_xlsx("author.xlsx") 

df2 <- df1 %>% 
  left_join(auth) %>% 
  filter(zmenit_na != "NA") %>%  # exclude records with NAs
  dplyr::rename("author" = "zmenit_na", "grid_small" = "code") %>% 
  dplyr::select(record_id, name_lat, year, half_dec, source, 
         grid_small, author, hab.map) %>% 
  distinct(grid_small, year, author, hab.map, name_lat, .keep_all = T) %>% # filter only one record per author*year*grid_cell*record origin
  mutate(grid_year_auth = str_c(grid_small, year, author, hab.map, sep = "_"))

# join number of species per author*grid cell*year*record origin
n.spec <- df2 %>% 
  group_by(grid_year_auth) %>% 
  dplyr::summarize(n = n())

df3 <- left_join(df2, n.spec) %>% 
  dplyr::rename("Nr_of_spec" = "n") %>% 
  dplyr::select(name_lat, grid_small, half_dec, grid_year_auth, hab.map, Nr_of_spec, author) %>% 
  filter(!is.na(grid_year_auth))

# exclude grid cells with only one visit
df3.nvis <- df3 %>% 
  group_by(grid_small) %>% 
  summarize(n.vis = length(unique(grid_year_auth))) %>% 
  filter(n.vis == 1)

# exclude species with less than 30 occurrences 
spe_30 <- df3 %>% 
  group_by(name_lat) %>% 
  summarize(nr_occ = length(unique(grid_year_auth))) %>% 
  filter(nr_occ >= 30)

# list of species with < 30 occurrences (Appendix S1)
spe_n_mod <- df3 %>% 
  group_by(name_lat) %>% 
  summarize(nr_occ = length(unique(grid_year_auth))) %>% 
  filter(nr_occ < 30) %>% 
  write_csv("species_n_mod.csv")

# species recorded only in one half-decade
spe_decade <- df3 %>% 
  group_by(name_lat) %>% 
  summarize(nr_decades = length(unique(half_dec))) %>% 
  filter(nr_decades > 1)

# filtering
df3 <- df3 %>% 
  semi_join(spe_30) %>% 
  semi_join(spe_decade) %>% 
  anti_join(df3.nvis)

# Appendix S2, numbers of records in each half-decade and numbers of decades with particular species' records
rec.plot <- df3 %>% 
  ggplot()+
  geom_bar(aes(half_dec), fill = 4, color = 1)+
  labs(y = "Number of records")+
  scale_x_continuous(breaks = seq(1920, 2010, 10), 
                     labels = c("1920s", "1930s", "1940s", "1950s", "1960s", 
                                "1970s", "1980s", "1990s", "2000s", "2010s"))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1100000), breaks = seq(0, 1000000, 250000))

dec.plot <- spe_decade %>% ggplot()+
  geom_bar(aes(nr_decades), color = 1, fill = 4)+
  labs(y = "Number of species", x = "Number of half-decades")+
  theme_bw()+
  scale_x_continuous(breaks = seq(1, 12, 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1700))

rec.plot + dec.plot + plot_annotation(tag_levels = "A")
ggsave("number_of_records_half_dec.png", width = 10, height = 4)


length(unique(df3$name_lat)) # number of species
length(unique(df3$grid_year_auth)) # number of visits
length(unique(df3$half_dec)) # number of half-decades
length(unique(df3$author)) # number of authors
length(unique(df3$grid_small)) # number of grid cells
length(unique(cz.grid$code)) # total number of grid cells in the Czech Republic

df3 <- df3 %>% 
  select(-author)

# export to csv
write_csv(df3, "pladias_occ_data_half_dec_20230621.csv")

# numbers of species records
freq <- df3 %>% 
  group_by(name_lat) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  mutate(name_lat = str_replace_all(name_lat, " ", "_"))

write_csv(freq, "pladias_frequency_20230621.csv")

# numbers of records per half-decade
rec_dec <- df3 %>% 
  group_by(name_lat, half_dec) %>% 
  count(name = "pres") %>% 
  write_csv("results/records_decades.csv")

# plot numbers of records per half-decade (Appendix S3)
cz.grid <- read_sf("shps/pladias_qq4_CZ_2017_WGS84.shp") %>% 
  left_join(df3, by = c("code" = "grid_small")) %>% 
  group_by(code, half_dec) %>% 
  summarize(n = n()) 

cz.grid1 <- cz.grid %>%
  as.data.frame() %>% 
  select(-geometry) %>%  
  pivot_wider(names_from = half_dec, values_from = n) %>% 
  select(-"NA") %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) %>% 
  pivot_longer(cols = (-code), names_to = "half_dec", values_to = "n")

cz.grid2 <- read_sf("shps/pladias_qq4_CZ_2017_WGS84.shp") %>% 
  left_join(cz.grid1) %>% 
  select(code, geometry, half_dec, n)

dens_map <- ggplot(cz.grid2)+
  theme_void()+
  scale_color_manual(values = 1)+
  geom_sf(mapping=aes(fill=n))+
  scale_fill_stepsn(breaks = c(0, 1, 50, 100, 500, 1000, 5000), 
                    colors = c("#ffea00","#ffba08", "#e85d04", "#9d0208"), 
                    na.value = "#ffea00", trans = "log1p")+
  theme(legend.position = "bottom",
        panel.background = element_blank(), plot.title = element_text(hjust = 0.05), 
        legend.key.width = unit(2, "cm"))+
  facet_wrap(~half_dec, ncol = 3) +
  labs(fill = "Number of records")

ggsave("map_half_decs_records.png", width = 8, height = 8)


# export csv for each species for occupancy model  ------------------------

df3$occ<-1

# modifying the variables to factors with 1:n sites etc.
plad.mod2<-setDF(dcast(setDT(df3), grid_year_auth+Nr_of_spec+half_dec+grid_small+hab.map
                     ~name_lat, fun=length))

plad.mod2$grid_year_auth <- as.integer(factor(plad.mod2$grid_year_auth))
plad.mod2$half_dec <- as.integer(factor(plad.mod2$half_dec))
plad.mod2$grid_small <- as.integer(factor(plad.mod2$grid_small))

specname <- colnames(plad.mod2)[6:ncol(plad.mod2)]

# write csv for each species
for(sp in specname){
  df <- plad.mod2 %>% 
    select(1:5, all_of(sp))
  colnames(df)[6] <- "value" 
  write_csv(df, paste0('species/', str_replace_all(sp, " ", "_"), '.csv'))
}
