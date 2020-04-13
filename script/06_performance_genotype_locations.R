library("tidyverse")
library("magrittr")
library("raster")
library("gosset")
library("PlackettLuce")
library("partykit")
library("qvcalc")


#...................................
#...................................
# read data 
list.files("data")

# farmer rankings
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()


# environmental indices
ind <- "data/environmental_indices.csv"
ind %<>% 
  read.csv()

#.....................................
#.....................................
# create PlackettLuce rankings ####
G <- rank_PL(data = df,
             items = "genotype",
             input = "farmer_rank",
             id = "id",
             grouped = TRUE)

#.................................................
#.................................................
# Set parameters ####
keep <- grepl("NT_|DT_", names(ind))

data <- cbind(G, ind[keep])

n <- nrow(data)

fit <- pltree(G ~  minNT_veg + maxNT_rep,
              data = data,
              npseudo = 5,
              gamma = TRUE, 
              alpha = 0.01)
plot(fit)

items <- dimnames(coef(fit))[[2]]

coeffs <- as_tibble(t(itempar(fit)))
names(coeffs) <- paste0("n",1:ncol(coeffs))
coeffs <- bind_cols(items = items, coeffs)

# ..................................
# ..................................
# read the passport data ####
pass <- read_csv("data/passport_data_durumwheat.csv")

pass %<>%
  mutate(items = genotype) %>%
  inner_join(. , coeffs, by = "items") %>%
  dplyr::select(. , items, lon , lat , n1, n2, adm1) %>% 
  filter(!is.na(lon))

# get the elevation from origin of items
# read the elevation data
# source http://gisweb.ciat.cgiar.org/TRMM/SRTM_Resampled_250m/

r <- stack("data/SRTM_NE_250m.tif")

e <- extract(r, pass[c("lon","lat")], buffer = 500)

names(e) <- paste0(pass$items,".")

e <- unlist(e, use.names = TRUE)

e <- bind_cols(items = str_split(names(e), pattern = "[.]", 
                                 simplify = TRUE)[,1], 
                         elev = e)
# ..............................
# ..............................
# create groups and calculate the difference ####

pass %>% 
  mutate(group = ifelse(n1 >= mean(n1), 1, 2),
         group = ifelse(n2 >= mean(n2), 2, 1)) -> 
  df

df

# add elevation data
df %<>% 
  inner_join(., e , by = "items", all.y = TRUE)

str(df)

ttest <- with(df,
              t.test(elev ~ group))
ttest


# ..............................
# ..............................
# export outputs ####
output <- "output/heat_tolerance/"
dir.create(output, 
           recursive = TRUE,
           showWarnings = FALSE)

capture.output(ttest,
               file = paste0(output, "t-test_genotype_groups.txt"))


capture.output(fit,
               summary(fit),
               file = paste0(output, "PL-fit_genotype_groups.txt"))



df %<>% 
  rename(genotype = items) %>% 
  group_by(genotype) %>% 
  summarise(elev = mean(elev),
            lon = unique(lon),
            lat = unique(lat),
            group = unique(group)) 

# add missing genotypes

x <- bind_cols(genotype = items[!items %in% df$genotype],
               elev = rep(NA,4),
               lon = rep(NA,4),
               lat = rep(NA,4),
               group = rep(NA,4))

df <- bind_rows(df, x)

write_csv(df, paste0(output, "genotype_groups.csv"))


