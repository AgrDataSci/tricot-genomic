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
G <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)

#.................................................
#.................................................
# Set parameters ####
keep <- grepl("NT_", names(ind))

data <- cbind(G, ind[keep])

n <- nrow(data)

fit <- pltree(G ~ ., 
              data = data, 
              npseudo = 15)



items <- dimnames(coef(fit))[[2]]

coeffs <- as_tibble(t(coef(fit)))
names(coeffs) <- paste0("node",1:ncol(coeffs))
coeffs <- bind_cols(items = items, coeffs)

plot(coeffs$node1)
plot(coeffs$node2)

# # Check if elevation play a role in the performance of varieties in Ethiopia ####
# nodes <- partykit::nodeids(tree, terminal = TRUE)
# 
# #get the estimates from nodes
# node <- list()
# 
# for(i in seq_along(nodes)){
#   node[[i]] <- tree[[ nodes[i] ]]$node$info$object %>%
#     itempar(., vcov = TRUE, alias = FALSE) %>%
#     qvcalc::qvcalc.itempar(.) %>%
#     extract2(2) %>%
#     as.matrix() %>%
#     as_tibble() %>%
#     mutate(items = items)
# }
# 
# node <- bind_cols(items = node[[1]]$items, 
#                   n1 = node[[1]]$estimate, 
#                   n2 = node[[2]]$estimate )

#read the passport data
pass <- read_csv("data/passport_data_durumwheat.csv")

pass %<>%
  mutate(items = genotype) %>%
  inner_join(. , coeffs, by = "items") %>%
  dplyr::select(. , items, lon , lat , node1, node2, adm1) %>% 
  filter(!is.na(lon))

# get the elevation from origin of items
# read the elevation data
# source http://gisweb.ciat.cgiar.org/TRMM/SRTM_Resampled_250m/

r <- stack("data/SRTM_NE_250m.tif")

e <- extract(r, pass[c("lon","lat")], buffer = 500)

names(e) <- paste0(pass$items,".")

e <- unlist(e, use.names = TRUE)

e <- as.data.frame(cbind(items = str_split(names(e), pattern = "[.]", simplify = T)[,1], elev = e))

rownames(e) <- 1:nrow(e)

#add to the main data
pass %>%
  inner_join(., e , by = "items", all.y = TRUE) %>%
  mutate(dif = node1 - node2,
         group = ifelse(node1 >= mean(node1) & node2 <= mean(node2), 1,
                        ifelse(node1 <= mean(node1) & node2 >= mean(node2), 2, 1)),
         elev = as.numeric(as.character(elev))) ->
  df

ttest <- with(df,
              t.test(elev ~ group))
ttest
