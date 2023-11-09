01f-timeseries-NODEBNGMpaper
================
Compiled at 2023-11-09 10:46:23 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "9506fa44-f8a3-401c-aa4c-950659e05f3f")
```

The purpose of this document is …

``` r
library("conflicted")
library(data.table)
library(dplyr)
library(ggplot2)
library(phyloseq)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Read Data and make Phyloseq objects

``` r
folderpath_data <-
  "input_data/NODEBNGM/"

files <- c("3DLV",
           "AFR1",
           "AFR2",
           "AFR3",
           "RPS")

for(ts_name in files) {
  tmp <-
    fread(paste0(folderpath_data, "TS_", ts_name, ".csv"), header = T) %>%
    .[, SampleID := sprintf("ID-%03d", t)] %>%
    tibble::column_to_rownames("SampleID") %>%
    rename(Time = t)
  
  assign(paste0("dt_", ts_name),
         tmp)
}

dt_HL <-
  fread(paste0(folderpath_data, "TS_HL.csv"), header = T) %>%
    .[, SampleID := sprintf("ID-%03d", as.numeric(row.names(.)))] %>%
    tibble::column_to_rownames("SampleID") %>%
    rename(Time = Year)

files <- c(files, "HL")

# convert to phyloseq after reading the data
for(ts_name in files) {
  tmp_otu <-
    get(paste0("dt_", ts_name)) %>% 
    subset(select = -Time) %>% 
    t() %>% 
    otu_table(taxa_are_rows = T)
  
  tmp_sample <-
    get(paste0("dt_", ts_name)) %>% 
    subset(select = Time) %>% 
    sample_data()

  assign(paste0("ps_", ts_name),
         phyloseq(tmp_otu, tmp_sample))
}

rm(tmp, tmp_otu, tmp_sample)
```

### Ushio Dataset

``` r
# special case (including more sample info): Ushio
dt_Ushio_raw <-
  fread(paste0(folderpath_data, "TS_Ushio.csv")) %>%
  .[, SampleID := sprintf("ID-%03d", time_step)] %>%
  tibble::column_to_rownames("SampleID")

# otu table
dt_Ushio <-
  dt_Ushio_raw %>% select(
    "Aurelia.sp",
    "Engraulis.japonicus",
    "Plotosus.lineatus",
    "Sebastes.inermis",
    "Trachurus.japonicus",
    "Girella.punctata",
    "Pseudolabrus.sieboldi",
    "Halichoeres.poecilopterus",
    "Halichoeres.tenuispinnis",
    "Chaenogobius.gulosus",
    "Pterogobius.zonoleucus",
    "Tridentiger.trigonocephalus",
    "Siganus.fuscescens",
    "Sphyraena.pinguis",
    "Rudarius.ercodes"
  )

otu_Ushio <-
  copy(dt_Ushio) %>% 
  t() %>% 
  otu_table(taxa_are_rows = T)

# sample info
samples_Ushio <-
  dt_Ushio_raw %>%
  select("date_tag", "surf.t", "bot.t", "Y", "M", "D", "time_step")
samples_Ushio$Day <-
  as.Date(paste0(samples_Ushio$Y, "-", samples_Ushio$M, "-", 
                 samples_Ushio$D), 
          format = "%Y-%m-%d")
# Mutate the Time column to be the integer day starting with 1
samples_Ushio <- samples_Ushio %>% 
  mutate(Time = as.numeric(difftime(Day, min(Day), units = "days")) + 1)

# taxonomic table
tax_Ushio <- 
  fread(paste0(folderpath_data, "tax_table_Ushio.csv")) %>% 
  as.matrix()
rownames(tax_Ushio) <- tax_Ushio[, "Species"]


# make phyloseq object
ps_Ushio <-
  phyloseq(otu_Ushio, 
           tax_table(tax_Ushio),
           sample_data(samples_Ushio))

# add Ushio to list of files
files <- c(files, "Ushio")
```

### 3DLV Data

Only keep values in the interval \[0, 60\].

``` r
ps_3DLV <- subset_samples(ps_3DLV, Time <= 60)
```

### Overview over the phyloseq objects

``` r
for(ts_name in files) {
  cat(ts_name, "\n")
  print(get(paste0("ps_", ts_name)))
  cat("\n")
}
```

    ## 3DLV 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3 taxa and 61 samples ]
    ## sample_data() Sample Data:       [ 61 samples by 1 sample variables ]
    ## 
    ## AFR1 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3 taxa and 67 samples ]
    ## sample_data() Sample Data:       [ 67 samples by 1 sample variables ]
    ## 
    ## AFR2 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3 taxa and 67 samples ]
    ## sample_data() Sample Data:       [ 67 samples by 1 sample variables ]
    ## 
    ## AFR3 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3 taxa and 41 samples ]
    ## sample_data() Sample Data:       [ 41 samples by 1 sample variables ]
    ## 
    ## RPS 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3 taxa and 648 samples ]
    ## sample_data() Sample Data:       [ 648 samples by 1 sample variables ]
    ## 
    ## HL 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2 taxa and 91 samples ]
    ## sample_data() Sample Data:       [ 91 samples by 1 sample variables ]
    ## 
    ## Ushio 
    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 15 taxa and 285 samples ]
    ## sample_data() Sample Data:       [ 285 samples by 9 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]

## Plot Datasets

``` r
for(ts_name in files) {
  ps_tmp <-
      psmelt(get(paste0("ps_", ts_name)))
  plt_tmp <-
    ggplot(ps_tmp, aes(x = Time, y = Abundance, col = OTU)) +
    geom_line() +
    labs(title = ts_name,
         subtitle = paste0("number of time points: ", uniqueN(ps_tmp[,1])),
         x = "Time")
  print(plt_tmp)
}
```

![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->

## Get relative abundances

### for TS_RPS

``` r
# calculcate relative abundances
ps_RPS_rel_counts <-
  transform_sample_counts(ps_RPS, function(x) x / sum(x))

# bar plot phyloseq objects
plot_bar(ps_RPS_rel_counts, x = "Time", fill = "OTU") +
  geom_bar(aes(color = OTU, fill = OTU), stat = "identity", position = "stack") +
  labs(title = "Ushio relative Abundances",
       x = "Time")
```

![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Save Phyloseq Objects

``` r
for (ts_name in c(files, "RPS_rel_counts")) {
  saveRDS(get(paste0("ps_", ts_name)),
          path_target(paste0("ps_", ts_name, ".rds")))
}
```

## Save TS as csv file

``` r
for (ts_name in c(files, "RPS_rel_counts")) {
  
  # get the tmp phyloseq object
  ps_obj <- get(paste0("ps_", ts_name))
  
  if(taxa_are_rows(ps_obj)) {
    otu_tmp <- t(otu_table(ps_obj))
  } else {
    otu_tmp <- otu_table(ps_obj)
  }
  # combine count data with time information
  ts_obj <-
    cbind(sample_data(ps_obj)[, "Time"],
          otu_tmp)
  # print(head(ts_obj))
  
  # save time series as csv file
  write.csv(
    ts_obj,
    path_target(paste0("ts_", ts_name, ".csv")),
    row.names = F
  )
    
}
```

## Files written

These files have been written to the target directory,
`data/01f-timeseries-NODEBNGMpaper`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 16 × 4
    ##    path                  type         size modification_time  
    ##    <fs::path>            <fct> <fs::bytes> <dttm>             
    ##  1 ps_3DLV.rds           file        2.04K 2023-11-09 10:46:31
    ##  2 ps_AFR1.rds           file         2.2K 2023-11-09 10:46:31
    ##  3 ps_AFR2.rds           file        1.63K 2023-11-09 10:46:31
    ##  4 ps_AFR3.rds           file        1.49K 2023-11-09 10:46:31
    ##  5 ps_HL.rds             file        1.44K 2023-11-09 10:46:31
    ##  6 ps_RPS.rds            file        7.99K 2023-11-09 10:46:31
    ##  7 ps_RPS_rel_counts.rds file       18.22K 2023-11-09 10:46:31
    ##  8 ps_Ushio.rds          file       12.23K 2023-11-09 10:46:31
    ##  9 ts_3DLV.csv           file        3.41K 2023-11-09 10:46:31
    ## 10 ts_AFR1.csv           file        3.91K 2023-11-09 10:46:31
    ## 11 ts_AFR2.csv           file        3.81K 2023-11-09 10:46:31
    ## 12 ts_AFR3.csv           file        2.37K 2023-11-09 10:46:31
    ## 13 ts_HL.csv             file        1.56K 2023-11-09 10:46:31
    ## 14 ts_RPS.csv            file       12.57K 2023-11-09 10:46:31
    ## 15 ts_RPS_rel_counts.csv file       36.97K 2023-11-09 10:46:31
    ## 16 ts_Ushio.csv          file       13.14K 2023-11-09 10:46:31
