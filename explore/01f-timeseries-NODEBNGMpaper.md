01f-timeseries-NODEBNGMpaper
================
Compiled at 2023-09-28 10:04:33 UTC

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
  "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/Literatur/Code/NODEBNGM/data/"

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

``` r
# special case (including more sample info): Ushio
dt_Ushio_raw <-
  fread(paste0(folderpath_data, "TS_Ushio.csv")) %>%
  .[, SampleID := sprintf("ID-%03d", time_step)] %>%
  tibble::column_to_rownames("SampleID")

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

samples_Ushio <-
  dt_Ushio_raw %>%
  select("date_tag", "surf.t", "bot.t", "Y", "M", "D", "time_step")
samples_Ushio$Time <-
  as.Date(paste0(samples_Ushio$Y, "-", samples_Ushio$M, "-", 
                 samples_Ushio$D), 
          format = "%Y-%m-%d")

ps_Ushio <-
  phyloseq(otu_Ushio, sample_data(samples_Ushio))

files <- c(files, "Ushio")
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
    ## otu_table()   OTU Table:         [ 3 taxa and 101 samples ]
    ## sample_data() Sample Data:       [ 101 samples by 1 sample variables ]
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
    ## sample_data() Sample Data:       [ 285 samples by 8 sample variables ]

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

![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->

## Get relative abundances

### for TS_RPS

``` r
# calculcate relative abundances
ps_RPS_rel <-
  transform_sample_counts(ps_RPS, function(x) x / sum(x))

# bar plot phyloseq objects
plot_bar(ps_RPS_rel, x = "Time", fill = "OTU") +
  geom_bar(aes(color = OTU, fill = OTU), stat = "identity", position = "stack") +
  labs(title = "Ushio relative Abundances",
       x = "Time")
```

![](01f-timeseries-NODEBNGMpaper_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Save Phyloseq Objects

``` r
for (ts_name in files) {
  saveRDS(get(paste0("ps_", ts_name)),
          path_target(paste0("ps_", ts_name, ".rds")))
}

# save relative counts for RPS
saveRDS(ps_RPS_rel,
        path_target("ps_RPS_rel_counts.rds"))
```

## Files written

These files have been written to the target directory,
`data/01f-timeseries-NODEBNGMpaper`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 8 × 4
    ##   path                  type         size modification_time  
    ##   <fs::path>            <fct> <fs::bytes> <dttm>             
    ## 1 ps_3DLV.rds           file        3.12K 2023-09-28 10:04:42
    ## 2 ps_AFR1.rds           file         2.2K 2023-09-28 10:04:42
    ## 3 ps_AFR2.rds           file        1.63K 2023-09-28 10:04:42
    ## 4 ps_AFR3.rds           file        1.49K 2023-09-28 10:04:42
    ## 5 ps_HL.rds             file        1.44K 2023-09-28 10:04:42
    ## 6 ps_RPS.rds            file        7.99K 2023-09-28 10:04:42
    ## 7 ps_RPS_rel_counts.rds file       18.22K 2023-09-28 10:04:42
    ## 8 ps_Ushio.rds          file       11.06K 2023-09-28 10:04:42
