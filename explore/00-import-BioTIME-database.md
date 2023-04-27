00-import-BioTIME-database
================
Compiled at 2023-04-27 08:12:24 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "92e0b108-b72d-4d77-bf7a-4e46e4cde787")
```

The purpose of this document is to import and have a first look into the
BioTIME database from the paper *BioTIME: A database of biodiversity
time series for the Anthropocene*.

``` r
library("conflicted")
library(dplyr)
library(data.table)
library(ggplot2)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Data

#### Read files

``` r
#### read the file
dt_fullquery <- 
  fread("data/BioTIME/BioTIMEquery_24_06_2021.csv")

# focus on one study
study_id = 63
```

#### Focus on only one study (ID = 63)

``` r
# extract one study
dt_query <- 
  dt_fullquery[STUDY_ID == study_id]

knitr::kable(head(dt_query))
```

|    V1 | STUDY_ID | DAY | MONTH | YEAR | SAMPLE_DESC                                     | PLOT | ID_SPECIES | LATITUDE | LONGITUDE | sum.allrawdata.ABUNDANCE | sum.allrawdata.BIOMASS | GENUS     | SPECIES        | GENUS_SPECIES            |
|------:|---------:|----:|------:|-----:|:------------------------------------------------|:-----|-----------:|---------:|----------:|-------------------------:|-----------------------:|:----------|:---------------|:-------------------------|
| 78181 |       63 |  NA |    NA | 1959 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       1414 |  52.4221 | -0.180928 |                      130 |                      0 | Libellula | quadrimaculata | Libellula quadrimaculata |
| 78182 |       63 |  NA |    NA | 1961 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       1414 |  52.4221 | -0.180928 |                       45 |                      0 | Libellula | quadrimaculata | Libellula quadrimaculata |
| 78183 |       63 |  NA |    NA | 1962 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       1309 |  52.4221 | -0.180928 |                       12 |                      0 | Ischnura  | elegans        | Ischnura elegans         |
| 78184 |       63 |  NA |    NA | 1962 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       1414 |  52.4221 | -0.180928 |                       95 |                      0 | Libellula | quadrimaculata | Libellula quadrimaculata |
| 78185 |       63 |  NA |    NA | 1962 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       2258 |  52.4221 | -0.180928 |                       19 |                      0 | Sympetrum | striolatum     | Sympetrum striolatum     |
| 78186 |       63 |  NA |    NA | 1963 | 52.422100\_-0.180928_WoodwaltonFen_20smallPonds | NA   |       1309 |  52.4221 | -0.180928 |                       26 |                      0 | Ischnura  | elegans        | Ischnura elegans         |

### plot time series for selected study

``` r
# ## group species
# rf_grouped <- as.data.table(rf) %>% 
#   # .[, c("Genus", "Species") := tstrsplit(Species, " ")] %>% 
#   .[, .(Abundance_sum = sum(Abundance)),
#         by = c("Year", "Species")]

ggplot(rf, aes(Year, Abundance, col = Species)) +
  geom_line() #+
```

<img src="00-import-BioTIME-database_files/figure-gfm/unnamed-chunk-5-1.png" width="100%" />

``` r
  # theme(legend.position = "none")
```

### values of selected time series

``` r
rf_time <- as.data.table(rf) %>% 
  # .[, c("Genus", "Species") := tstrsplit(Species, " ")] %>% 
  dcast(., Species ~ Year, value.var = "Abundance")

knitr::kable(rf_time)
```

| Species                  | 1959 | 1961 | 1962 | 1963 | 1964 | 1965 | 1966 | 1967 | 1968 | 1969 | 1970 | 1971 | 1972 | 1973 | 1974 | 1975 | 1976 | 1977 | 1978 | 1979 | 1980 | 1981 | 1982 | 1983 | 1984 | 1985 | 1986 | 1987 | 1988 |
|:-------------------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| Coenagrion puella        |   NA |   NA |   NA |    1 |    9 |   14 |   16 |    3 |    1 |    4 |   16 |   30 |   14 |   12 |   19 |   51 |   11 |    7 |  137 |  204 |   67 |  197 |  183 |  149 |  102 |  141 |  147 |  178 |   52 |
| Ischnura elegans         |   NA |   NA |   12 |   26 |   41 |   14 |   30 |   33 |   15 |   25 |   18 |   38 |    7 |   25 |   10 |   25 |   13 |    9 |   30 |   27 |   22 |   27 |   43 |   45 |   32 |   22 |   15 |   12 |   16 |
| Lestes sponsa            |   NA |   NA |   NA |   NA |    1 |    1 |   36 |   16 |   57 |   63 |   65 |  111 |   17 |   41 |   91 |   21 |   15 |   61 |   42 |   31 |    9 |  155 |  188 |   62 |   89 |   66 |   93 |   69 |   92 |
| Libellula quadrimaculata |  130 |   45 |   95 |  210 |   13 |    2 |    5 |    1 |    2 |    3 |    9 |    2 |   NA |    1 |   NA |    2 |    4 |    4 |    6 |   13 |   16 |   20 |   17 |   13 |   15 |   12 |   19 |   13 |   14 |
| Sympetrum striolatum     |   NA |   NA |   19 |   18 |   16 |   12 |   16 |   12 |   17 |    7 |    8 |    4 |    4 |    9 |    8 |    2 |    6 |    9 |   11 |    8 |    4 |    4 |    6 |    6 |    4 |   11 |    7 |   12 |   19 |

## BioTime Metadata

#### Read file

``` r
#### read the file
dt_biotimeMeta <-
  fread("data/BioTIME/biotimeMetadata_24_06_2021.csv")
```

<!-- #### List of all columns -->
<!-- ```{r, results = "asis", echo = FALSE} -->
<!-- # list all column names of Metadata -->
<!-- cat(paste("*", colnames(dt_biotimeMeta)), sep="\n") -->
<!-- ``` -->

#### First rows of table (only including important columns)

Select some important columns to get a better overview which studies are
included in the dataset.

Some of the Studies in the BioTIME data are listed in the table below
(only studies with at least 10 datapoints):

``` r
knitr::kable(dt_biotimeMeta_small %>% head(n = 20))
```

| STUDY_ID | REALM       | CLIMATE         | TAXA                      | DATA_POINTS | START_YEAR | END_YEAR | NUMBER_OF_SPECIES | NUMBER_OF_SAMPLES | TOTAL | ABUNDANCE_TYPE   | BIOMASS_TYPE |
|---------:|:------------|:----------------|:--------------------------|------------:|-----------:|---------:|------------------:|------------------:|------:|:-----------------|:-------------|
|       18 | Terrestrial | Temperate       | Terrestrial plants        |          29 |       1923 |     1973 |                98 |               542 |  8034 | Count            | NA           |
|       33 | Marine      | Temperate       | Marine plants             |          18 |       1992 |     2009 |               170 |                 1 | 17841 | Count            | NA           |
|       39 | Terrestrial | Temperate       | Birds                     |          45 |       1970 |     2015 |                52 |                45 |   959 | Density          | NA           |
|       41 | Terrestrial | Temperate       | Birds                     |          10 |       1923 |     1940 |                56 |                10 |   418 | Count            | NA           |
|       45 | Marine      | Tropical        | Fish                      |          10 |       2006 |     2015 |               338 |              1105 | 47282 | Count            | Size         |
|       46 | Terrestrial | Temperate       | Birds                     |          47 |       1928 |     1979 |                29 |                47 |   528 | Count            | NA           |
|       47 | Terrestrial | Temperate       | Birds                     |          26 |       1952 |     1977 |                13 |                35 |   392 | Count            | NA           |
|       51 | Terrestrial | Temperate       | Birds                     |          14 |       1964 |     1977 |                 4 |                 1 |    56 | MeanCount        | NA           |
|       52 | Terrestrial | Polar           | Mammals                   |          13 |       1968 |     1980 |                 3 |                 5 |   144 | MeanCount        | NA           |
|       53 | Terrestrial | Temperate       | Mammals                   |          10 |       1966 |     1976 |                 5 |                 1 |    44 | MeanCount        | NA           |
|       54 | Terrestrial | Tropical        | Terrestrial invertebrates |          24 |       1991 |     2014 |                19 |               798 | 21702 | Count            | NA           |
|       56 | Terrestrial | Temperate       | Mammals                   |          20 |       1989 |     2008 |                28 |              8026 | 16657 | Count            | Weight       |
|       57 | Freshwater  | Temperate       | Fish                      |          32 |       1981 |     2012 |                76 |               258 | 10892 | Count            | NA           |
|       58 | Terrestrial | Tropical        | Birds                     |          18 |       1991 |     2008 |                31 |                 6 |  1171 | Count            | NA           |
|       59 | Terrestrial | Temperate       | Mammals                   |          26 |       1977 |     2002 |                29 |                 2 |   427 | Count            | NA           |
|       63 | Terrestrial | Temperate       | Terrestrial invertebrates |          29 |       1959 |     1988 |                 5 |                 1 |   132 | Count            | NA           |
|       67 | Terrestrial | Temperate       | Birds                     |          24 |       1983 |     2006 |                68 |               417 | 15448 | Count            | NA           |
|       68 | Marine      | Temperate       | Marine invertebrates      |          13 |       1988 |     2000 |                15 |              2569 | 38535 | Presence/Absence | NA           |
|       69 | Marine      | Temperate       | Birds                     |          10 |       1994 |     2003 |                27 |              9864 | 22694 | Presence/Absence | NA           |
|       71 | Marine      | Polar/Temperate | Marine plants             |          10 |       1993 |     2003 |               412 |              3250 | 35323 | Count            | NA           |

### Detailed Info about some columns

##### DATA_POINTS

Number of years covered by study

``` r
ggplot(dt_biotimeMeta, aes(DATA_POINTS)) +
  geom_histogram(binwidth = 5) +
  scale_x_continuous(breaks = seq(0,100,10))
```

![](00-import-BioTIME-database_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

##### AB_BIO

Does study contain abundance, biomass or both (A, B or AB)

``` r
knitr::kable(dt_biotimeMeta$AB_BIO %>% table(useNA = "ifany"))
```

| .   | Freq |
|:----|-----:|
|     |    1 |
| A   |  270 |
| AB  |   64 |
| B   |   46 |

##### ABUNDANCE_TYPE

``` r
knitr::kable(dt_biotimeMeta$ABUNDANCE_TYPE %>% table(useNA = "ifany"))
```

| .                | Freq |
|:-----------------|-----:|
| Count            |  262 |
| Density          |   40 |
| MeanCount        |    7 |
| Presence/Absence |   28 |
| NA               |   44 |

##### BIOMASS_TYPE

``` r
knitr::kable(dt_biotimeMeta$BIOMASS_TYPE %>% table(useNA = "ifany"))
```

| .      | Freq |
|:-------|-----:|
| Cover  |   34 |
| Size   |    3 |
| Volume |    3 |
| Weight |   68 |
| NA     |  273 |

##### REALM

``` r
dt_biotimeMeta$REALM %>% unique()
```

    ## [1] "Terrestrial" "Marine"      "Freshwater"

##### CLIMATE

``` r
dt_biotimeMeta$CLIMATE %>% unique()
```

    ## [1] "Temperate"          "Tropical"           "Polar"             
    ## [4] "Polar/Temperate"    "Temperate/Tropical" "Global"

##### TAXA

``` r
dt_biotimeMeta$TAXA %>% unique()
```

    ##  [1] "Terrestrial plants"        "Marine plants"            
    ##  [3] "Birds"                     "Fish"                     
    ##  [5] "Mammals"                   "Terrestrial invertebrates"
    ##  [7] "Marine invertebrates"      "Benthos"                  
    ##  [9] "All"                       "Freshwater invertebrates" 
    ## [11] "Freshwater plants"         "Reptiles"                 
    ## [13] "Amphibians"                "Fungi"

##### NUMBER_OF_SPECIES

Number of species in study

``` r
range(dt_biotimeMeta$NUMBER_OF_SPECIES)
```

    ## [1]    3 5624

``` r
# ggplot(dt_biotimeMeta, aes(NUMBER_OF_SPECIES)) +
#   geom_histogram(binwidth = 10)

ggplot(dt_biotimeMeta[NUMBER_OF_SPECIES < 200], aes(NUMBER_OF_SPECIES)) +
  geom_histogram(binwidth = 1)
```

![](00-import-BioTIME-database_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

##### NUMBER_OF_SAMPLES

Number of unique samples in study

``` r
dt_biotimeMeta$NUMBER_OF_SAMPLES %>% range()
```

    ## [1]      1 171843

``` r
# ggplot(dt_biotimeMeta, aes(NUMBER_OF_SAMPLES)) +
#   geom_histogram(binwidth = 100)

ggplot(dt_biotimeMeta[NUMBER_OF_SAMPLES < 100], aes(NUMBER_OF_SAMPLES)) +
  geom_histogram(binwidth = 1)
```

![](00-import-BioTIME-database_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

##### NUMBER_LAT_LONG

``` r
dt_biotimeMeta$NUMBER_LAT_LONG %>% range()
```

    ## [1]     1 93906

``` r
knitr::kable(dt_biotimeMeta$NUMBER_LAT_LONG %>% table() %>% head(n = 10)) # mostly 1
```

| .   | Freq |
|:----|-----:|
| 1   |  170 |
| 2   |   14 |
| 3   |    4 |
| 4   |    5 |
| 5   |    2 |
| 6   |    7 |
| 7   |    5 |
| 8   |    3 |
| 9   |    3 |
| 10  |    5 |

## Files written

These files have been written to the target directory,
`data/00-import-BioTIME-database`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: path <fs::path>, type <fct>, size <fs::bytes>,
    ## #   modification_time <dttm>
