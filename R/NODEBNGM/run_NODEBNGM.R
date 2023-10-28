# script to start m_main_general.r several times with different input files

## load packages
library(data.table)
library(tidyverse)

## load NODEBNGM functions
source("f_NODEBNGM.r")
source("m_main_general.r")

#
# ### example run
# data_name <- "3DLV_test_20-50"
#
# data_file <-
#   # "ALR_denom-7-other_ts_donorB_Genus_10_most_abundant_rel_counts.csv"
#   "ts_3DLV_20-50.csv"
#   # "ts_Ushio.csv"
#
# folderpath <-
#   # "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/Python/ALR_transformation/ALR_transformed_data/"
#   "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/final_datasets/"
#
# run_NODEBNGM(data_name, data_file, folderpath)


### list possible data files

synthetic_data <-
  list(
    c('miaSim', 'miaSim_GLV_4species_new.csv'),
    c('miaSim_noise_0-005', 'ts_miaSim_GLV_4species_new_noise_0-005.csv'),
    c('miaSim_noise_0-01', 'ts_miaSim_GLV_4species_new_noise_0-01.csv'),
    c('miaSim_noise_0-02', 'ts_miaSim_GLV_4species_new_noise_0-02.csv'),
    c('miaSim_noise_0-04', 'ts_miaSim_GLV_4species_new_noise_0-04.csv'),
    c('3DLV', 'ts_3DLV.csv'),
    c('VanderPol', 'ts_VanderPol.csv'),
    c('VanderPol_noise_0-1', 'ts_VanderPol_noise_0-1.csv'),
    c('VanderPol_noise_0-2', 'ts_VanderPol_noise_0-2.csv'),
    c('VanderPol_noise_0-5', 'ts_VanderPol_noise_0-5.csv'),
    c('VanderPol_noise_1', 'ts_VanderPol_noise_1.csv')
  )

compositional_data <-
  list(
    c('donorA_ALR', 'ALR_denom-8-other_ts_donorA_Genus_10_most_abundant_rel_counts.csv'),
    c('donorB_ALR', 'ALR_denom-7-other_ts_donorB_Genus_10_most_abundant_rel_counts.csv'),
    c('male_ALR', 'ALR_denom-5-other_ts_male_Genus_10_most_abundant_rel_counts.csv'),
    c('female_ALR', 'ALR_denom-10-other_ts_female_Genus_10_most_abundant_rel_counts.csv'),
    c('Silverman_all_ALR', 'ALR_denom-7-other_ts_Silverman_Vall_all_Genus_10_most_abundant_rel_counts.csv'),
    c('Silverman_daily_ALR', 'ALR_denom-7-other_ts_Silverman_Vall_daily_Genus_10_most_abundant_rel_counts.csv'),
    c('Silverman_hourly_ALR', 'ALR_denom-7-other_ts_Silverman_Vall_hourly_Genus_10_most_abundant_rel_counts.csv'),
    c('Bucci_ALR', 'ALR_denom-3-Clostridium-scindens_ts_bucci_subject_all_rel_counts_denoised.csv')
  )

other_ts_data <-
  list(
    c('BioTIME_study_339_Genus_10', 'ts_study_339_Genus_10_most_abundant.csv'),
    c('BioTIME_study_339_Species_15', 'ts_study_339_Species_15_most_abundant.csv'),
    c('BioTIME_study_363_Genus_10', 'ts_study_363_Genus_10_most_abundant.csv'),
    c('BioTIME_study_363_Species_15', 'ts_study_363_Species_15_most_abundant.csv'),
    c('BioTIME_study_39_Genus_10', 'ts_study_39_Genus_10_most_abundant.csv'),
    c('BioTIME_study_39_Species_15', 'ts_study_39_Species_15_most_abundant.csv'),
    c('BioTIME_study_478_Genus_10', 'ts_study_478_Genus_10_most_abundant.csv'),
    c('BioTIME_study_478_Species_15', 'ts_study_478_Species_15_most_abundant.csv'),
    c('Ushio', 'ts_Ushio.csv')
  )

all_files <- c(synthetic_data, other_ts_data)
ALR_files <- compositional_data


### start runs

for(run in 2:7) {
  for(dataset in synthetic_data) {
    folderpath <-
      "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/explore/data/final_datasets/"
    data_name <- dataset[1]
    data_file <- dataset[2]

    run_NODEBNGM(data_name, data_file, folderpath, run)
  }

  for(dataset in ALR_files){
    folderpath <-
      "C:/Users/Maria/Documents/Masterstudium/Masterarbeit/MScThesis/Python/ALR_transformation/ALR_transformed_data/"
    data_name <- dataset[1]
    data_file <- dataset[2]

    run_NODEBNGM(data_name, data_file, folderpath, run)
  }
}

