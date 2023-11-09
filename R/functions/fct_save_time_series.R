# function to save time series data (phyloseq as csv file) for further analysis with
# * deepmod
# * clv
# * NODEBNGM

save_ps_as_csv <- function(ps, filepath, t_name = "Time") {
  # combine count data with time information
  ts_obj <-
    cbind(sample_data(ps)[, t_name],
          t(otu_table(ps)))

  # save time series as csv file
  write.csv(
    ts_obj,
    filepath,
    row.names = F
  )
}
