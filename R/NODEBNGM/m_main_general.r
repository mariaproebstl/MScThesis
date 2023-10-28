##################
## m_main_AFR.r ##
##################

## goal:
## - use Bayesian neural gradient matching (BNGM) to fit neural ordinary differential equation model (NODE) to time series
## - analyse NODEs to derive interactions and contributions between variables (sensitivities, Geber method)

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

#
###

##############
## INITIATE ##
##############

## goal: load data, functions

run_NODEBNGM <- function(data_name, data_file, folderpath, run) {

## make output directory
pathToOut <<- paste0("output/out_", data_name,
                     "/out_", data_name, "_run", sprintf("%02d", run))
system(paste("mkdir",pathToOut))

#
###

#################
## TIME SERIES ##
#################

## goal: load and format time series

filename <<- data_file
  # "ALR_denom-7-other_ts_donorB_Genus_10_most_abundant_rel_counts.csv"
  # "miaSim_GLV_4species_new.csv"
  # "ts_Ushio.csv"

## load data
TS <<- read.table(paste0(folderpath, filename), sep=",",header=T)
colnames(TS)[1] <<- "t"
TS$t <<- order(TS$t)

n_taxa <<- ncol(TS) - 1

## extract time steps and columns of interest
selected_time_steps <<- 1:nrow(TS)
selected_columns  <<- colnames(TS)
TS <<- TS[selected_time_steps, selected_columns]

## normalise time series
normalise <<- function(x) (x-min(x))/(max(x)-min(x)) * 10
TS[,-1] <<- apply(TS[,-1], 2, normalise)

## set 0s to small value to avoid NAs
for(i in 2:ncol(TS)) TS[,i][which(TS[,i] < 0.005)] <<- 0.005

## save TS as csv
write.csv(
  TS,
  paste0(pathToOut, "/TS_", data_name, ".csv"),
  row.names = F
)

## visualise time series
# pdf(paste(pathToOut,"/fig_time_series.pdf",sep=""))
png(paste(pathToOut,"/fig_time_series.png",sep=""), width = 3000, height = 3000)
plotTimeSeries(TS)
dev.off()

#
###


###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model (i.e. interpolate each variable in time series and compute interpolated temporal derivative to approximate temporal dynamics)

## parameters of observation model
N       <<- ncol(TS) - 1       # number of variables
K_o     <<- 3                  # number of ensemble elements
W_o     <<- rep(30, N)         # number of neurons in observation model (i.e. single layer perceptron)
N_o     <<- W_o * 3            # total number of parameters in observation model
rho     <<- 1                  # proportion of best samples to reject (to refine quality of fit if necessary)
alpha_i <<- 1                  # upsampling interpolation factor (2 doubles the number of points in interpolated time series)

## train observation model
runtime_o <<- system.time({
    model_o = trainModel_o(TS, alpha_i, N_o, K_o, rho)
})[3]
Yhat_o      <<- model_o$Yhat_o
ddt.Yhat_o  <<- model_o$ddt.Yhat_o
Omega_o     <<- model_o$Omega_o

## visualise observation model fit
pdf(paste(pathToOut,"/fig_predictions_o.pdf",sep=""))
# png(paste(pathToOut,"/fig_predictions_o.png",sep=""), width = 640, height = 640)
plotModel_o(TS, alpha_i, Yhat_o, ddt.Yhat_o)
dev.off()

## save results
save(Yhat_o,     file=paste(pathToOut,"/","Yhat_o.RData", sep=""))
save(ddt.Yhat_o, file=paste(pathToOut,"/","ddt.Yhat_o.RData", sep=""))
save(Omega_o,    file=paste(pathToOut,"/","Omega_o.RData", sep=""))
write(runtime_o, file=paste(pathToOut,"/","runtime_o.txt", sep=""))


### save results as csv
# preparation: calculating the mean
E.Yhat_o.all <- data.table(Time = TS$Time)
E.ddt.Yhat_o.all <- data.table(Time = TS$Time)

for(i in 1:n_taxa) {
  name_i = names(Yhat_o)[i]
  E.Yhat_o       = apply(Yhat_o[[i]],2,mean)
  E.ddt.Yhat_o   = apply(ddt.Yhat_o[[i]],2,mean)
  E.Yhat_o.all[, new := E.Yhat_o]
  colnames(E.Yhat_o.all) <- c(head(colnames(E.Yhat_o.all), -1), name_i)
  E.ddt.Yhat_o.all[, new := E.ddt.Yhat_o]
  colnames(E.ddt.Yhat_o.all) <- c(head(colnames(E.ddt.Yhat_o.all), -1), name_i)
}

# save as csv files
fwrite(
  E.Yhat_o.all,
  paste0(pathToOut, "/Yhat.csv")
)
fwrite(
  E.ddt.Yhat_o.all,
  paste0(pathToOut, "/Yhat_ddt.csv"),
  row.names = F
)
###


#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model (i.e. explain the per-capita growth rate of the populations calculated as 1/Y*dY/dt as a function of the states Y(t))

## load model_o
load(file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
load(file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
load(file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

## parameters of process model
K_p   <<- 3                      # number of models to fit
W_p   <<- rep(10, N)             # number of neurons in single layer perceptron (SLP)
N_p   <<- 2 * W_p * (2+N)        # number of parameters in process model
sd1_p <<- 0.1                    # standard deviation of model likelihood
sd2_p <<- as.list(rep(0.1,N))    # standard deviation of prior distributions

## train process model
runtime_p <<- system.time({
    model_p = trainModel_p(Yhat_o, ddt.Yhat_o, N_p, sd1_p, sd2_p, K_p, trainSplit=2/3)
})[3]
Yhat_p     <<- model_p$Yhat_p
ddx.Yhat_p <<- model_p$ddx.Yhat_p
Geber_p    <<- model_p$Geber_p
Omega_p    <<- model_p$Omega_p

## visualise process model
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
# png(paste(pathToOut,"/fig_predictions_p.png",sep=""), width = 640, height = 640)
plotModel_p(TS, alpha_i, Yhat_p, ddx.Yhat_p, Geber_p)
dev.off()

## store results
save(Yhat_p,     file=paste(pathToOut,"/","Yhat_p.RData", sep=""))
save(ddx.Yhat_p, file=paste(pathToOut,"/","ddx.Yhat_p.RData", sep=""))
save(Geber_p,    file=paste(pathToOut,"/","Geber_p.RData", sep=""))
save(Omega_p,    file=paste(pathToOut,"/","Omega_p.RData", sep=""))
write(runtime_p, file=paste(pathToOut,"/","runtime_p.txt", sep=""))

#
###


##############
## ANALYSIS ##
##############

## load results process model
load(paste(pathToOut,"/","Yhat_p.RData", sep=""))
load(paste(pathToOut,"/","ddx.Yhat_p.RData", sep=""))
load(paste(pathToOut,"/","Geber_p.RData", sep=""))
load(paste(pathToOut,"/","Omega_p.RData", sep=""))

## compute Jacobian and contribution matrix
MSq <<- function(x) mean(x^2)
prop <<- function(x) x/sum(x)
J <<- t(matrix(unlist(lapply(ddx.Yhat_p,function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,mean))),ncol=ncol(TS)-1)) ## average across samples then average across time steps
C = t(matrix(unlist(lapply(Geber_p,   function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,MSq))),ncol=ncol(TS)-1)) ## average across samples then take mean square across time steps
C = t(apply(C,1,prop))

## thresholding
J <<- J*(C>0.1)
C = C*(C>0.1)

## dynamical interaction plot (v1)
pdf(paste(pathToOut,"/fig_DIN.pdf",sep=""), width = 10, height = 10)
# png(paste(pathToOut,"/fig_DIN.png",sep=""), width = 640, height = 640)
plotDIN(J, C, colnames(TS)[-1])
dev.off()


## prepare plot of table J and C
names <-
  data.frame(var = paste0("x", 1:n_taxa), names = colnames(TS[1:n_taxa+1]))
colnames(J) <- names$var
colnames(C) <- names$var

# table J (effectsMat)
dat1 <-
  as.data.frame(J) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(
    Var1 = factor(paste0("x", Var1), levels = paste0("x", (1:n_taxa))),
    Var2 = factor(Var2, levels = paste0("x", rev(1:n_taxa)))
  )

# create and save plot
plt_J <- ggplot(dat1, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2() +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(axis.title=element_blank())
ggsave(paste(pathToOut,"/effectsMat.pdf",sep=""),
       plt_J,
       width = (n_taxa +1), height = n_taxa)

# table C (effectsMat)
dat2 <-
  as.data.frame(C) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(
    Var1 = factor(paste0("x", Var1), levels = paste0("x", (1:n_taxa))),
    Var2 = factor(Var2, levels = paste0("x", rev(1:n_taxa)))
  )

# create and save plot
plt_J <- ggplot(dat2, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(axis.title=element_blank())
ggsave(paste(pathToOut,"/weightsMat.pdf",sep=""),
       plt_J,
       width = (n_taxa +1), height = n_taxa)

# save J as csv files
colnames(J) = colnames(TS)[-1]
write.csv(
  J,
  paste0(pathToOut, "/effectsMat.csv"),
  row.names = F
)
write.csv(
  names,
  paste0(pathToOut, "/Names_xi.csv"),
  row.names = F
)

## save weightsMat C
colnames(C) = colnames(TS)[-1]
write.csv(
  C,
  paste0(pathToOut, "/weightsMat.csv"),
  row.names = F
)

###
}

#
###

# ####################################
# ## CROSS-VALIDATION PROCESS MODEL ##
# ####################################
#
# ## goal: perform cross-validation of process model by predicting the second half of the time series (of the interpolated p.c. growth rates)
# ##       the cross-validation is performed to identify the best value for regularisation on the nonlinear part of the model
#
# ## load model_o
# load(file=paste(pathToOut,"/","Yhat_o.RData", sep=""))
# load(file=paste(pathToOut,"/","ddt.Yhat_o.RData", sep=""))
# load(file=paste(pathToOut,"/","Omega_o.RData", sep=""))
#
# ## parameters for cross-validation
# K_p             = 1                            # number of models to fit per folds and regularisation parameter
# max_fold        = 2/3                          # beginning of test set
# folds           = list(c(0,1/2) * max_fold,
#                        c(1/2,1) * max_fold)    # proportion of the data in each fold
# crossValParVect = seq(0.005,0.05,0.005)        # vector of values of the regularisation parameter
#
# ## run cross-validation
# resultsCrossVal_p = crossVal_p(TS, alpha_i, Yhat_o, ddt.Yhat_o, N_p, sd1_p, sd2_p, K_p, folds, crossValParVect)
#
# ## visualise
# pdf(paste(pathToOut,"/fig_crossVal_p.pdf",sep=""))
# # png(paste(pathToOut,"/fig_crossVal_p.png",sep=""), width = 640, height = 640)
# plotCrossVal_p(resultsCrossVal_p)
# dev.off()
#
# ## store results
# save(resultsCrossVal_p,file=paste(pathToOut,"/","crossVal_p.RData"   ,sep=""))
#
# #
# ###

