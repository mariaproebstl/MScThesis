Simulation of the Van der Pol oscillator
================
Compiled at 2023-10-06 13:24:36 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "ed3de31b-cc20-4300-90cc-e98faa8c0c62")
```

The purpose of this document is to simulate the in “SPARSE
RECONSTRUCTION OF ORDINARY DIFFERENTIAL EQUATIONS WITH INFERENCE” by
Sara Venkatraman et al. described Van der Pol ODE.

``` r
library("conflicted")
library(data.table)
library(dplyr)
library(ggplot2)
library(phyloseq)

library(deSolve)
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Simulate Van der Pol ODE

### RK4 Solver

``` r
# Define the RK4 solver function
rk4 <- function(f, g, x1_0, x2_0, t0, tf, h) {
  n <- round((tf - t0) / h)  # Number of time steps
  t <- seq(t0, tf, by = h)  # Time points
  x1 <- numeric(n + 1)
  x2 <- numeric(n + 1)
  
  x1[1] <- x1_0
  x2[1] <- x2_0
  
  for (i in 1:n) {
    k1_x1 <- h * f(x1[i], x2[i])
    k1_x2 <- h * g(x1[i], x2[i])
    
    k2_x1 <- h * f(x1[i] + 0.5 * k1_x1, x2[i] + 0.5 * k1_x2)
    k2_x2 <- h * g(x1[i] + 0.5 * k1_x1, x2[i] + 0.5 * k1_x2)
    
    k3_x1 <- h * f(x1[i] + 0.5 * k2_x1, x2[i] + 0.5 * k2_x2)
    k3_x2 <- h * g(x1[i] + 0.5 * k2_x1, x2[i] + 0.5 * k2_x2)
    
    k4_x1 <- h * f(x1[i] + k3_x1, x2[i] + k3_x2)
    k4_x2 <- h * g(x1[i] + k3_x1, x2[i] + k3_x2)
    
    x1[i + 1] <- x1[i] + (k1_x1 + 2 * k2_x1 + 2 * k3_x1 + k4_x1) / 6
    x2[i + 1] <- x2[i] + (k1_x2 + 2 * k2_x2 + 2 * k3_x2 + k4_x2) / 6
  }
  
  return(data.frame(t, x1, x2))
}
```

### Van der Pol ODE system

``` r
# Define the system of differential equations
f <- function(x1, x2) {
  return(x2)
}

g <- function(x1, x2) {
  return(x1 + mu * (1 - x1^2) * x2)
}

# Set initial conditions and time span
mu = 2
x1_0 <- 1.0
x2_0 <- 0.0
t_start <- 0
t_end <- 15
h <- 0.05  # Step size
```

### Solve the ODE

``` r
# Solve the system using RK4
solution <- rk4(f, g, x1_0, x2_0, t_start, t_end, h)
```

### Plot the time series

``` r
# Create a data frame for plotting
data <- data.frame(t = solution$t, x1 = solution$x1, x2 = solution$x2)

# Plot the time series
ggplot(data, aes(x = t)) +
  geom_line(aes(y = x1, color = "x1(t)")) +
  geom_line(aes(y = x2, color = "x2(t)")) +
  labs(x = "Time", y = "Value", color = "Variable") +
  scale_color_manual(values = c("x1(t)" = "blue", "x2(t)" = "red")) +
  theme_minimal()
```

![](01g-timeseries-simulation-VdP_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Directly save simulated ODE

``` r
# save time series as csv file
write.csv(
  solution,
  path_target(paste0("ts_VanderPol.csv")),
  row.names = F
)
```

## Files written

These files have been written to the target directory,
`data/01g-timeseries-simulation-VdP`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 1 × 4
    ##   path             type         size modification_time  
    ##   <fs::path>       <fct> <fs::bytes> <dttm>             
    ## 1 ts_VanderPol.csv file        11.9K 2023-10-06 13:24:41
