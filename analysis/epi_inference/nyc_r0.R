# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

########################################################################
####################################EpiNow2####################################
########################################################################

library(EpiNow2)

# set current working directory
setwd("/Users/jonathanpekar/Desktop/Wertheim/mpox_nyc/nyc_mpox_phylogeography/analysis/epi_inference/") # set to your working directory
mpox_data <- read.csv(file = 'case_counts.forEpiNow.csv')
reported_cases_mpox <- data.table::data.table()
reported_cases_mpox$date <- seq(as.Date("2022-05-04"), as.Date("2023-10-07"), by="days")
reported_cases_mpox$confirm <- mpox_data$CaseNum
# subset to start at row 40
reported_cases_mpox <- reported_cases_mpox[29:196,]
generation_time <- Gamma(mean = 12.6, sd = 5.7, max= 25)
incubation_period <- Gamma(mean = 9.1, sd = 3, max = 20)
fit_epinow <- estimate_infections(data = reported_cases_mpox,
                                    generation_time = generation_time_opts(generation_time),
                                    delays = delay_opts(incubation_period),
                                    rt = rt_opts(prior = list(mean = 2, sd = 0.5)),
                                    horizon = 0,
                                    stan = stan_opts(cores = 4, 
                                          chains=4, 
                                          samples = 6400, # iterations = # samples / # chains + warmup
                                          warmup = 400,
                                          control = list(adapt_delta = 0.999)),
                                    CrIs = c(0.2, 0.5, 0.95))

# examine model fit
knitr::kable(summary(fit_epinow))
plot(fit_epinow)
plot(fit_epinow, "growth_rate")


# save the model
save(fit_epinow, file = "fit_epinow.RData")
# save estimates_mpox without index
write.csv(fit_epinow$summarised, file = "estimates_epinow.csv", row.names = FALSE)


########################################################################################################
# Eales et al., Epidemics (2022): Appropriately smoothing prevalence data to inform estimates of growth rate and reproduction number
# https://github.com/mrc-ide/reactidd/blob/master/vignettes/TemporalMethodsPaper/PHE_rounds1-7_analysis.rmd
########################################################################################################

knitr::opts_chunk$set(echo = TRUE)
library('reactidd')

# a few functions weren't installed for some reason, so I'm sourcing them manually. Install these from the github repo above. 
source('/Users/jonathanpekar/Desktop/Wertheim/mpox_nyc/other_studies/reactidd/R/plot_p_spline_igr.R')
source('/Users/jonathanpekar/Desktop/Wertheim/mpox_nyc/other_studies/reactidd/R/plot_p_spline_R.R')

# load data
mpox_data <- mpox_data[,1:2]
# reverse the order of the data of mpox_data so that it is in descending order
mpox_data <- mpox_data[rev(rownames(mpox_data)),]

colnames(mpox_data) <- c("date", "n_cases")
mpox_data$date <- as.Date(mpox_data$date, format = "%Y-%m-%d") # Adjust the format as per your date format

# define dates
# min_date_r1 <- as.Date("2022-06-01") 
min_date_r1 <- as.Date("2022-05-15") # needed to rerun the model with a different start date for plotting purposes
max_date_r7 <- as.Date("2022-11-15")

# fitting bayesian p-spline model 
fit_reactidd <- stan_p_spline_phe(X = mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$date,
                                   Y= mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$n_cases,
                                   target_dist_between_knots = 5,
                                   spline_degree = 3,
                                   iter = 5000,
                                   warmup = 1000,
                                   cores = 4,
                                   chains = 4,
                                   control = list(adapt_delta=0.999, max_treedepth = 20))


# plot the model fit with 95%CI and 50%CI
p_spline_plot <- plot_p_spline_phe(X = mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$date,
                                   Y= mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$n_cases,
                                   p_spline_fit = fit_reactidd, 
                                    target_dist_between_knots = 5,
                                  spline_degree = 3,
                                   ylim = 100.0)

print(p_spline_plot[[1]])


# We can also calculate and plot the instanteous growth rate over the study period
p_spline_igr <- plot_p_spline_igr(X = mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$date,
                                   p_spline_fit = fit_reactidd, 
                                    target_dist_between_knots = 5,
                                    spline_degree = 3,
                                  ylim = 0.15,
                                  link_function = "log")

print(p_spline_igr[[1]])


# and also the rolling two week average Reproduction number
p_spline_Rt <- plot_p_spline_R(X = mpox_data[mpox_data$date>=min_date_r1 & mpox_data$date<= max_date_r7,]$date,
                                   p_spline_fit = fit_reactidd, 
                                    target_dist_between_knots = 5,
                                    spline_degree = 3,
                               link_function = "log")

print(p_spline_Rt[[1]])

# add a column specifying that p_spline_igr is growth_rate, that is the length of p_spline_igr[[2]]$x, and call that column "variable"
p_spline_igr_results <- data.frame(p_spline_igr[[2]], variable = rep("growth_rate", length(p_spline_igr[[2]]$x)))
p_spline_Rt_results <- data.frame(p_spline_Rt[[2]], variable = rep("Rt", length(p_spline_Rt[[2]]$x)))
# combine rows of p_spline_igr and p_spline_Rt
combined_df <- rbind(p_spline_igr_results, p_spline_Rt_results)

# save model
save(fit_reactidd, file = "fit_reactidd.RData")
# save results
write.csv(combined_df, file = "estimates_reactidd.csv", row.names = FALSE)