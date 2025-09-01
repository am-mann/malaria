library(tidyverse)
library(rstan)
library(shinystan)
library(gridExtra)
library("bayesplot")
library("tidybayes")
library(matrixStats)
library("loo")
library(ncdf4)
library(here)
options(mc.cores = parallel::detectCores())


# Project root and output dir (use "2-param-fit" if this script is for the 2-param version)
root_dir   <- here::here()
fit_out_dir <- here("4-param-fit")
dir.create(fit_out_dir, recursive = TRUE, showWarnings = FALSE)

# Keep other outputs (CSV/PNG) alongside the RDS files
path <- fit_out_dir

beta_h = 0.022  #0.04; #
beta_v = 0.48  #0.09; #
mu_h = 1/(60.3*12)  # 0.000056 #
mu_v = 0.033*30   #1/(15/30)  15 days
sigma_h  = 0.1*30  # 1/(7/30) #14 days  # rate of transition from exposed to symptomatic infectious
gamma_h = 0.0035*30 # 1/(7/30)   15 days    # recovery rate

regions <- c("Centre", "Adamaoua", "ExtremeNord", "Nord", "Littoral", "Sud", "Ouest", "NordOuest", "SudOuest", "Est")

L_h_list <- c(Centre = 10327,  # Centre      (Central)
              Adamaoua = 3080, # Adamaoua 
              ExtremeNord= 11885, #ExtremeNord 
              Nord = 8806,   # Nord        (North)
              Littoral = 8352,   # Littoral    (Littoral)
              Sud = 1067, 
              Ouest= 1968,   # Ouest       (West)
              NordOuest= 2570,   # NordOuest   (North-West)
              SudOuest= 3294,   # SudOuest    (South-West)
              Est= 2171)   # Est         (East)

TotalPop_list <- c(
  Centre     = 4723372,  # Central (Centre)
  Adamaoua   = 1309666,  # Adamawa (Adamaoua)
  Est        = 1120567,  # East (Est)
  ExtremeNord= 4595669,  # Far North (Extreme-Nord)
  Littoral   = 3887698,  # Littoral
  Nord       = 2856872,  # North (Nord)
  NordOuest  = 2246302,  # Northwest (Nord-Ouest)
  Ouest      = 2089297,  # West (Ouest)
  Sud        = 805741,   # South (Sud)
  SudOuest   = 1857169   # Southwest (Sud Ouest)
)

nmonths = 12

# Organizing the cases data by region
#ncfname <- paste(data_path2, "cameroon_Adamaoua_monthly", ".nc", sep="")
#dname <- "t2m" 
#ncin <- nc_open(ncfname)
#temp1 <- ncvar_get(ncin,"t2m")
#convert <- rep(273.15, 36)
#temp <- temp1[50,50,]-rep(273.15, 36)

for (region in regions) {
  TotalPop <- TotalPop_list[[region]]
  
  L_h <- L_h_list[region]
  # Organizing the cases data by region
  data_path  <- here("RegionalData")
  Raw_Cases = read.csv(paste(data_path, "/Region ", region,".csv", sep="") )  
  
temp1 <- read_csv(paste(data_path, "/Cameroon_", region, "_daily.csv", sep=""), skip=10, col_names = c("YEAR", "MO", "DY", "T2M"))

# Combine YEAR, MO, and DY into a single date column
temp2 <- temp1 %>%
  mutate(Date = as.Date(paste(YEAR, MO, DY, sep = "-"), format = "%Y-%m-%d"))

# Calculate weekly averages
weekly_data <- temp2 %>%
  group_by(Week = floor_date(Date, "week")) %>%
  summarise(Weekly_Avg_Temp = mean(T2M, na.rm = TRUE))

# Calculate monthly averages
monthly_data <- temp2 %>%
  group_by(Month = floor_date(Date, "month")) %>%
  summarise(Monthly_Avg_Temp = mean(T2M, na.rm = TRUE))
temp <- monthly_data$Monthly_Avg_Temp
temp<- temp[1:nmonths]


# less than 5 years
Raw_Cases_L5 <- Raw_Cases %>%
	filter(Age_group == "Age_L5")

Raw_Cases_G5 <- Raw_Cases %>%
	filter(Age_group == "Age_G5")

ReportCases <- Raw_Cases_L5$Cases + Raw_Cases_G5$Cases
#Dates <- as.Date(
#  unique(
#    Raw_Cases$Date[format(as.Date(Raw_Cases$Date), "%Y") %in% c("2019", "2020", "2021")]
#  )
#)
Dates = as.Date(
  unique(
    Raw_Cases$Date[format(as.Date(Raw_Cases$Date), "%Y") == "2019"]
  ))
MonthlyCases <- ReportCases
MonthlyCases <- MonthlyCases[1:nmonths]

# initial condition
X0 = rep(0, 7)  # initializing the initial population vector
X0[1] = TotalPop  # initial susceptible population
X0[2] = 1 #(MonthlyCases[1]/sigma_h)/2   # Eh
X0[3] = 1 #(MonthlyCases[1]/sigma_h)/2   # Ih
X0[4] = 0        #Rh
X0[5] = 5e3  # Sv0 is assumed
X0[6] = 1    # Iv0 is fitted

initial_infected = sum(X0[2:3])

state = c(X0[2], X0[3])/sum(X0[2:3])
X0[2] = state[1]   # Eh
X0[3] = state[2]   # Ih

y0_vars = X0



# Preparing the parameters for stan
t0 = 0								# Initial time
# n_months = nrow(ReportCases)     # number of weeks
n_months = length(MonthlyCases)     # number of weeks
ts = 1:n_months

data_seir = list(
	L_h = L_h, 
	beta_h = beta_h,
	beta_v = beta_v,
	mu_h = mu_h,
	sigma_h  = sigma_h,
	gamma_h = gamma_h,
	mu_v = mu_v,
	PIE = pi,
	TotalPop = TotalPop,

	t0 = t0,
	n_months = n_months,
	ts = ts,
	inference = 1,
	doprint = 0,
	
	
	# Data to fit
	Reported_cases = MonthlyCases,
	temperature = temp,
	
	# Priors   Note: I have not used them in the computation yet, just created place holders
	p_i0 = c(log(5000), 1),   
	#p_t_min = c(11.7, 5), #
	p_q = c((0.1), 0.5),#
	p_Iv0 = c(log(500), 0.5),#
	#p_t_max = c(42.3, 5),  #  1/12
	p_wh = c(log(0.2), 0.5),  #
	p_delta = c(log(0.02), 0.5), #
	#p_tc = c(0.000828, 0.1), #
	p_phi = c(5), # parameters for the prior for phi 
	#p_A = c(0.0011, 0.001),  #from mordecai 2019 (for Anopheles gambiae!)
	#p_B = c(-0.049, 0.005),  
	#p_C = c(0.609, 0.5), 
	#p_mu_v_q = c(1, 0.000005),
	#p_q_MDR = c((0.0000952), 0.1),  #from mordecai 2019 (for Anopheles gambiae!)
	#p_t_min_MDR = c(6.35,10),  
	#p_t_max_MDR = c(25, 10), 
	#p_death = c(log(1),10),
	p_birth = c(5000,300),
	p_q_mu = c(2,5),
	#p_q_kappa = c(0.000111, 0.001),  #from mordecai, 2013
	#p_t_min_b = c(log(5),100), #from mordecai, 2019 (for Aedes aegypti)
	#p_t_max_b = c(log(35),100), #from mordecai, 2019 (for Aedes aegypti)
	
	# initial condition
	y0_vars = y0_vars
)

MalariaModel_Full <- stan_model(here("Adamaoua.stan"))

# fit_seir <- sampling(AgeModel_Full,data = data_seir, iter = 5, init = function() inti_cond(data_seir),  chains = 4)
# fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 100,  chains = 4)

 fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 1000,  init = function() list( i0 = 5000, q = 0.05, q_mu = 2, Iv0 = 500, w_h = 0.5, delta = 0.1, q_birth = 1), chains = 4)
#fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 100, chains = 1, 
#                     init = function() list(i0=100, t_min=10, Iv0=1000))
# fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 100, chains = 1, init = function() list(i0 = 100, t_min = 10, p = 0.5, q = 1, Iv0 = 1000, t_max = 35, w_h = 0.5, delta = 0.1, q_2 = 0.1, phi = 5))
# check_hmc_diagnostics(sampling(MalariaModel_Full, data = data_seir, iter = 100, chains = 1, init = function() list(i0=100, t_min=10, Iv0=1000)))
# fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 100, chains = 1)

 ################### saving stan object ###################
 path = here()
 saveRDS(fit_seir, file = file.path(fit_out_dir, paste0("StanOutput_", region, ".RDS")))
 fit_sir_negbin <- readRDS(file.path(fit_out_dir, paste0("StanOutput_", region, ".RDS")))

 
 # saving the sampling output
 # Extracting the posteriors  from the output of the sampling and saving it in a file
 Extract <- rstan::extract(fit_sir_negbin)
 Post  <- data.frame(i0 = Extract$i0, birth = Extract$birth, q_b = Extract$q, q_mu = Extract$q_mu,
                     Iv_0 = Extract$Iv0,  delta = Extract$delta,
                     w_h = Extract$w_h, phi = Extract$phi)
 write.table(Post, file= paste(path, "Par_PostDist_", region, ".csv", sep=""))
 
 # parameters for the prior for phi 
 Post <- read.table(file= paste(path, "Par_PostDist_", region, ".csv", sep="") )

 
 
 # checking inference
 pars = c( "i0", "q", "birth","Iv0", 'delta', 'w_h', 'phi', 'q_mu') # specifying parameters of interest
 print(fit_seir, pars = pars)
 # stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
 # traceplot(fit_seir, pars = pars)
 # pairs(fit_seir, pars = pars)
 
 densityplot <- stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
 densityplot
 # saving the projecttion plot in a .png file
 png(filename= paste(path, region, "_Densityplot.png", sep=""))
 plot(densityplot)
 dev.off()
 
 
 tplot <- traceplot(fit_seir, pars = pars)
 tplot
 # saving the projecttion plot in a .png file
 png(filename= paste(path, region, "_Traceplot.png", sep=""))
 plot(tplot)
 dev.off()
 
 
 # pairsplot <- pairs(fit_seir, pars = pars)
 # pairsplot
 # saving the projecttion plot in a .png file
 # png(filename= "Fit_Summary/Adamaoua_pairsplot.png")
 # plot(pairsplot)
 # dev.off()

 
 # All age groups
 smr_pred <- cbind(as.data.frame(summary(fit_seir, pars = "pred_cases", probs = c(0.05, 0.25, 0.5, 0.75, 0.95))$summary), ts, MonthlyCases,  days = as.Date(Dates))
 colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
 
 Reg1_CI <- ggplot(smr_pred, mapping = aes(x = days)) +
   geom_ribbon(aes(ymin = X5., ymax = X95.), fill ="steelblue4",  alpha = 0.35) +
   geom_ribbon(aes(ymin = X25., ymax = X75.), fill ="steelblue4",  alpha = 0.45) +
   geom_line(mapping = aes(x = days, y = X50.), col = "navy", size = 2.5) +
   geom_point(mapping = aes(y = MonthlyCases, x=days), size = 5.5) +
   labs(x = "Time (Months)", y = "Reported cases") +  theme_bw() + theme(text = element_text(size = 45)) +
   ggtitle(region) +  #theme(plot.title = element_text(hjust = 0.5)) +
   theme(plot.title=element_text(hjust=0.08, vjust=-0.95, margin=margin(t=50,b=-30), size=45)) + 
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", size=2, fill=NA), 
         axis.line = element_line(colour = "black")) +
   theme(axis.text.x=element_text(angle = 45, hjust = 1.0))
 Reg1_CI
 
 # saving the projecttion plot in a .png file
 png(filename= paste(path, region, "_MalariaFit_F.png", sep=""),  width = 700, height = 600)
 plot(Reg1_CI)
 dev.off()
}

for(region in regions){
    cat(region,"\n")
    File_Path = "/Users/amymann/Documents/Malaria Modeling/R/fits/jul23-fits/"
    region_data <-read.table(paste(File_Path, "Par_PostDist_", region, ".csv", sep=""))
    region_Mat = data.matrix(region_data)
    region_quant_L <- data.frame(y_pred_0.50 = colQuantiles(region_Mat, probs = 0.50),
                                 y_pred_0.05 = colQuantiles(region_Mat, probs = 0.05),
                                 y_pred_0.95 = colQuantiles(region_Mat, probs = 0.95)) 
    print(region_quant_L)
    cat("\n")
    
}

 