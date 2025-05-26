

# panel models for Lesk and Mankin (2026)

rm(list = ls())
library(data.table)
library(stargazer)
library(Hmisc)
library(lfe)
library(cowplot)
library(RColorBrewer)
library(pracma)
library(car)
library(dplyr)
library(performance)
setwd('C:/Users/corey/Documents/research/pexq/')

# load basin IDs for SE clustering and basin-level models
basins <- fread(file = 'data/basin_id_df.csv')


############################

#P gini effect on TWS, CPC:
d <- fread(file = "pexq_df.csv")
d$basin_id <- basins$basin_id

df <- d %>%
  filter(!water_year %in% c("2002","2017","2022"))

mP =  felm(data = df, formula = LWE ~ P|gridindex + water_year)
mT = felm(data = df, formula = LWE ~ T|gridindex + water_year)
mPT = felm(data = df, formula = LWE ~ P + T|gridindex + water_year)
mPTG = felm(data = df, formula = LWE ~ P + T + Gp_anoms|gridindex + water_year)
mPTGP = felm(data = df, formula = LWE ~ P + T + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year)
stargazer(mP,mT,mPT,mPTG,mPTGP, type = "text",out='TWS_basin_cluster.txt')

# basin clustering

dfc <- df %>%
  filter(!is.na(basin_id))

mPTG = felm(data = dfc, formula = LWE ~ P + T + Gp_anoms|gridindex + water_year)
mPTGP = felm(data = dfc, formula = LWE ~ P + T + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year)
mPTGc = felm(data = dfc, formula = LWE ~ P + T + Gp_anoms|gridindex + water_year|0|basin_id)
mPTGPc = felm(data = dfc, formula = LWE ~ P + T + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year|0|basin_id)
stargazer(mPTG,mPTGP,mPTGc,mPTGPc, type = "text",out='TWS_basin_cluster.txt')#,out='table1.txt')

# covariance among coefficients
vcov_matrix <- vcov(mPTGPc)
vcov_matrix

#save
coefs <- summary(mPTGPc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P", "Estimate"],
  se_P = coefs["P", "Cluster s.e."],
  pval_P = coefs["P", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_anoms", "Estimate"],
  se_Gp = coefs["Gp_anoms", "Cluster s.e."],
  pval_Gp = coefs["Gp_anoms", "Pr(>|t|)"],
  
  coef_GpP = coefs["Gp_anoms:Pmean", "Estimate"],
  se_GpP = coefs["Gp_anoms:Pmean", "Cluster s.e."],
  pval_GpP = coefs["Gp_anoms:Pmean", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[3,4]
)

write.csv(model_info, "CPC_regression_results.csv")

####################
# basin-scale models :
# here, need to run separately for different data products
# to do so, simply change the dataset suffix when loading in the panel file 

d <- fread(file = "pexq_df_GPCP.csv")
#d$basin_id <- basins$basin_id

df <- d %>%
  filter(!water_year %in% c("2002","2017","2022"))

#get unique non-NA basin ids
basin_ids <- unique(basins$basin_id)
basin_ids <- basin_ids[-1]

#holder df w interaction:
model_results_interaction <- data.frame(
  basin_id = integer(),
  coef_P = numeric(),
  se_P = numeric(),
  pval_P = numeric(),
  
  coef_T = numeric(),
  se_T = numeric(),
  pval_T = numeric(),
  
  coef_Gp = numeric(),
  se_Gp = numeric(),
  pval_Gp = numeric(),
  
  coef_Gp_Pmean = numeric(),
  se_Gp_Pmean = numeric(),
  pval_Gp_Pmean = numeric(),
  
  R2 = numeric(),
  model_pval = numeric(),
  stringsAsFactors = FALSE
)

#holder df w/o interaction:
model_results <- data.frame(
  basin_id = integer(),
  coef_P = numeric(),
  se_P = numeric(),
  pval_P = numeric(),
  
  coef_T = numeric(),
  se_T = numeric(),
  pval_T = numeric(),
  
  coef_Gp = numeric(),
  se_Gp = numeric(),
  pval_Gp = numeric(),
  
  R2 = numeric(),
  model_pval = numeric(),
  stringsAsFactors = FALSE
)

#loop basins
for (b in basin_ids){
print(b)
#subset data to basin  
dfc <- df %>%
  filter(basin_id==b)

if (nrow(dfc) == 0) {
  print("No data for basin_id:")
  next
}

if (length(unique(dfc$gridindex))==1) {
  print("1 pt basin:")
  next
}
#fit model
basinPTGPc = felm(data = dfc, formula = LWE ~ P_GPCP + T + Gp_GPCP|gridindex + water_year)

# Extract coefficients, standard errors, p-values, R2, and model p-value
coefs <- summary(basinPTGPc)$coefficients

model_info <- data.frame(
  basin_id = b,
  coef_P = coefs["P_GPCP", "Estimate"],
  se_P = coefs["P_GPCP", "Std. Error"],
  pval_P = coefs["P_GPCP", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Std. Error"],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_GPCP", "Estimate"],
  se_Gp = coefs["Gp_GPCP", "Std. Error"],
  pval_Gp = coefs["Gp_GPCP", "Pr(>|t|)"],
  
  R2 = summary(basinPTGPc)$r2adj,
  model_pval = as.numeric(summary(basinPTGPc)$P.fstat['p.F'])
)

# Append to results df
model_results <- rbind(model_results, model_info)

}

write.csv(model_results, "basin_models_no_interaction_GPCP.csv")

###########################
#Other obs precip products:

###
#GPCP precip gini
d = fread(file = "pexq_df_GPCP.csv")

df <- d %>%
  filter(!water_year %in% c("2002","2017","2022"))

mP =  felm(data = df, formula = LWE ~ P_GPCP|gridindex + water_year)
mT = felm(data = df, formula = LWE ~ T|gridindex + water_year)
mPT = felm(data = df, formula = LWE ~ P_GPCP + T|gridindex + water_year)
mPTG = felm(data = df, formula = LWE ~ P_GPCP + T + Gp_GPCP|gridindex + water_year)
mPTGP = felm(data = df, formula = LWE ~ P_GPCP + T + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year)
stargazer(mP,mT,mPT,mPTG,mPTGP, type = "text")#,out='table1.txt')

#clustering
dfc <- df %>%
  filter(!is.na(basin_id))

mPTGc = felm(data = dfc, formula = LWE ~ P_GPCP + T + Gp_GPCP|gridindex + water_year|0|basin_id)
mPTGPc = felm(data = dfc, formula = LWE ~ P_GPCP + T + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year|0|basin_id)
stargazer(mPTGc,mPTGPc, type = "text")#,out='table1.txt')

# covariance among coefficients
vcov_matrix <- vcov(mPTGPc)

coefs <- summary(mPTGPc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P_GPCP", "Estimate"],
  se_P = coefs["P_GPCP", "Cluster s.e."],
  pval_P = coefs["P_GPCP", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_GPCP", "Estimate"],
  se_Gp = coefs["Gp_GPCP", "Cluster s.e."],
  pval_Gp = coefs["Gp_GPCP", "Pr(>|t|)"],
  
  coef_GpP = coefs["Gp_GPCP:Pmean_GPCP", "Estimate"],
  se_GpP = coefs["Gp_GPCP:Pmean_GPCP", "Cluster s.e."],
  pval_GpP = coefs["Gp_GPCP:Pmean_GPCP", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[3,4]
)

write.csv(model_info, "GPCP_regression_results.csv")


###
#GPCC precip gini
d = fread(file = "pexq_df_GPCC.csv")

df <- d %>%
  filter(!water_year %in% c("2002","2017","2021"))

mP =  felm(data = df, formula = LWE ~ P_GPCC|gridindex_x + water_year_x)
mT = felm(data = df, formula = LWE ~ T|gridindex_x + water_year_x)
mPT = felm(data = df, formula = LWE ~ P_GPCC + T|gridindex_x + water_year_x)
mPTG = felm(data = df, formula = LWE ~ P_GPCC + T + Gp_GPCC|gridindex_x + water_year_x)
mPTGP = felm(data = df, formula = LWE ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex_x + water_year_x)
stargazer(mP,mT,mPT,mPTG,mPTGP, type = "text")#,out='table1.txt')

#clustering
dfc <- df %>%
  filter(!is.na(basin_id))

mPTGc = felm(data = dfc, formula = LWE ~ P_GPCC + T + Gp_GPCC|gridindex + water_year|0|basin_id)
mPTGPc = felm(data = dfc, formula = LWE ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year|0|basin_id)
stargazer(mPTGc,mPTGPc, type = "text")#,out='table1.txt')


# covariance among coefficients
vcov_matrix <- vcov(mPTGPc)
vcov_matrix

coefs <- summary(mPTGPc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P_GPCC", "Estimate"],
  se_P = coefs["P_GPCC", "Cluster s.e."],
  pval_P = coefs["P_GPCC", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_GPCC", "Estimate"],
  se_Gp = coefs["Gp_GPCC", "Cluster s.e."],
  pval_Gp = coefs["Gp_GPCC", "Pr(>|t|)"],
  
  coef_GpP = coefs["Gp_GPCC:Pmean_GPCC", "Estimate"],
  se_GpP = coefs["Gp_GPCC:Pmean_GPCC", "Cluster s.e."],
  pval_GpP = coefs["Gp_GPCC:Pmean_GPCC", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[3,4]
)

write.csv(model_info, "GPCC_regression_results.csv")

#################
# Fig. 3: link between Gp and downwelling surface shortwave:

d = fread(file = "data/pexq_df_CPC_SW.csv")

sPTGP = felm(data = d, formula = SW ~ P + T + Gp + Gp:Pmean|gridindex + water_year|0|basin_id)
stargazer(sPTGP, type = "text")#,out='table1.txt')
vcov_CPC = vcov(sPTGP)[3,4]
coefs_CPC <- summary(sPTGP)$coefficients

d = fread(file = "data/pexq_df_GPCP_SW.csv")

sPTGP = felm(data = d, formula = SW ~ P_GPCP + T + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year|0|basin_id)
stargazer(sPTGP, type = "text")#,out='table1.txt')
vcov_GPCP = vcov(sPTGP)[3,4]
coefs_GPCP <- summary(sPTGP)$coefficients

d = fread(file = "data/pexq_df_GPCC_SW.csv")

sPTGP = felm(data = d, formula = SW ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year|0|basin_id)
stargazer(sPTGP, type = "text")#,out='table1.txt')
vcov_GPCC = vcov(sPTGP)[3,4]
coefs_GPCC <- summary(sPTGP)$coefficients


model_info <- data.frame(
  
  coef_Gp_CPC = coefs_CPC["Gp", "Estimate"],
  se_Gp_CPC = coefs_CPC["Gp", "Cluster s.e."],
  pval_Gp_CPC = coefs_CPC["Gp", "Pr(>|t|)"],
  coef_GpP_CPC = coefs_CPC["Gp:Pmean", "Estimate"],
  se_GpP_CPC = coefs_CPC["Gp:Pmean", "Cluster s.e."],
  pval_GpP_CPC = coefs_CPC["Gp:Pmean", "Pr(>|t|)"],
  cov_CPC = vcov_CPC,
  
  coef_Gp_GPCP = coefs_GPCP["Gp_GPCP", "Estimate"],
  se_Gp_GPCP = coefs_GPCP["Gp_GPCP", "Cluster s.e."],
  pval_Gp_GPCP = coefs_GPCP["Gp_GPCP", "Pr(>|t|)"],
  coef_GpP_GPCP = coefs_GPCP["Gp_GPCP:Pmean_GPCP", "Estimate"],
  se_GpP_GPCP = coefs_GPCP["Gp_GPCP:Pmean_GPCP", "Cluster s.e."],
  pval_GpP_GPCP = coefs_GPCP["Gp_GPCP:Pmean_GPCP", "Pr(>|t|)"],
  cov_GPCP = vcov_GPCP,
  
  coef_Gp_GPCC = coefs_GPCC["Gp_GPCC", "Estimate"],
  se_Gp_GPCC = coefs_GPCC["Gp_GPCC", "Cluster s.e."],
  pval_Gp_GPCC = coefs_GPCC["Gp_GPCC", "Pr(>|t|)"],
  coef_GpP_GPCC = coefs_GPCC["Gp_GPCC:Pmean_GPCC", "Estimate"],
  se_GpP_GPCC = coefs_GPCC["Gp_GPCC:Pmean_GPCC", "Cluster s.e."],
  pval_GpP_GPCC = coefs_GPCC["Gp_GPCC:Pmean_GPCC", "Pr(>|t|)"],
  cov_GPCC = vcov_GPCC
)

write.csv(model_info, "SW_Gp_regression_results.csv")

###############
# Fig. 3: Partioning distribution versus shortwave effect (control for SW in main model)
# CPC: TWS = f(Gp,SW)
d = fread(file = "data/pexq_df_CPC_SW.csv")

df <- d %>%
  filter(!water_year %in% c("2002","2017","2021"))

mPTGPc = felm(data = df, formula = LWE ~ P + T +SW + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year|0|basin_id)
mPTGPSc = felm(data = df, formula = LWE ~ P + T + SW +SW:Pmean + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year|0|basin_id)

stargazer(mPTGPc,mPTGPSc,type='text')

vcov_matrix <- vcov(mPTGPSc)
vcov_matrix[4,6]
coefs <- summary(mPTGPSc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P", "Estimate"],
  se_P = coefs["P", "Cluster s.e."],
  pval_P = coefs["P", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_anoms", "Estimate"],
  se_Gp = coefs["Gp_anoms", "Cluster s.e."],
  pval_Gp = coefs["Gp_anoms", "Pr(>|t|)"],
  
  coef_GpP = coefs["Pmean:Gp_anoms", "Estimate"],
  se_GpP = coefs["Pmean:Gp_anoms", "Cluster s.e."],
  pval_GpP = coefs["Pmean:Gp_anoms", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[4,6]
)

write.csv(model_info, "CPC_regression_results_SWint.csv")

# GPCP: TWS = f(Gp,SW)
d = fread(file = "pexq_df_GPCP_SW.csv")
df <- d %>%
  filter(!water_year %in% c("2002","2017","2021"))
dfc <- df %>%
  filter(!is.na(basin_id))

mPTGPc = felm(data = dfc, formula = LWE ~ P_GPCP + T + SW + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year|0|basin_id)
mPTGPSc = felm(data = dfc, formula = LWE ~ P_GPCP + T + SW + SW:Pmean_GPCP + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year|0|basin_id)

stargazer(mPTGPc,mPTGPSc, type = "text")#,out='table1.txt')

# covariance among coefficients
vcov_matrix <- vcov(mPTGPSc)

coefs <- summary(mPTGPSc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P_GPCP", "Estimate"],
  se_P = coefs["P_GPCP", "Cluster s.e."],
  pval_P = coefs["P_GPCP", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_GPCP", "Estimate"],
  se_Gp = coefs["Gp_GPCP", "Cluster s.e."],
  pval_Gp = coefs["Gp_GPCP", "Pr(>|t|)"],
  
  coef_GpP = coefs["Pmean_GPCP:Gp_GPCP", "Estimate"],
  se_GpP = coefs["Pmean_GPCP:Gp_GPCP", "Cluster s.e."],
  pval_GpP = coefs["Pmean_GPCP:Gp_GPCP", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[4,6]
)

write.csv(model_info, "GPCP_regression_results_SWint.csv")

# GPCC: TWS = f(Gp,SW)
d = fread(file = "pexq_df_GPCC_SW.csv")
df <- d %>%
  filter(!water_year %in% c("2002","2017","2021"))

#clustering
dfc <- df %>%
  filter(!is.na(basin_id))


mPTGPc = felm(data = dfc, formula = LWE ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year|0|basin_id)
mPTGPSc = felm(data = dfc, formula = LWE ~ P_GPCC +T + SW + SW:Pmean_GPCC + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year|0|basin_id)

stargazer(mPTGPc,mPTGPSc, type = "text")#,out='table1.txt')


# covariance among coefficients
vcov_matrix <- vcov(mPTGPSc)
vcov_matrix

#setup data:
coefs <- summary(mPTGPSc)$coefficients

model_info <- data.frame(
  coef_P = coefs["P_GPCC", "Estimate"],
  se_P = coefs["P_GPCC", "Cluster s.e."],
  pval_P = coefs["P_GPCC", "Pr(>|t|)"],
  
  coef_T = coefs["T", "Estimate"],
  se_T = coefs["T", "Cluster s.e."],
  pval_T = coefs["T", "Pr(>|t|)"],
  
  coef_Gp = coefs["Gp_GPCC", "Estimate"],
  se_Gp = coefs["Gp_GPCC", "Cluster s.e."],
  pval_Gp = coefs["Gp_GPCC", "Pr(>|t|)"],
  
  coef_GpP = coefs["Pmean_GPCC:Gp_GPCC", "Estimate"],
  se_GpP = coefs["Pmean_GPCC:Gp_GPCC", "Cluster s.e."],
  pval_GpP = coefs["Pmean_GPCC:Gp_GPCC", "Pr(>|t|)"],
  
  R2 = summary(mPTGPc)$r2adj,
  model_pval = as.numeric(summary(mPTGPc)$P.fstat['p.F']),
  vcov = vcov_matrix[4,6]
)

write.csv(model_info, "GPCC_regression_results_SWint.csv")


################
# Fig 3. Gp effect on ET
dfc <- d %>%
  filter(Pmean_CPC > 0 & Pmean_CPC < 1000)

ePTGP = felm(data = d, formula = ET_GLEAM ~ P_CPC + T + Gp_CPC +Gp_CPC:Pmean_CPC |gridindex + water_year)
stargazer(ePTGP, type = "text")

# ET basin clustering
d = fread(file = "pexq_df_ET_allPproducts.csv")

df <- d %>%
  filter(!water_year %in% c("2002","2017","2022"))

dfc <- df %>%
  filter(!is.na(basin_id))

ePTGP_CPC = felm(data = dfc, formula = ET_GLEAM ~ P + T + Gp_anoms + Gp_anoms:Pmean|gridindex + water_year|0|basin_id)
ePTGP_GPCC = felm(data = dfc, formula = ET_GLEAM ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year|0|basin_id)
ePTGP_GPCP = felm(data = dfc, formula = ET_GLEAM ~ P_GPCP + T + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year|0|basin_id)

stargazer(ePTGP_CPC,ePTGP_GPCC,ePTGP_GPCP, type = "text")#,out='table1.txt')
vcov_CPC <- vcov(ePTGP_CPC)[3,4]
vcov_GPCC <- vcov(ePTGP_GPCC)[3,4]
vcov_GPCP <- vcov(ePTGP_GPCP)[3,4]

coefs <- summary(ePTGP_CPC)$coefficients
coefs_GPCP <- summary(ePTGP_GPCP)$coefficients
coefs_GPCC <- summary(ePTGP_GPCC)$coefficients

model_info <- data.frame(
  
  coef_Gp_CPC = coefs["Gp_anoms", "Estimate"],
  se_Gp_CPC = coefs["Gp_anoms", "Cluster s.e."],
  coef_GpP_CPC = coefs["Gp_anoms:Pmean", "Estimate"],
  se_GpP_CPC = coefs["Gp_anoms:Pmean", "Cluster s.e."],
  vcov_CPC = vcov_CPC,
  
  coef_Gp_GPCP = coefs_GPCP["Gp_GPCP", "Estimate"],
  se_Gp_GPCP = coefs_GPCP["Gp_GPCP", "Cluster s.e."],
  coef_GpP_GPCP = coefs_GPCP["Gp_GPCP:Pmean_GPCP", "Estimate"],
  se_GpP_GPCP = coefs_GPCP["Gp_GPCP:Pmean_GPCP", "Cluster s.e."],
  vcov_GPCP = vcov_GPCP,
  
  coef_Gp_GPCC = coefs_GPCC["Gp_GPCC", "Estimate"],
  se_Gp_GPCC = coefs_GPCC["Gp_GPCC", "Cluster s.e."],
  coef_GpP_GPCC = coefs_GPCC["Gp_GPCC:Pmean_GPCC", "Estimate"],
  se_GpP_GPCC = coefs_GPCC["Gp_GPCC:Pmean_GPCC", "Cluster s.e."],
  vcov_GPCC = vcov_GPCC
)

# ET longer time period
d = fread(file = "ET_df_1981-2022.csv")

ePTGP_CPC = felm(data = df, formula = ET_GLEAM ~ P_CPC + T + Gp_CPC + Gp_CPC:Pmean_CPC|gridindex + water_year)#|0|basin_id)
ePTGP_GPCC = felm(data = d, formula = ET_GLEAM ~ P_GPCC + T + Gp_GPCC + Gp_GPCC:Pmean_GPCC|gridindex + water_year)#|0|basin_id)
ePTGP_GPCP = felm(data = d, formula = ET_GLEAM ~ P_GPCP + T + Gp_GPCP + Gp_GPCP:Pmean_GPCP|gridindex + water_year)#|0|basin_id)

vcov_CPC <- vcov(ePTGP_CPC)

stargazer(ePTGP_CPC,type='text')
write.csv(model_info, "ET_regression_results.csv")


##################
# Fig. 1: hydroclimate drivers of Gp?

d = fread(file = "data/panels/metrics_std_df.csv")
mPTGP95 = felm(data = d, formula = Gp ~ P + T + WD95 + DD|gridindex + water_year)
mPTGP99_CPC = felm(data = d, formula = Gp ~ P + T + WD99 + DD|gridindex + water_year)

stargazer(mPTGP99_CPC, type = "text")#,out='table1.txt')
coefs_CPC <- summary(mPTGP99_CPC)$coefficients

d = fread(file = "data/panels/metrics_std_df_GPCC.csv")
mPTGP95 = felm(data = d, formula = Gp_GPCC ~ P_GPCC + T + WD95 + DD|gridindex + water_year|0|basin_id)
mPTGP99_GPCC = felm(data = d, formula = Gp_GPCC ~ P_GPCC + T + WD99 + DD|gridindex + water_year|0|basin_id)

stargazer(mPTGP99_GPCC, type = "text")#,out='table1.txt')
coefs_GPCC <- summary(mPTGP99_GPCC)$coefficients

d = fread(file = "data/panels/metrics_std_df_GPCP.csv")
mPTGP95 = felm(data = d, formula = Gp_GPCP ~ P_GPCP + T + WD95 + DD|gridindex + water_year|0|basin_id)
mPTGP99_GPCP = felm(data = d, formula = Gp_GPCP ~ P_GPCP + T + WD99 + DD|gridindex + water_year|0|basin_id)

stargazer(mPTGP99_GPCP, type = "text")#,out='table1.txt')
coefs_GPCP <- summary(mPTGP99_GPCP)$coefficients

model_info <- data.frame(
  coef_P_CPC = coefs_CPC["P", "Estimate"],
  se_P_CPC = coefs_CPC["P", "Std. Error"],
  coef_T_CPC = coefs_CPC["T", "Estimate"],
  se_T_CPC = coefs_CPC["T", "Std. Error"],
  coef_DD_CPC = coefs_CPC["DD", "Estimate"],
  se_DD_CPC = coefs_CPC["DD", "Std. Error"],
  coef_WD_CPC = coefs_CPC["WD99", "Estimate"],
  se_WD_CPC = coefs_CPC["WD99", "Std. Error"],
  
  coef_P_GPCC = coefs_GPCC["P_GPCC", "Estimate"],
  se_P_GPCC = coefs_GPCC["P_GPCC", "Cluster s.e."],
  coef_T_GPCC = coefs_GPCC["T", "Estimate"],
  se_T_GPCC = coefs_GPCC["T", "Cluster s.e."],
  coef_DD_GPCC = coefs_GPCC["DD", "Estimate"],
  se_DD_GPCC = coefs_GPCC["DD", "Cluster s.e."],
  coef_WD_GPCC = coefs_GPCC["WD99", "Estimate"],
  se_WD_GPCC = coefs_GPCC["WD99", "Cluster s.e."],
  
  coef_P_GPCP = coefs_GPCP["P_GPCP", "Estimate"],
  se_P_GPCP = coefs_GPCP["P_GPCP", "Cluster s.e."],
  coef_T_GPCP = coefs_GPCP["T", "Estimate"],
  se_T_GPCP = coefs_GPCP["T", "Cluster s.e."],
  coef_DD_GPCP = coefs_GPCP["DD", "Estimate"],
  se_DD_GPCP = coefs_GPCP["DD", "Cluster s.e."],
  coef_WD_GPCP = coefs_GPCP["WD99", "Estimate"],
  se_WD_GPCP = coefs_GPCP["WD99", "Cluster s.e."]
)

write.csv(model_info, "regression_results_drivers_of_Gp.csv")

# note that whether you use 95th or 99th percentile for 'wet days' doesn't matter much

# correlation between P and WD
correlation <- cor(d$WD95, d$WD99, use = "complete.obs")
# Print the result
print(correlation)
