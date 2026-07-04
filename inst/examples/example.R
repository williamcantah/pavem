library(readxl)
library(pavem)

data_file <- "africal.xlsx"

raw <- read_excel(data_file, sheet = "panel_data")

ecowas <- raw[raw$country_code %in% c("GHA", "NGA", "SLE", "SEN", "BEN"), ]

vars <- c("inflation", "gdp_growth", "interest_rate")

ecowas <- prep(ecowas, id = "country", time = "year", vars = vars)

csd_results <- csdt(ecowas, id = "country", time = "year", vars = vars)
unit_results <- urts(ecowas, id = "country", time = "year", vars = vars)
coin_results <- coint(ecowas, id = "country", time = "year",
                      y = "inflation", x = "interest_rate", method = "ped")

pvar_fit <- pvar(ecowas, id = "country", time = "year", vars = vars)
pvar_diag <- pvdi(pvar_fit)
pvar_irf <- pvir(pvar_fit, n.ahead = 10, draws = 200)
pvar_fevd <- pvfd(pvar_fit, n.ahead = 10)
pipl(pvar_irf, file = "pvar_irf.pdf")
pfpl(pvar_fevd, file = "pvar_fevd.pdf")

vecm_fit <- vecm(ecowas, id = "country", time = "year",
                 vars = vars, lags = 2, rank = 1)
vecm_diag <- vedi(vecm_fit)
vecm_irf <- veir(vecm_fit, n.ahead = 10, draws = 200)
vecm_fevd <- vefd(vecm_fit, n.ahead = 10)
vipl(vecm_irf, file = "vecm_irf.pdf")
vfpl(vecm_fevd, file = "vecm_fevd.pdf")
vecm_pub <- vpub(vecm_fit, n.ahead = 10, draws = 200)

print(csd_results)
print(unit_results)
print(coin_results$table)
print(pvar_diag)
print(head(pvar_irf$table))
print(head(pvar_fevd$table))
print(vecm_diag)
print(head(vecm_irf$table))
print(head(vecm_fevd))
