# pavem

`pavem` provides beginner-friendly functions for Panel VAR and Panel VECM
workflows using existing CRAN packages.

## Important GitHub Upload Note

Upload the **contents of this `pavem_github` folder** to your GitHub repository
root. The repository root must contain:

```text
DESCRIPTION
NAMESPACE
R/
man/
inst/
README.md
LICENSE
```

Do not upload this package inside another folder such as `outputs/pavem`, or
`remotes::install_github("williamcantah/pavem")` may install an old or wrong
copy.

## Installation

Restart R, then run:

```r
remove.packages("pavem")
install.packages("remotes")
remotes::install_github("williamcantah/pavem", force = TRUE, upgrade = "never")
library(pavem)

packageVersion("pavem")
c("pipl", "pfpl", "vipl", "vfpl") %in% getNamespaceExports("pavem")
```

The version should be `0.2.0`, and all four export checks should return `TRUE`.

## Main Functions

| Function | Purpose |
|---|---|
| `prep()` | Prepare and sort panel data |
| `csdt()` | Cross-sectional dependence tests |
| `urts()` | Panel unit-root tests |
| `coint()` | Panel cointegration tests |
| `pvls()` | Panel VAR lag selection |
| `pvar()` | Estimate Panel VAR |
| `pvdi()` | Panel VAR diagnostics |
| `pvir()` | Panel VAR IRF table with confidence intervals |
| `pvfd()` | Panel VAR FEVD table |
| `pipl()` | Plot Panel VAR IRFs with confidence intervals |
| `pfpl()` | Plot Panel VAR FEVDs |
| `vecm()` | Estimate Panel VECM |
| `vedi()` | Panel VECM diagnostics |
| `veir()` | Panel VECM IRF table with confidence intervals |
| `vefd()` | Panel VECM FEVD table |
| `vipl()` | Plot Panel VECM IRFs with confidence intervals |
| `vfpl()` | Plot Panel VECM FEVDs |
| `vpub()` | Publication-ready Panel VECM results |

All exported function names have five characters or fewer.

## Example

```r
library(readxl)
library(pavem)

data <- read_excel("africa_pvar_pvecm_simulated_panel.xlsx", sheet = "panel_data")

data <- data[data$country_code %in% c("GHA", "NGA", "SLE", "SEN", "BEN"), ]

vars <- c("inflation", "gdp_growth", "interest_rate")

data <- prep(data, id = "country", time = "year", vars = vars)

csdt(data, id = "country", time = "year", vars = vars)
urts(data, id = "country", time = "year", vars = vars)
coint(data, id = "country", time = "year",
      y = "inflation", x = "interest_rate", method = "ped")

fit <- pvar(data, id = "country", time = "year", vars = vars)
pvdi(fit)
pvir(fit, n.ahead = 10, draws = 200)$table
pvfd(fit, n.ahead = 10)$table
pipl(fit, n.ahead = 10, draws = 200, file = "pvar_irf.pdf")
pfpl(fit, n.ahead = 10, file = "pvar_fevd.pdf")

vfit <- vecm(data, id = "country", time = "year",
             vars = vars, lags = 2, rank = 1)
vedi(vfit)
veir(vfit, n.ahead = 10, draws = 200)$table
vefd(vfit, n.ahead = 10)
vipl(vfit, n.ahead = 10, draws = 200, file = "vecm_irf.pdf")
vfpl(vfit, n.ahead = 10, file = "vecm_fevd.pdf")
vpub(vfit)
```
