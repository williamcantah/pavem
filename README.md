# pavem

`pavem` provides beginner-friendly functions for Panel VAR and Panel VECM
workflows using existing CRAN packages.

The package is designed for teaching, applied macroeconomics, central banking,
and policy training settings where users need one-command wrappers for common
panel time-series tasks.

## Installation

After uploading this folder to GitHub under `williamcantah/pavem`, users can
install it with:

```r
install.packages("remotes")
remotes::install_github("williamcantah/pavem")
```

If R says a function such as `pipl()` or `vipl()` cannot be found, restart R and
reinstall the latest version:

```r
remove.packages("pavem")
remotes::install_github("williamcantah/pavem", force = TRUE)
library(pavem)

"pipl" %in% getNamespaceExports("pavem")
"vipl" %in% getNamespaceExports("pavem")
```

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
| `pvir()` | Panel VAR IRFs with confidence intervals |
| `pvfd()` | Panel VAR FEVD tables |
| `pipl()` | Plot Panel VAR IRFs with confidence intervals |
| `pfpl()` | Plot Panel VAR FEVDs |
| `vecm()` | Estimate Panel VECM |
| `vels()` | Panel VECM lag selection |
| `vedi()` | Panel VECM diagnostics |
| `veir()` | Panel VECM IRFs with confidence intervals |
| `vefd()` | Panel VECM FEVD tables |
| `vipl()` | Plot Panel VECM IRFs with confidence intervals |
| `vfpl()` | Plot Panel VECM FEVDs |
| `vpub()` | Publication-ready Panel VECM results |

All function names are five characters or fewer.

## Quick Example

```r
library(readxl)
library(pavem)

data <- read_excel("africa_pvar_pvecm_simulated_panel.xlsx", sheet = "panel_data")

data <- prep(
  data = data,
  id = "country",
  time = "year",
  vars = c("inflation", "gdp_growth", "interest_rate", "m2_growth")
)

csdt(data, id = "country", time = "year",
     vars = c("inflation", "gdp_growth", "interest_rate"))

urts(data, id = "country", time = "year",
     vars = c("inflation", "gdp_growth", "interest_rate"))

coint(data, id = "country", time = "year",
      y = "inflation", x = "interest_rate", method = "ped")

vars <- c("inflation", "gdp_growth", "interest_rate")

lag_pvar <- pvls(data, id = "country", time = "year", vars = vars, lags = 1:3)
lag_pvar$table
lag_pvar$selected

fit <- pvar(data, id = "country", time = "year",
            vars = vars, lags = "auto")

pvdi(fit)
pvir(fit, n.ahead = 10, draws = 200)$table
pvfd(fit, n.ahead = 10)$table
pipl(fit, n.ahead = 10, draws = 200, file = "pvar_irf.pdf")
pfpl(fit, n.ahead = 10, file = "pvar_fevd.pdf")

lag_vecm <- vels(data, id = "country", time = "year", vars = vars, lags = 1:4)
lag_vecm$table
lag_vecm$selected

vfit <- vecm(data, id = "country", time = "year",
             vars = vars, lags = "auto", rank = 1)

vedi(vfit)
veir(vfit, n.ahead = 10, draws = 200)$table
vefd(vfit, n.ahead = 10)
vipl(vfit, n.ahead = 10, draws = 200, file = "vecm_irf.pdf")
vfpl(vfit, n.ahead = 10, file = "vecm_fevd.pdf")
vpub(vfit)
```

## Notes

- `pvar()` uses `panelvar::pvargmm()`.
- `vecm()` uses `pvars::pvarx.VEC()`.
- `coint(method = "ped")` uses `pco::pedroni99()`.
- `coint(method = "jo")` uses `pvars::pcoint.JO()`.
- Panel VECM IRF confidence intervals are obtained by resampling countries,
  re-estimating the Panel VECM, and taking bootstrap quantiles.
