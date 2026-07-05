\name{pavem}
\alias{pavem}
\alias{prep}
\alias{csdt}
\alias{urts}
\alias{coint}
\alias{pvls}
\alias{pvar}
\alias{pipl}
\alias{pvdi}
\alias{pvir}
\alias{pvfd}
\alias{pfpl}
\alias{vecm}
\alias{vedi}
\alias{veir}
\alias{vefd}
\alias{vels}
\alias{vipl}
\alias{vfpl}
\alias{vpub}
\title{Simple Panel VAR and Panel VECM Workflows}
\description{
Beginner-friendly wrappers for Panel VAR and Panel VECM analysis using existing
CRAN packages. The functions return plain tables and model objects for
cross-sectional dependence tests, panel unit-root tests, panel cointegration
tests, PVAR estimation, Panel VECM estimation, diagnostics, IRFs, FEVD, and
publication-ready summaries.
}
\usage{
prep(data, id, time, vars = NULL)
csdt(data, id, time, vars, test = "cd", alpha = 0.05)
urts(data, id, time, vars, tests = c("ips", "madwu"),
     exo = "intercept", lags = 1, alpha = 0.05)
coint(data, id, time, y = NULL, x = NULL, vars = NULL,
      method = c("ped", "jo"), lags = 1, type = "Case3", alpha = 0.05)
pvls(data, id, time, vars, lags = 1:3, transformation = "fod",
     steps = "onestep", system = TRUE, max_instr = 4,
     min_instr = 2L, collapse = TRUE, prefer = "MBIC")
pvar(data, id, time, vars, lags = 1, transformation = "fod",
     steps = "onestep", system = TRUE, max_instr = 4,
     min_instr = 2L, collapse = TRUE, prefer = "MBIC")
pvdi(obj)
pvir(obj, n.ahead = 10, draws = 200, ci = 0.95, cores = 1)
pvfd(obj, n.ahead = 10)
pipl(x, n.ahead = 10, draws = 200, ci = 0.95, file = NULL,
     col = "steelblue", shade = grDevices::rgb(0.2, 0.4, 0.8, 0.20))
pfpl(x, n.ahead = 10, file = NULL,
     cols = c("steelblue", "goldenrod", "darkseagreen",
              "firebrick", "mediumpurple", "gray50"))
vels(data, id, time, vars, lags = 1:4, type = "const", prefer = "BIC")
vecm(data, id, time, vars, lags = 2, rank = 1,
     type = "Case3", order = vars, prefer = "BIC")
vedi(obj)
veir(obj, n.ahead = 10, draws = 200, seed = 123)
vefd(obj, n.ahead = 10)
vipl(x, n.ahead = 10, draws = 200, seed = 123, file = NULL,
     col = "firebrick", shade = grDevices::rgb(0.7, 0.3, 0.1, 0.20))
vfpl(x, n.ahead = 10, file = NULL,
     cols = c("steelblue", "goldenrod", "darkseagreen",
              "firebrick", "mediumpurple", "gray50"))
vpub(obj, n.ahead = 10, draws = 200)
}
\arguments{
\item{data}{A panel data frame.}
\item{id}{Name of the cross-section identifier column.}
\item{time}{Name of the time identifier column.}
\item{vars}{Character vector of variable names.}
\item{test}{Cross-sectional dependence test type passed to \code{plm::pcdtest}.}
\item{alpha}{Significance level used for decisions.}
\item{tests}{Panel unit-root tests passed to \code{plm::purtest}.}
\item{exo}{Deterministic component for panel unit-root tests.}
\item{lags}{Lag order or lag orders.}
\item{prefer}{Preferred lag-selection criterion. For PVAR, use \code{"MBIC"}, \code{"MAIC"}, or \code{"MQIC"}. For Panel VECM, use \code{"BIC"}, \code{"AIC"}, \code{"HQ"}, or \code{"FPE"}.}
\item{y}{Dependent variable for Pedroni cointegration test.}
\item{x}{Regressor variable for Pedroni cointegration test, or a fitted model / IRF / FEVD result object for plotting functions.}
\item{method}{Cointegration method: \code{"ped"} or \code{"jo"}.}
\item{type}{Deterministic case for \code{pvars} functions.}
\item{transformation}{Transformation passed to \code{panelvar::pvargmm}.}
\item{steps}{GMM step option passed to \code{panelvar::pvargmm}.}
\item{system}{Whether to use system instruments in \code{panelvar::pvargmm}.}
\item{max_instr}{Maximum dependent-variable instruments.}
\item{min_instr}{Minimum dependent-variable instruments.}
\item{collapse}{Whether to collapse instruments.}
\item{obj}{Object returned by \code{pvar()} or \code{vecm()}.}
\item{n.ahead}{IRF/FEVD horizon.}
\item{draws}{Bootstrap draws.}
\item{ci}{Confidence interval level.}
\item{cores}{Number of cores for \code{panelvar::bootstrap_irf}.}
\item{rank}{Panel VECM cointegration rank.}
\item{order}{Cholesky ordering for Panel VECM identification.}
\item{seed}{Random seed for Panel VECM IRF bootstrap.}
\item{file}{Optional PDF path for plot output.}
\item{col}{Line color for IRF plots.}
\item{shade}{Confidence-band shading color for IRF plots.}
\item{cols}{Bar colors for FEVD plots.}
}
\value{
Data frames, lists, or fitted model objects depending on the function. IRF and
FEVD functions include table outputs suitable for export or publication.
}
\examples{
\dontrun{
library(readxl)
library(pavem)

d <- read_excel("africa_pvar_pvecm_simulated_panel.xlsx", sheet = "panel_data")
d <- prep(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))

csdt(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
urts(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
coint(d, "country", "year", y = "inflation", x = "interest_rate")

pf <- pvar(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
pvls(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
pfa <- pvar(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"), lags = "auto")
pvir(pf)$table
pvfd(pf)$table
pipl(pf, file = "pvar_irf.pdf")
pfpl(pf, file = "pvar_fevd.pdf")

vf <- vecm(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
vels(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"))
vfa <- vecm(d, "country", "year", c("inflation", "gdp_growth", "interest_rate"), lags = "auto")
vpub(vf)
vipl(vf, file = "vecm_irf.pdf")
vfpl(vf, file = "vecm_fevd.pdf")
}
}
