Mpox delay distributions
================

# Setup

``` r
library("epidist")
library("dplyr")
```

# Data

We use data from<sup>\[[1](#ref-marziano2024incub_i)\]</sup>:

``` r
case_id <- seq_len(15)
incubation_period <- as.integer(c(5, 6, 7, 8, 8, 8, 9, 9, 10, 10, 11, 12, 13, 14, 14))
ip <- data.frame(case_id, incubation_period)
```

We first assume an offset of 4 days as
in<sup>\[[1](#ref-marziano2024incub_i)\]</sup>:

``` r
ip_offset <- ip |>
  mutate(incubation_period = incubation_period - 4L)
```

Next we prepare for fitting. We don’t know the primary event time so we
can only account for censoring here, not truncation

``` r
ip_offset_epidist <- ip_offset |>
  mutate(
    delay_daily = incubation_period,
    delay_lwr = incubation_period - 1,
    delay_upr = incubation_period + 1,
    censored = "interval"
  )
```

# Models

As in the paper, fit lognormal, weibull and gamma, here with and without
censoring adjustment.

``` r
lognormal <- naive_delay(data = ip_offset_epidist)
censored_lognormal <- censoring_adjusted_delay(data = ip_offset_epidist)
weibull <- naive_delay(
  brms::bf(delay_daily ~ 1, shape ~ 1), ip_offset_epidist, family = "weibull"
)
censored_weibull <- censoring_adjusted_delay(
  brms::bf(delay_lwr | cens(censored, delay_upr) ~ 1, shape ~ 1), ip_offset_epidist, family = "weibull"
)
gamma <- naive_delay(
  brms::bf(delay_daily ~ 1, shape ~ 1) ~ 1, data = ip_offset_epidist, family = "gamma"
)
censored_gamma <- censoring_adjusted_delay(
  brms::bf(delay_lwr | cens(censored, delay_upr) ~ 1, shape ~ 1), data = ip_offset_epidist, family = "gamma"
)
```

Model comparison

``` r
loo(
  lognormal,
  censored_lognormal,
  weibull,
  censored_weibull,
  gamma,
  censored_gamma
)
```

    ## Warning: Not all models have the same y variable. ('yhash' attributes do not
    ## match)

    ## Output of model 'lognormal':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -38.3 2.5
    ## p_loo         1.8 0.7
    ## looic        76.6 4.9
    ## ------
    ## MCSE of elpd_loo is 0.1.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 0.8]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Output of model 'censored_lognormal':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -27.7 2.4
    ## p_loo         1.6 0.5
    ## looic        55.3 4.8
    ## ------
    ## MCSE of elpd_loo is 0.0.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 0.7]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Output of model 'weibull':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -37.9 1.9
    ## p_loo         1.4 0.3
    ## looic        75.8 3.8
    ## ------
    ## MCSE of elpd_loo is 0.0.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 0.6]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Output of model 'censored_weibull':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -27.5 1.9
    ## p_loo         1.4 0.3
    ## looic        55.0 3.8
    ## ------
    ## MCSE of elpd_loo is 0.0.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 0.7]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Output of model 'gamma':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -38.0 2.2
    ## p_loo         1.6 0.4
    ## looic        76.0 4.4
    ## ------
    ## MCSE of elpd_loo is 0.0.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 0.7]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Output of model 'censored_gamma':
    ## 
    ## Computed from 4000 by 15 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo    -27.5 2.2
    ## p_loo         1.5 0.4
    ## looic        54.9 4.4
    ## ------
    ## MCSE of elpd_loo is 0.0.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 0.7]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## Model comparisons:
    ##                    elpd_diff se_diff
    ## censored_gamma       0.0       0.0  
    ## censored_weibull     0.0       0.6  
    ## censored_lognormal  -0.2       0.3  
    ## weibull            -10.4       0.6  
    ## gamma              -10.5       0.1  
    ## lognormal          -10.8       0.4

Estimates using censored gamma

``` r
draws <- exp(as_draws_df(censored_gamma)$Intercept)
signif(mean(draws), 2)
```

    ## [1] 9.7

``` r
signif(quantile(draws, c(0.025, 0.975)), 2)
```

    ##  2.5% 97.5% 
    ##   8.2  11.0

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-marziano2024incub_i" class="csl-entry">

1\. Marziano, V., Guzzetta, G., Longini, I., & Merler, S. (2024).
*Incubation period, serial interval, generation time and reproduction
number of mpox clade i*. <https://doi.org/10.1101/2024.05.10.24307157>

</div>

</div>
