# Spain Mortality Modelling and Longevity Risk

Original academic project by **Claudio Raffaele Salerno**, dated **2 July 2024**,
exploring Spanish male and female mortality through stochastic mortality models
and actuarial applications in R.

**This repository preserves the supplied original project.** The R code,
R Markdown, datasets, Word document and submitted PDF are unchanged. Only this
README and a Git ignore file have been added for presentation on GitHub.

## Project overview

The project compares six model families: **Lee–Carter, Renshaw–Haberman,
Age–Period–Cohort, Cairns–Blake–Dowd, PLAT and M7**. It includes parameter plots,
residual diagnostics, log-likelihood and BIC comparisons, mortality projections,
simulation, semiparametric bootstrap, life expectancy and period annuity
illustrations at a fixed 2% interest rate.

| Material | Original specification |
| --- | --- |
| Main R Markdown and submitted PDF | Fitted ages 0–100; estimation years 1950–2019 |
| Main projection horizon | 67 years, 2020–2086 |
| Main simulation configuration | 1,000 paths and 20 bootstrap replications |
| Separate R script | Alternative version using estimation years through 2021 and 100 paths |
| Supplied mortality snapshot | Spain, 1908–2021; file headers dated 18 September 2023 |
| Main R tools | StMoMo, demography, fanplot, ggplot2 and R Markdown |

The two source versions are retained as supplied and are not presented as
identical executions. BIC and log-likelihood comparisons concern in-sample
fit; the project does not establish holdout forecasting performance.

## Read the report

The main portfolio document is the [original submitted PDF](Progetto%20spagna%20FS.pdf).
It contains the methods, plots and results of the academic work, in Italian.
The [original Word document](Progetto%20-%20spagna.docx) is also included.

## Files

| File | Contents |
| --- | --- |
| `Progetto RMD spagna(1).Rmd` | Original R Markdown analysis corresponding to the submitted report |
| `mortalità spagnola codice.R` | Original standalone R script with the alternative estimation endpoint |
| `Functions.R` | Original model definition and actuarial helper functions |
| `Mx_1x1_spagna.txt` | Supplied age/year mortality rates |
| `Exposures_1x1_spagna.txt` | Supplied exposures to risk |
| `Deaths_1x1_spagna.txt` | Supplied deaths data |
| `Population_spagna.txt` | Supplied population data |
| `Progetto spagna FS.pdf` | Original submitted report |
| `Progetto - spagna.docx` | Original Word report |

## Execution context

The code is preserved in its original form, including filenames, package
requirements and workspace assumptions. It was not revised for cross-platform
execution or independently rerun for this repository. The supplied sources
refer to `Functions.r`, while the uploaded helper filename is `Functions.R`;
that case difference matters on case-sensitive filesystems. The original
package requirements are listed in the R Markdown source.

The repository is an archive of the academic work as submitted, rather than
a claim that the original calculations have undergone a new methodological
validation. Interpret results within the original project's assumptions.

## Data attribution

The standalone source identifies the **Human Mortality Database (HMD)** as
the data source. The supplied files retain their original bytes and headers.
Their original download date was not supplied.

HMD. Human Mortality Database. Max Planck Institute for Demographic Research
(Germany), University of California, Berkeley (USA), and French Institute for
Demographic Studies (France). Available at [mortality.org](https://www.mortality.org/).

See the [Spain country page](https://www.mortality.org/Country/Country?cntr=ESP)
and the [HMD user agreement](https://www.mortality.org/Data/UserAgreement)
for source documentation and data terms. No new license is assigned to the
supplied academic code or reports.
