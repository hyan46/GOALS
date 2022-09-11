## Reference: 
Here is the code for the paper: "Jiang, Feiyu, Zifeng Zhao, and Xiaofeng Shao. "Modelling the COVID‐19 infection trajectory: A piecewise linear quantile trend model." Journal of the Royal Statistical Society: Series B (Statistical Methodology) (2021)." The paper can be found [here](https://rss.org.uk/RSS/media/File-library/Events/Discussion%20meetings/Jiang-et-al-_JRSSB_discussion-paper_16-Jul-21.pdf).

If you have any questions about the code, please feel free to contact Jiang, Feiyu with email jiangfy@fudan.edu.cn and Shao, Xiaofeng with email xshao@illinois.edu.

## Code

### Methods

- [Util_BIC.R](methods/Util_BIC.R)          contains general R functions associated with codes regarding multi-scanning M-GOALS in Section S5.2 in the supplement.
- [Util_pool.R](methods/Util_pool.R)         contains general R functions used for (M-)GOALS: it produces subsample testing statistics defined at each time point. 
- [Util_SN.R](methods/Util_SN.R)           contains general R functions used for generation of DGPs in simulations and implementation of (M-)GOALS: 1. global testing part and 2. local scanning part
- [SN_MultipleSlope_CriticalValue.csv](SN_MultipleSlope_CriticalValue.csv)           contains the critical values (90%,95%,99%) of the pivotal distribution for each trimming parameter pairs of M-GOALS

### Simulation
[simulation](simulation) folder contains GOALS for simulations part regarding DGP 1-3:
- [Util_single.R](simulation/Util_single.R)                        contains general R functions used to implement GOALS (single quantile). 
- [Simu_code_SingleTune.R](simulation/Simu_code_SingleTune.R)     is used to reproduce results in Table S1-S3 in the supplement, with fixed tuning parameters (epsilon=0.1, delta=0.02)，critical values (90%,95%,99%) for GOALS (single quantile) are tabulated in the code.
- [Simu_code_MultiBIC.R](simulation/Simu_code_MultiBIC.R)       is used to reproduce results in Table S4-S5 in the supplement with multiple tuning parameters combined by multi-bandwidth BIC.

### COVID Prediction
realdata folder contains GOALS for COVID-19 data analysis.
- [CDC_Forecast folder](CDC_Forecast)           contains predictions based on CDC website
- [covid_app.R](realdata/covid_app.R)                   is used to reproduce all the insample COVID-19 data analysis 
- [prediction.R](realdata/prediction.R)                  is used to reproduce all the out-of-sample forecasting results for U.S. daily new cases
- [Util_noncrossing.R](realdata/Util_noncrossing.R)	contains general R functions regarding quantile noncrossing constrain in Bondell et al. (2010).
- [owid-covid-data.csv](realdata/owid-covid-data.csv)           is the reported COVID related data download from "Our World in Data" Website.

### Reproducing Results

#### Figures
Code needed for reproducing application results in the main text.

- Figure 1-2:    [covid_app.R](realdata/covid_app.R)
- Figure 3:      [prediction.R](realdata/prediction.R)

#### Tables
Code needed for reproducing application results in the supplement

- Table S1-S3, Figure S1:     [Simu_code_SingleTune.R](simulation/Simu_code_SingleTune.R)
- Table S4-S5:                [Simu_code_MultiBIC.R](simulation/Simu_code_MultiBIC.R) 
- Table S6, Figure S2-S6:     [covid_app.R](realdata/covid_app.R)
- Table S7-S8, FigureS7-S8:   [prediction.R](realdata/prediction.R)
