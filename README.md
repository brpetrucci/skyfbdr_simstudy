README.Rmd
================
2025-04-01

## Skyline treeless FBD range simulation study

This project aims to validate the treeless version of the skyline FBD
range model (Stadler et al 2018, Warnock et al 2020). To do so, we
simulated stratigraphic ranges using the R package paleobuddy (do
Rosario Petrucci et al 2022), and applied FBD range to it using RevBayes
(Hoehna et al 2016).

### Simulation

The simulation code can be found in `simulation/sim.R`. That includes
both accuracy and coverage simulations.

For coverage simulations, we drew the model parameters (speciation rate,
extinction rate, fossil-sampling rate, and origin time) from prior
distributions (the same ones we use in the analysis code). This included
3 values for each rate, since this is a skyline model. Stages were
divided evenly between the present (0) and the origin time. 1000
simulations were ran, with rejection sampling to throw out replicates
where there were less than 5 or more than 20 species sampled (drawing
new values from the prior each time). They were then saved as fossil
counts (saved in `simulation/coverage/fossil_counts`, representing the
`k` parameter in FBD range), ranges (in `simulation/coverage/ranges`),
and stages (in `simulation/coverage/times`, note that 0 and the origin
time are not included) files. Also saved in `simulation/coverage` are
`seeds.RData`, keeping all the seeds randomly drawn to run the
simulations, `sim_list.RData`, keeping all the BD simulations that
generated the fossil records, and `true_vals.tsv`, listing the values
taken from the prior to run the simulations.

For accuracy simulations, we selected a number of axes of variation to
test the model’s strengths and limitations. `lambda`, `mu` and `psi`
might be constant (`"C"`), increase (`"LI3"` and `"LI2"`) or decrease
(`"LD3"` and `"LD2"`) linearly, or vary as a “triangle”, increasing
(`"UT"`) or decreasing (`"DT"`) in the middle interval. The
abbreviations in parentheses were used to identify how each rate varies,
as seen in the key below. Age also varied, and every rate scenario was
ran with `age = 10`, `20`, and `30`. Simulations were saved following
sampling model 3, i.e. with presence-absence matrices describing the
presence or absence of each species on each interval. Extant singletons
(extant species for which no fossils were sampled) were excluded. These
configurations form the bulk of the study, with 57 total simulation
scenarios, each ran with 100 reps. We also ran some simulations
diverging from these scenarios so as to explore different features of
the model:

- In scenario 58, `mu` and `psi` were set to have a constant rate, but
  will be analyzed in RevBayes with parameters for each stage. In this
  way, we will check whether the model can detect that there is no
  time-heterogeneity in rates.
- In scenario 59, `mu` was set to vary (`"UT"`), but will be analyzed in
  RevBayes assuming constant `mu`. This will allow us to check how the
  model reacts to parameter misspecification.
- In scenarios 60 and 61, all rates vary with 3 stages, but will be
  analyzed with parameters for 5 or 2 (respectively). We will then check
  how model estimates react to incorrect stage specification.
- In scenario 63, all rates vary and we include extant singletons in the
  saved ranges. By comparing with scenario 62 (which has the same
  parameters, but does not save extant singletons), we will check how
  extant singletons change our estimates.
- Finally, scenarios 64-67 will test the how the fossil data format
  impact analyses. Instead of saving fossil counts as presence-absence
  matrices, they are given to the model either as total fossil counts
  per interval (model 1), or not given at all (model 2, note that fossil
  counts were still saved). We also included scenarios to test whether
  uncertainty in models 1 and 2 impacts estimates, with ranges being
  saved with time bins as opposed to true values in that case.

For a review of simulation scenarios, see the key below. The names
column for the first 57 simulation is named as follows: a letter to
indicate which rate varies (`l`, `m` or `p`), the mode of variation
(e.g. `DT`), and underscore, and the age of the simulation (`10`, `20`
or `30`). These will also be the names of the directories within
`simulation/accuracy` where simulations are saved. The internal
structure of those directories is the same as in the case of the
coverage simulations. The key is saved under
`simulation/accuracy/key.tsv`, and within the `simulation/accuracy/refs`
directory, `Rev` files are saved to conveniently set parameters for the
RevBayes analyses (following the column names of `key` closely).

``` r
key
```

    ##                 names lambda  mu psi age model   unc extant_singles diff_bins
    ## 1             base_10      C   C   C  10     3 FALSE          FALSE         0
    ## 2             lLI3_10    LI3   C   C  10     3 FALSE          FALSE         0
    ## 3             lLI2_10    LI2   C   C  10     3 FALSE          FALSE         0
    ## 4             lLD3_10    LD3   C   C  10     3 FALSE          FALSE         0
    ## 5             lLD2_10    LD2   C   C  10     3 FALSE          FALSE         0
    ## 6              lUT_10     UT   C   C  10     3 FALSE          FALSE         0
    ## 7              lDT_10     DT   C   C  10     3 FALSE          FALSE         0
    ## 8             mLI3_10      C LI3   C  10     3 FALSE          FALSE         0
    ## 9             mLI2_10      C LI2   C  10     3 FALSE          FALSE         0
    ## 10            mLD3_10      C LD3   C  10     3 FALSE          FALSE         0
    ## 11            mLD2_10      C LD2   C  10     3 FALSE          FALSE         0
    ## 12             mUT_10      C  UT   C  10     3 FALSE          FALSE         0
    ## 13             mDT_10      C  DT   C  10     3 FALSE          FALSE         0
    ## 14            pLI3_10      C   C LI3  10     3 FALSE          FALSE         0
    ## 15            pLI2_10      C   C LI2  10     3 FALSE          FALSE         0
    ## 16            pLD3_10      C   C LD3  10     3 FALSE          FALSE         0
    ## 17            pLD2_10      C   C LD2  10     3 FALSE          FALSE         0
    ## 18             pUT_10      C   C  UT  10     3 FALSE          FALSE         0
    ## 19             pDT_10      C   C  DT  10     3 FALSE          FALSE         0
    ## 20            base_20      C   C   C  20     3 FALSE          FALSE         0
    ## 21            lLI3_20    LI3   C   C  20     3 FALSE          FALSE         0
    ## 22            lLI2_20    LI2   C   C  20     3 FALSE          FALSE         0
    ## 23            lLD3_20    LD3   C   C  20     3 FALSE          FALSE         0
    ## 24            lLD2_20    LD2   C   C  20     3 FALSE          FALSE         0
    ## 25             lUT_20     UT   C   C  20     3 FALSE          FALSE         0
    ## 26             lDT_20     DT   C   C  20     3 FALSE          FALSE         0
    ## 27            mLI3_20      C LI3   C  20     3 FALSE          FALSE         0
    ## 28            mLI2_20      C LI2   C  20     3 FALSE          FALSE         0
    ## 29            mLD3_20      C LD3   C  20     3 FALSE          FALSE         0
    ## 30            mLD2_20      C LD2   C  20     3 FALSE          FALSE         0
    ## 31             mUT_20      C  UT   C  20     3 FALSE          FALSE         0
    ## 32             mDT_20      C  DT   C  20     3 FALSE          FALSE         0
    ## 33            pLI3_20      C   C LI3  20     3 FALSE          FALSE         0
    ## 34            pLI2_20      C   C LI2  20     3 FALSE          FALSE         0
    ## 35            pLD3_20      C   C LD3  20     3 FALSE          FALSE         0
    ## 36            pLD2_20      C   C LD2  20     3 FALSE          FALSE         0
    ## 37             pUT_20      C   C  UT  20     3 FALSE          FALSE         0
    ## 38             pDT_20      C   C  DT  20     3 FALSE          FALSE         0
    ## 39            base_30      C   C   C  30     3 FALSE          FALSE         0
    ## 40            lLI3_30    LI3   C   C  30     3 FALSE          FALSE         0
    ## 41            lLI2_30    LI2   C   C  30     3 FALSE          FALSE         0
    ## 42            lLD3_30    LD3   C   C  30     3 FALSE          FALSE         0
    ## 43            lLD2_30    LD2   C   C  30     3 FALSE          FALSE         0
    ## 44             lUT_30     UT   C   C  30     3 FALSE          FALSE         0
    ## 45             lDT_30     DT   C   C  30     3 FALSE          FALSE         0
    ## 46            mLI3_30      C LI3   C  30     3 FALSE          FALSE         0
    ## 47            mLI2_30      C LI2   C  30     3 FALSE          FALSE         0
    ## 48            mLD3_30      C LD3   C  30     3 FALSE          FALSE         0
    ## 49            mLD2_30      C LD2   C  30     3 FALSE          FALSE         0
    ## 50             mUT_30      C  UT   C  30     3 FALSE          FALSE         0
    ## 51             mDT_30      C  DT   C  30     3 FALSE          FALSE         0
    ## 52            pLI3_30      C   C LI3  30     3 FALSE          FALSE         0
    ## 53            pLI2_30      C   C LI2  30     3 FALSE          FALSE         0
    ## 54            pLD3_30      C   C LD3  30     3 FALSE          FALSE         0
    ## 55            pLD2_30      C   C LD2  30     3 FALSE          FALSE         0
    ## 56             pUT_30      C   C  UT  30     3 FALSE          FALSE         0
    ## 57             pDT_30      C   C  DT  30     3 FALSE          FALSE         0
    ## 58               free     DT   C   C  30     3 FALSE          FALSE         0
    ## 59           not_free     DT  UT   C  30     3 FALSE          FALSE         0
    ## 60           stages_5     DT  UT  DT  30     3 FALSE          FALSE         5
    ## 61           stages_2     DT  UT  DT  30     3 FALSE          FALSE         2
    ## 62  no_extant_singles     DT  UT  DT  30     3 FALSE          FALSE         0
    ## 63 yes_extant_singles     DT  UT  DT  30     3 FALSE           TRUE         0
    ## 64          m1_no_unc     DT  UT  DT  30     1 FALSE          FALSE         0
    ## 65             m1_unc     DT  UT  DT  30     1  TRUE          FALSE         0
    ## 66          m2_no_unc     DT  UT  DT  30     2 FALSE          FALSE         0
    ## 67             m2_unc     DT  UT  DT  30     2  TRUE          FALSE         0

### Analysis

`analysis/coverage` contains
