# landmark-test
<!-- Authors: [Vera Weisbecker](v.weisbecker@uq.edu.au), [Thomas Guillerme](mailto:guillert@tcd.ie)... -->

This repository contains all the code and data used in the manuscript [Link to final published pdf will be here]().

<!-- To cite the paper:  -->
<!-- > Thomas Guillerme \& Martin Brazeau. 2018. Influence of different modes of morphological character correlation on phylogenetic tree inference -->

<!-- To cite this repo:  -->
<!-- > Thomas Guillerme \& Martin Brazeau. 2018. Influence of different modes of morphological character correlation on phylogenetic tree inference -->

# Data

<!-- All the data used in the manuscript is archive on [Figshare](https://figshare.com/s/7a8fde8eaa39a3d3cf56). -->



# Analyses

The tables and figures present in the manuscript are all reproducible through the following scripts:

## 01 - Data preparation

This script contains all the procedure to modify and prepare the data in a format usable in the analysis below,
This script is available [here in Rmd](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/01-Data_preparation.Rmd) or [here in html](https://rawgit.com/TGuillerme/landmark-test/master/Analysis/01-Data_preparation.html).

## 02 - Landmark region difference

This vignette describes in details the test used to compare different regions of landmarks in geometric morphometrics.
This example uses the `plethodon` dataset from the `geomorph` package as an illustration.
This script is available [here in Rmd](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/02-Landmark_region_difference.Rmd) or [here in html](https://rawgit.com/TGuillerme/landmark-test/master/Analysis/02-Landmark_region_difference.html).

## 03 - Landmark test analysis and 04 - landmark test analysis rarefied

These scripts runs the full landmark test analysis on each tested partition for each dataset (cranium and mandible) with and without rarefaction.
These script take some time to run (~30 min each) and will automatically save data in `../Data/Results/`.
If the results are already computed, it is possible to run the script faster by commenting out the lines containing:

```{r}
results <- pipeline.test(species, dataset, "../Data/Processed/", [...])
save(results, file = paste0("../Data/Results/", species, [...]))
```

And un-commenting the lines containing:

```{r}
load(file = paste0("../Data/Results/", species, "_", [...]))
```

The test run for either the 100\% confidence interval or the 95\% confidence intervals.
To switch between both, you can search and replace `CI = 0.95` and `CI95` by respectively `CI = 1` and `CI100` or the other way around.

This script is available [here in Rmd](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/03-Landmark_test_analysis.Rmd) ([here for the rarefied](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/04-Landmark_test_analysis_rarefied.Rmd)).

## 05 - Results summary

This script produces the figures @@@ and @@@ in the manuscript and the tables @@@ in the supplementary materials.
This script is available [here in Rmd](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/05-Results_summary.Rmd) or [here in html](https://rawgit.com/TGuillerme/landmark-test/master/Analysis/05-Results_summary.html).

## 06 - 


## 07 - 


## 08 - 




# Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("2018-06-15")
```

The analysis were run on two different machines.
For reproducibility purposes, the output of `devtools::session_info()` used to perform the analyses in the publication for the analysis `03-Landmark_test_analysis.Rmd`, `04-Landmark_test_analysis_rarefied.Rmd` and `05-Results_summary.Rmd` is available [here](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/Session_info-2018-06-15_machine1.txt).
The `devtools::session_info()` used to perform the analyses in the publication for the analysis `06-Main_Analyses.Rmd`, `07-heatplots.Rmd` and `08-Figures.Rmd` is available [here](https://github.com/TGuillerme/landmark-test/blob/master/Analysis/Session_info-2018-06-15_machine2.txt).
The other analysis where run on both machines.
