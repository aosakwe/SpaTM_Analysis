# SpaTM Analysis Scripts

This repository contains the analysis code used to produce the results in our recent manuscript titled [_SpaTM: Topic Models for Inferring Spatially Informed Transcriptional Programs_](https://www.biorxiv.org/content/10.1101/2025.01.24.634726v1). The code used for each figure are stored in Quarto notebooks.

The SpaTM R package is stored in a separate repository and [can be found here](https://github.com/li-lab-mcgill/SpaTM/tree/main).

**IMPORTANT** The scripts were run using version 0.2 of ```SpaTM```. Thus, to run these scripts, you must use the corresponding version stored within a separate branch [(archive/v0.2)](https://github.com/li-lab-mcgill/SpaTM/tree/archive/v0.2). To install the specified branch, you can run the following:
```
devtools::install_github("li-lab-mcgill/SpaTM", ref = "archive/v0.2")
```

## Datasets used
The DLPFC 10X Visium dataset was acquired from the SpatialLIBD project. It can be directly loaded into R by using the ```spatialLIBD``` package as follows:
```
## Load the package
library("spatialLIBD")

## Download the spot-level data
spe <- fetch_data(type = "spe")
```
More information on the above data can be found through [this portal](https://research.libd.org/spatialLIBD/) created by the authors of the original study.

The Major Depressive Disorder (MDD) snRNA-seq dataset was acquired through [a separate study](https://www.nature.com/articles/s41467-023-38530-5).

## SpaTM Workflow
![SpaTM Workflow](https://github.com/aosakwe/SpaTM_Analysis/blob/main/SpaTM.png)
