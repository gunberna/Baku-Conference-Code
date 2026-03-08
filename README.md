# Baku-Conference-Code
The code for the ITTA 2026 conference paper

## Requirements
Run the following at the start of your R session to install any missing packages and load them:
```r
required_packages <- c("MASS", "mclust", "e1071", "ggplot2", "gridExtra", 
                        "clusterGeneration", "ppclust", "cluster", "factoextra", 
                        "fcvalid", "dplyr", "pastecs", "patchwork")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)
```

## Structure
1. Kwon Mahalanobis Function (line 34)
2. Outlier Experiment (line 90)
   - Visualization (line 127)
   - Comparison between indices (line 189)
3. Monte Carlo Stress Tests
   - 30 Trial Stress Test (line 266)
   - ⚠️ 50 Trial Index Comparison (line 314) — takes a long time to run
   - Stress Test (line 427)
   - Visualization (line 468)
4. Kwon Calculations
   - Manual Calculation (line 533)
   - Mahalanobis Calculation (line 573)
   - Performance Comparison (line 624)
5. Mahalanobis vs Euclidean Visualization (line 665)
