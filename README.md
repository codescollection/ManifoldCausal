# ManifoldCausal: Causal Inference with Manifold‑Valued Treatments and Outcomes

[![R](https://img.shields.io/badge/R-%3E%3D%203.5-blue)](https://www.r-project.org)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Dependencies: R (>= 3.5), shapes, ggplot2, dplyr, tidyr, viridis.


## Overview

**ManifoldCausal** implements the manifold kernel instrumental variable (MKIV) estimator and the geodesic local average treatment effect (G‑LATE) proposed in *Shape Causal Inference with Manifold Endogeneity:
A Nonparametric IV Approach*. It provides an R framework for causal inference when treatments and/or outcomes reside on non‑Euclidean shape manifolds (e.g., road networks, anatomical contours). The package avoids ad‑hoc Euclidean approximations by working directly with Procrustes shape distances and reproducing kernel Hilbert spaces on Riemannian manifolds. *(under construction...)*

## Key Features

- **Shape alignment and distance matrix computation** – Joint Procrustes registration and full Procrustes distance matrices from landmark configurations.
- **Manifold kernel IV (MKIV) estimation** – Two‑stage kernel IV that respects the intrinsic geometry of shape spaces.
- **Automatic regularization selection** – Leave‑one‑out cross‑validation for the ridge penalty.
- **Structural function visualisation** – Plot \(\hat{g}(d)\) along shape principal components with observed data.
- **Geodesic local average treatment effect (G‑LATE)** – Weighted Fréchet mean estimation for shape‑valued outcomes and a binary instrument.
- **Built‑in example data** – Road‑network shapes from four Chinese cities (Beijing, Chongqing, Hong Kong, Macau) aligned to 150 landmarks.
- **Reproducible vignette** – Step‑by‑step reproduction of the empirical analysis linking road morphology to GDP per capita.

## Installation

The development version of ManifoldCausal can be installed directly from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("codescollection/ManifoldCausal", build_vignettes = TRUE)




Quick Start
library(ManifoldCausal)

# Load built-in road shape data (k = 150 landmarks, 4 cities)
data(road_shape_array)

# Step 1: Align shapes and compute Procrustes distance matrix
align <- shape_alignment(road_shape_array)

# Step 2: Define outcome Y (log GDP per capita) and instrument Z
Y <- c(12.38, 11.57, 12.92, 13.19)   # Beijing, Chongqing, HK, Macau
Z <- c(0.8, 1.2, -0.6, -1.0)         # composite geographic‑historical IV

# Step 3: Build Gaussian kernel matrices
sigma_D <- median(align$dist[upper.tri(align$dist)])
K_D <- gauss_kernel(align$dist, sigma_D)
K_Z <- gauss_kernel(as.matrix(dist(Z)), median(dist(Z)))

# Step 4: Fit the MKIV estimator
fit <- mkiv_fit(Y, K_D, K_Z, lambda = 0.1, rho1 = 0.1)

# Step 5: Cross‑validate lambda
best_lambda <- cv_mkiv(Y, K_D, K_Z)
fit_opt <- mkiv_fit(Y, K_D, K_Z, lambda = best_lambda, rho1 = 0.1)

# Step 6: Plot the structural function
plot_structural(align$gpa, fit_opt$gamma, Y, road_shape_array,
                sigma_D = sigma_D,
                cities = c("Beijing","Chongqing","Hong Kong","Macau"))

```



| Function | Description |
|----------|-------------|
| `shape_alignment(x)` | Procrustes alignment of a `k x 2 x n` array, returns aligned shapes, tangent scores, and distance matrix. |
| `gauss_kernel(D, sigma)` | Gaussian kernel from a distance matrix. |
| `mkiv_fit(Y, K_D, K_Z, lambda, rho1)` | Two‑stage manifold kernel IV estimation. Returns fitted values, structural coefficients `gamma`. |
| `cv_mkiv(Y, K_D, K_Z, lambda_grid, rho1)` | Leave‑one‑out cross‑validation for the regularization parameter. |
| `plot_structural(gpa, gamma, Y, shape_array, sigma_D, cities)` | Plots the estimated dose–response curve along PC1. |
| `glate(Y, T, Z)` | Estimates the geodesic local average treatment effect for shape‑valued outcomes and binary instrument. |




## Contributing

Contributions are welcome! Please submit bug reports, feature requests, or pull requests via the GitHub issue tracker.



## Citation

If you use ManifoldCausal in your research, please cite the accompanying paper:

Shape Causal Inference with Manifold Endogeneity: A Nonparametric IV Approach. (submitted)

please contact haoshiming@ynu.edu.cn if you have any questions concerning the use or methedologies of "ManifoldCausal".



## License

ManifoldCausal is released under the GNU General Public License v3.0.
